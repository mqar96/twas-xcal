import argparse, glob, json, os, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, t as tdist
from xcal.calibration import Panel
from xcal.core import ld_scores, estimate_c, satterthwaite_df

ap = argparse.ArgumentParser(description="End-to-end benchmark: c and intercept a estimated by LD-score regression, leave-one-chromosome-out")
ap.add_argument("genodir"); ap.add_argument("weightdir"); ap.add_argument("out")
ap.add_argument("--genomes", type=int, default=300)
ap.add_argument("--chromosomes", type=int, default=10)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--se-genomes", type=int, default=5, help="genomes per setting used to obtain block-jackknife SEs of (a, c)")
a = ap.parse_args()
POPS = ("YRI", "CEU", "CHB")
CELLS = [(1.0, 0.0, 1.0), (1.0, 0.1, 1.0), (1.0, 0.4, 1.0), (0.01, 0.1, 1.0), (0.01, 0.4, 1.0), (1.0, 0.1, 1.1), (1.0, 0.4, 1.1)]
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
PART = a.out + ".parts"
os.makedirs(PART, exist_ok=True)

files = sorted(glob.glob(f"{a.weightdir}/weights_*.npz"))
rows, crow = [], []
t0 = time.time()
for gidx, g in enumerate(POPS):
    if all(os.path.exists(f"{PART}/{g}_{ci}.json") for ci in range(len(CELLS))):
        for ci in range(len(CELLS)):
            pr = json.load(open(f"{PART}/{g}_{ci}.json"))
            rows.extend(pr["rows"]); crow.extend(pr["crow"])
        continue
    regs = []
    for wf in files:
        rid = int(os.path.basename(wf)[8:12])
        d = np.load(f"{a.genodir}/region_{rid:04d}.npz")
        if int(d["n_snps"]) < 20:
            continue
        P = Panel(d[f"{g}_ldref"])
        Rg = d[f"{g}_R_gwas"].astype(np.float64)
        sG = d[f"{g}_sd_gwas"]
        lam, U = np.linalg.eigh(Rg)
        genes = []
        for r in json.loads(str(np.load(wf)["recs"])):
            if r["n_w"] == 0:
                continue
            w = np.zeros(int(d["n_snps"])); w[r["idx"]] = r["w"]
            v = w * sG
            m = P.moments(v)
            if m["L"] <= 1e-12:
                continue
            Lb, Qb = P.jackknife(v, K=20, seed=rid)
            av = Rg @ v
            genes.append(dict(v=v, L=m["L"], Q=m["Q"], Lb=Lb, Qb=Qb, C0=float(v @ av), C1=float(av @ av)))
        regs.append(dict(rid=rid, chrom=len(regs) * a.chromosomes // 500, l=ld_scores(d[f"{g}_ldref"], d["pos"], window_bp=2_000_000),
                         Rg=Rg, F=U * np.sqrt(np.clip(lam, 0, None)), genes=genes))
    nreg = len(regs)
    for x in regs:
        x["chrom"] = min(x["chrom"], a.chromosomes - 1)
    ldall = np.concatenate([x["l"] for x in regs])
    chrom_snp = np.concatenate([np.full(len(x["l"]), x["chrom"]) for x in regs])
    for ci, (pi, c, a_true) in enumerate(CELLS):
        pf = f"{PART}/{g}_{ci}.json"
        if os.path.exists(pf):
            pr = json.load(open(pf)); rows.extend(pr["rows"]); crow.extend(pr["crow"])
            continue
        rng = np.random.default_rng([808, gidx, ci])
        r0, c0_ = len(rows), len(crow)
        zall = np.zeros((a.genomes, len(ldall)), np.float32)
        nums, meta = [], []
        o = 0
        for ri, x in enumerate(regs):
            p = len(x["l"])
            eps = np.sqrt(a_true) * (x["F"] @ rng.standard_normal((p, a.genomes))).T
            if c == 0:
                sig = 0.0
            elif pi == 1.0:
                sig = rng.normal(0, np.sqrt(c), (a.genomes, p)) @ x["Rg"]
            else:
                mask = rng.random((a.genomes, p)) < pi
                rr, cc = np.nonzero(mask)
                sig = sparse.csr_matrix((rng.normal(0, np.sqrt(c / pi), len(rr)), (rr, cc)), shape=(a.genomes, p)) @ x["Rg"]
            z = eps + sig
            zall[:, o:o + p] = z
            o += p
            for gg in x["genes"]:
                nums.append(z @ gg["v"]); meta.append((ri, x["chrom"], gg))
        nums = np.array(nums)
        chrom_gene = np.array([m[1] for m in meta]); reg_gene = np.array([m[0] for m in meta])
        L = np.array([m[2]["L"] for m in meta]); Q = np.array([m[2]["Q"] for m in meta])
        Lb = np.array([m[2]["Lb"] for m in meta]); Qb = np.array([m[2]["Qb"] for m in meta])
        Vt = np.array([a_true * m[2]["C0"] + c * m[2]["C1"] for m in meta])
        fits = {}
        for mode, fix in (("free", None), ("fixed", 1.0)):
            ch = np.zeros((a.genomes, a.chromosomes)); ah = np.zeros((a.genomes, a.chromosomes))
            for k in range(a.chromosomes):
                keep = chrom_snp != k
                for gi in range(a.genomes):
                    f = estimate_c(zall[gi, keep], ldall[keep], n_blocks=0, fix_intercept=fix)
                    ch[gi, k], ah[gi, k] = f["c"], f["intercept"]
            se = [estimate_c(zall[gi, chrom_snp != 0], ldall[chrom_snp != 0], n_blocks=200, fix_intercept=fix) for gi in range(a.se_genomes)]
            fits[mode] = (ch, ah, float(np.mean([s["se_c"] for s in se])), float(np.mean([s["se_intercept"] for s in se])))
            crow.append(dict(gwas_pop=g, pi=pi, c_true=c, a_true=a_true, intercept_mode=mode, genomes=a.genomes, mean_chat=ch.mean(),
                             sd_chat=ch.std(), mean_ahat=ah.mean(), sd_ahat=ah.std(), mean_jackknife_se_c=fits[mode][2],
                             mean_jackknife_se_a=fits[mode][3]))
        variants = [("LOCO free intercept: a_hat, c_hat, SEs propagated", "free", True, True),
                    ("LOCO free intercept: a_hat, c_hat, SEs not propagated", "free", True, False),
                    ("LOCO fixed intercept a=1: c_hat, SE propagated", "fixed", False, True),
                    ("true a and c", None, None, None)]
        for label, mode, use_a, prop in variants:
            if mode is None:
                cc = np.full((len(meta), a.genomes), c); aa = np.full_like(cc, a_true); se_c = se_a = 0.0
            else:
                ch, ah, se_c, se_a = fits[mode]
                cc = ch[:, chrom_gene].T
                aa = (np.maximum(ah[:, chrom_gene].T, 1.0) if use_a else np.ones_like(cc))
                if not prop:
                    se_c = se_a = 0.0
                if not use_a:
                    se_a = 0.0
            D = aa * L[:, None] + cc * Q[:, None]
            Db = aa[:, None, :] * Lb[:, :, None] + cc[:, None, :] * Qb[:, :, None]
            df = satterthwaite_df(np.moveaxis(Db, 1, -1), D, extra_var=(Q[:, None] * se_c) ** 2 + (L[:, None] * se_a) ** 2)
            stat = np.abs(nums) / np.sqrt(D)
            r = dict(gwas_pop=g, pi=pi, c_true=c, a_true=a_true, method=label, n_tests=int(stat.size), genomes=a.genomes,
                     se_c_used=se_c, se_a_used=se_a)
            bw_g = rng.multinomial(a.genomes, np.full(a.genomes, 1 / a.genomes), size=a.boot)
            bw_r = rng.multinomial(nreg, np.full(nreg, 1 / nreg), size=a.boot)
            for al in ALPHAS:
                thr = tdist.isf(al / 2, df)
                rej = stat > thr
                C = np.zeros((a.genomes, nreg)); N = np.zeros((a.genomes, nreg))
                np.add.at(C.T, reg_gene, rej.astype(float))
                np.add.at(N.T, reg_gene, np.ones_like(rej, float))
                bc = ((bw_g @ C) * bw_r).sum(1) / ((bw_g @ N) * bw_r).sum(1)
                r[f"count_{al:g}"] = int(rej.sum()); r[f"expected_{al:g}"] = stat.size * al
                r[f"ratio_{al:g}"] = rej.mean() / al
                r[f"lo_{al:g}"] = np.quantile(bc, .025) / al; r[f"hi_{al:g}"] = np.quantile(bc, .975) / al
                if pi == 1.0:
                    r[f"conditional_ratio_{al:g}"] = float(np.mean(2 * norm.sf(thr * np.sqrt(D / Vt[:, None])))) / al
            rows.append(r)
        json.dump(dict(rows=rows[r0:], crow=crow[c0_:]), open(pf + ".tmp", "w"), default=float)
        os.replace(pf + ".tmp", pf)
        print(g, pi, c, a_true, f"{time.time() - t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
pd.DataFrame(crow).to_csv(a.out.replace(".csv", "_c_estimates.csv"), index=False)
