import argparse, glob, json, os, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, t as tdist
from xcal.calibration import Panel
from xcal.core import ld_scores, estimate_c, satterthwaite_df

ap = argparse.ArgumentParser(description="End-to-end benchmark with c estimated by LD-score regression from each simulated genome")
ap.add_argument("genodir"); ap.add_argument("weightdir"); ap.add_argument("out")
ap.add_argument("--genomes", type=int, default=100)
a = ap.parse_args()
POPS = ("YRI", "CEU", "CHB")
CELLS = [(1.0, 0.0), (1.0, 0.1), (1.0, 0.4), (0.01, 0.1), (0.01, 0.4)]
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
rng = np.random.default_rng(808)

files = sorted(glob.glob(f"{a.weightdir}/weights_*.npz"))
rows, crow = [], []
t0 = time.time()
for g in POPS:
    regs = []
    for wf in files:
        rid = int(os.path.basename(wf)[8:12])
        d = np.load(f"{a.genodir}/region_{rid:04d}.npz")
        if int(d["n_snps"]) < 20:
            continue
        P = Panel(d[f"{g}_ldref"])
        l = ld_scores(d[f"{g}_ldref"], d["pos"], window_bp=2_000_000)
        Rg = d[f"{g}_R_gwas"].astype(np.float64)
        lam, U = np.linalg.eigh(Rg)
        F = U * np.sqrt(np.clip(lam, 0, None))
        genes = []
        for r in json.loads(str(np.load(wf)["recs"])):
            if r["n_w"] == 0:
                continue
            w = np.zeros(int(d["n_snps"])); w[r["idx"]] = r["w"]
            m = P.moments(w)
            if m["L"] <= 1e-12:
                continue
            Lb, Qb = P.jackknife(w, K=20, seed=rid)
            v = P.v(w); av = Rg @ v
            genes.append(dict(gene=r["gene"], v=v, L=m["L"], Q=m["Q"], Lb=Lb, Qb=Qb, C0=float(v @ av), C1=float(av @ av)))
        regs.append(dict(l=l, Rg=Rg, F=F, genes=genes))
    ldall = np.concatenate([x["l"] for x in regs])
    for pi, c in CELLS:
        zall = np.zeros((a.genomes, len(ldall)), np.float32)
        nums, meta = [], []
        o = 0
        for x in regs:
            p = len(x["l"])
            eps = (x["F"] @ rng.standard_normal((p, a.genomes))).T
            if c == 0:
                sig = 0.0
            elif pi == 1.0:
                sig = rng.normal(0, np.sqrt(c), (a.genomes, p)) @ x["Rg"]
            else:
                mask = rng.random((a.genomes, p)) < pi
                rr, cc = np.nonzero(mask)
                U_ = sparse.csr_matrix((rng.normal(0, np.sqrt(c / pi), len(rr)), (rr, cc)), shape=(a.genomes, p))
                sig = U_ @ x["Rg"]
            z = eps + sig
            zall[:, o:o + p] = z
            o += p
            for gg in x["genes"]:
                nums.append(z @ gg["v"]); meta.append(gg)
        nums = np.array(nums)
        fits = [estimate_c(zall[k], ldall, n_blocks=0) for k in range(a.genomes)]
        fits1 = [estimate_c(zall[k], ldall, n_blocks=0, fix_intercept=1.0) for k in range(a.genomes)]
        se_fit = estimate_c(zall[0], ldall, n_blocks=200)
        chat = np.array([f["c"] for f in fits]); chat1 = np.array([f["c"] for f in fits1])
        ahat = np.array([f["intercept"] for f in fits])
        crow.append(dict(gwas_pop=g, pi=pi, c_true=c, n_snps=len(ldall), genomes=a.genomes, mean_chat=chat.mean(), sd_chat=chat.std(),
                         mean_chat_fixed_intercept=chat1.mean(), sd_chat_fixed_intercept=chat1.std(),
                         mean_intercept=ahat.mean(), jackknife_se_example=se_fit["se_c"]))
        L = np.array([m["L"] for m in meta]); Q = np.array([m["Q"] for m in meta])
        Lb = np.array([m["Lb"] for m in meta]); Qb = np.array([m["Qb"] for m in meta])
        Vt = np.array([m["C0"] + c * m["C1"] for m in meta])
        ones = np.ones(a.genomes)
        for label, cc, aa in [("free intercept: a=1, c=c_hat", chat, ones), ("free intercept: a=a_hat, c=c_hat", chat, np.maximum(ahat, 1.0)),
                              ("fixed intercept a=1, c=c_hat", chat1, ones), ("true c", np.full(a.genomes, c), ones)]:
            D = aa[None, :] * L[:, None] + cc[None, :] * Q[:, None]
            Db = aa[None, None, :] * Lb[:, :, None] + cc[None, None, :] * Qb[:, :, None]
            df = satterthwaite_df(np.moveaxis(Db, 1, -1), D)
            stat = np.abs(nums) / np.sqrt(D)
            r = dict(gwas_pop=g, pi=pi, c_true=c, c_used=label, n_tests=int(stat.size))
            for al in ALPHAS:
                thr = tdist.isf(al / 2, df)
                r[f"count_{al:g}"] = int(np.count_nonzero(stat > thr))
                r[f"ratio_{al:g}"] = r[f"count_{al:g}"] / stat.size / al
                r[f"pred_ratio_{al:g}"] = float(np.mean(2 * norm.sf(thr * np.sqrt(D / Vt[:, None])))) / al
            rows.append(r)
        print(g, pi, c, f"chat={chat.mean():.4f}±{chat.std():.4f}", f"{time.time() - t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
pd.DataFrame(crow).to_csv(a.out.replace(".csv", "_c_estimates.csv"), index=False)
