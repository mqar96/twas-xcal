import argparse, glob, hashlib, json, os, time
import numpy as np
import pandas as pd
from bed_reader import open_bed
from xcal.core import panel_moments, jackknife_moments

ap = argparse.ArgumentParser(description="Real-genotype benchmark: 1000 Genomes phase 3 x MESA DivPop expression models (split-half panels)")
ap.add_argument("--geno", required=True, help="PLINK prefix template with {} for the super-population, e.g. /data/1000G_{}")
ap.add_argument("--models", required=True, help="directory with mesa_{AFA,CAU,HIS}_weights.csv.gz and mesa_*_extra.csv (from read_divpop.py)")
ap.add_argument("--out", required=True, help="output prefix, e.g. out/real")
ap.add_argument("--chroms", default=",".join(str(i) for i in range(1, 23)))
ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]
SETS = ["AFA", "CAU", "HIS"]
CFG = dict(version="0.2.0", window_bp=1_000_000, K=20, chunk=25_000_000, pad=5_000_000, split_seed=2026, pops=POPS, sets=SETS,
           models_md5={s: hashlib.md5(open(f"{a.models}/mesa_{s}_weights.csv.gz", "rb").read()).hexdigest() for s in SETS})
CFG_HASH = hashlib.sha256(json.dumps(CFG, sort_keys=True).encode()).hexdigest()[:16]
json.dump(dict(CFG, config_hash=CFG_HASH), open(f"{a.out}.config.json", "w"), indent=1)
chroms = [int(x) for x in a.chroms.split(",")]

bim = pd.read_csv(a.geno.format("EUR") + ".bim", sep="\t", header=None, names=["chr", "rsid", "cm", "pos", "a1", "a2"])
bim["idx"] = np.arange(len(bim))
famn = {p: sum(1 for _ in open(a.geno.format(p) + ".fam")) for p in POPS}
rng = np.random.default_rng(CFG["split_seed"])
halves = {}
for p in POPS:
    perm = rng.permutation(famn[p]); h = famn[p] // 2
    halves[p] = (np.sort(perm[:h]), np.sort(perm[h:]))
weights = {s: pd.read_csv(f"{a.models}/mesa_{s}_weights.csv.gz") for s in SETS}
extra = {s: pd.read_csv(f"{a.models}/mesa_{s}_extra.csv")[["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval"]] for s in SETS}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def std_panel(X):
    X = X.astype(np.float64)
    miss = X < 0
    if miss.any():
        X[miss] = np.nan
        X = np.where(np.isnan(X), np.nanmean(X, 0), X)
    mu = X.mean(0); sd = X.std(0)
    ok = sd > 0
    Z = np.zeros_like(X); Z[:, ok] = (X[:, ok] - mu[ok]) / sd[ok]
    return Z, sd, ok


def moments(panel, lo, hi, loc, wa):
    Z, sd, ok = panel
    Zw, sw, okw = Z[:, lo:hi], sd[lo:hi], ok[lo:hi]
    v = np.zeros(hi - lo); v[loc] = wa * sw[loc]
    return panel_moments(Zw, v, okw), Zw, sw, okw, v


t0 = time.time()
for ch in chroms:
    fout = f"{a.out}.chr{ch}.csv"
    if os.path.exists(fout):
        prev = pd.read_csv(fout, nrows=1)
        if len(prev) and prev.get("config_hash", pd.Series([None])).iloc[0] == CFG_HASH:
            continue
    rows = []
    bc = bim[bim.chr == ch].reset_index(drop=True)
    gpos = {s_: weights[s_].merge(bc[["rsid", "pos"]], on="rsid").groupby("gene").pos.min() for s_ in SETS}
    for cs in range(0, int(bc.pos.max()) + 1, CFG["chunk"]):
        ce = cs + CFG["chunk"]
        sel = (bc.pos >= cs - CFG["pad"]) & (bc.pos < ce + CFG["pad"])
        if sel.sum() == 0:
            continue
        b = bc[sel].reset_index(drop=True)
        idx = b.idx.values
        P = {}
        for p in POPS:
            with open_bed(a.geno.format(p) + ".bed", count_A1=True, num_threads=a.threads) as bed:
                X = bed.read(index=np.s_[:, idx], dtype="int8")
            A, B = halves[p]
            P[p] = {"full": std_panel(X), "A": std_panel(X[A]), "B": std_panel(X[B])}
            del X
        rs2loc = dict(zip(b.rsid, range(len(b))))
        blocks = {p: np.array_split(np.random.default_rng(ch * 1000 + cs // CFG["chunk"]).permutation(len(halves[p][0])), CFG["K"]) for p in POPS}
        for s in SETS:
            w = weights[s]
            w = w[w.rsid.isin(rs2loc) & w.gene.isin(gpos[s][(gpos[s] >= cs) & (gpos[s] < ce)].index)]
            for gene, gw in w.groupby("gene"):
                loc = np.array([rs2loc[r] for r in gw.rsid])
                a1 = b.a1.values[loc]; a2 = b.a2.values[loc]
                eff = gw.eff_allele.values; ref = gw.ref_allele.values
                sign = np.zeros(len(loc))
                for i in range(len(loc)):
                    if eff[i] == a1[i] and ref[i] == a2[i]: sign[i] = 1
                    elif eff[i] == a2[i] and ref[i] == a1[i]: sign[i] = -1
                    elif COMP.get(eff[i]) == a1[i] and COMP.get(ref[i]) == a2[i]: sign[i] = 1
                    elif COMP.get(eff[i]) == a2[i] and COMP.get(ref[i]) == a1[i]: sign[i] = -1
                keep = sign != 0
                if keep.sum() == 0:
                    continue
                loc = loc[keep]; wa = gw.weight.values[keep] * sign[keep]
                pos = b.pos.values[loc]
                lo = np.searchsorted(b.pos.values, pos.min() - CFG["window_bp"])
                hi = np.searchsorted(b.pos.values, pos.max() + CFG["window_bp"], side="right")
                Sl = loc - lo
                mEf = moments(P["EUR"]["full"], lo, hi, Sl, wa)[0]
                mEa = moments(P["EUR"]["A"], lo, hi, Sl, wa)[0]
                if mEf["L"] <= 1e-12 or mEa["L"] <= 1e-12:
                    continue
                for p in POPS:
                    mA, ZA, sA, okA, vA = moments(P[p]["A"], lo, hi, Sl, wa)
                    ZB, sB, okB = P[p]["B"][0][:, lo:hi], P[p]["B"][1][lo:hi], P[p]["B"][2][lo:hi]
                    mB_vA = panel_moments(ZB, vA, okB)
                    if p == "EUR":
                        mE, vE = mEa, np.zeros(hi - lo)
                        vE[Sl] = wa * P["EUR"]["A"][1][lo:hi][Sl]
                    else:
                        mE, vE = mEf, np.zeros(hi - lo)
                        vE[Sl] = wa * P["EUR"]["full"][1][lo:hi][Sl]
                    mB_vE = panel_moments(ZB, vE, okB)
                    mF = moments(P[p]["full"], lo, hi, Sl, wa)[0]
                    wwin = np.zeros(hi - lo); wwin[Sl] = wa
                    Lb, Qb = jackknife_moments(ZA, sA, np.where(okA, wwin, 0.0), blocks=blocks[p])
                    if min(mA["L"], mF["L"], mB_vA["L"], mB_vE["L"]) <= 1e-12:
                        continue
                    rows.append(dict(set=s, gene=gene, chr=ch, pop=p, n_w=len(loc), n_win=hi - lo,
                                     linA=mA["L"], unbA=mA["Q"], floorA=mA["floor_active"],
                                     linB_vA=mB_vA["L"], unbB_vA=mB_vA["Q"], linE=mE["L"], unbE=mE["Q"],
                                     linB_vE=mB_vE["L"], unbB_vE=mB_vE["Q"], q_full=mF["Q"] / mF["L"], qE=mE["Q"] / mE["L"],
                                     floorF=mF["floor_active"], jk_lin=";".join(f"{x:.8g}" for x in Lb),
                                     jk_unb=";".join(f"{x:.8g}" for x in Qb), config_hash=CFG_HASH))
        del P
    pd.DataFrame(rows).to_csv(fout + ".tmp", index=False)
    os.replace(fout + ".tmp", fout)
    print(f"chr{ch} done, rows {len(rows)}, {time.time() - t0:.0f}s", flush=True)

df = pd.concat([pd.read_csv(f"{a.out}.chr{c}.csv") for c in chroms], ignore_index=True)
for s in SETS:
    e = extra[s].set_index("gene")
    m = df.set == s
    df.loc[m, "pred_perf_pval"] = df.loc[m, "gene"].map(e["pred.perf.pval"])
    df.loc[m, "pred_perf_R2"] = df.loc[m, "gene"].map(e["pred.perf.R2"])
    df.loc[m, "genename"] = df.loc[m, "gene"].map(e["genename"])
df.to_csv(a.out + ".csv.gz", index=False)
print("done", len(df), "rows", df.groupby("set").gene.nunique().to_dict())
