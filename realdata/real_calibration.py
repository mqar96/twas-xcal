import sys, time, numpy as np, pandas as pd
from bed_reader import open_bed
G = "/mnt/user-data/uploads/Downloads/1000G_{}"
POPS = ["AFR", "EUR", "EAS", "SAS", "AMR"]
SETS = ["AFA", "CAU", "HIS"]
WIN = 1_000_000
K = 20
out = sys.argv[1]
chroms = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else list(range(1, 23))

bim = pd.read_csv(G.format("EUR") + ".bim", sep="\t", header=None, names=["chr", "rsid", "cm", "pos", "a1", "a2"])
bim["idx"] = np.arange(len(bim))
famn = {p: sum(1 for _ in open(G.format(p) + ".fam")) for p in POPS}
rng = np.random.default_rng(2026)
halves = {}
for p in POPS:
    perm = rng.permutation(famn[p]); h = famn[p] // 2
    halves[p] = (np.sort(perm[:h]), np.sort(perm[h:]))
weights = {s: pd.read_csv(f"mesa_{s}_weights.csv.gz") for s in SETS}
extra = {s: pd.read_csv(f"mesa_{s}_extra.csv")[["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval"]] for s in SETS}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def std_panel(X):
    X = X.astype(np.float32)
    miss = X < 0
    if miss.any():
        X[miss] = np.nan
        mu = np.nanmean(X, 0); X = np.where(np.isnan(X), mu, X)
    mu = X.mean(0); sd = X.std(0)
    ok = sd > 0
    Z = np.zeros_like(X); Z[:, ok] = (X[:, ok] - mu[ok]) / sd[ok]
    return Z, sd


def stats(Z, S_local, v):
    n, p = Z.shape
    u = Z[:, S_local] @ v
    lin = float(u @ u) / n
    Sv = Z.T @ u / n
    quad = float(Sv @ Sv)
    unb = max((quad - (p / n) * lin) / (1 + 1 / n), lin)
    return lin, unb


import os
t0 = time.time()
CHUNK = 25_000_000; PAD = 5_000_000
for ch in chroms:
    fout = f"{out}.chr{ch}.csv"
    if os.path.exists(fout):
        continue
    rows = []
    bc = bim[bim.chr == ch].reset_index(drop=True)
    gpos = {}
    for s_ in SETS:
        w_ = weights[s_]; w_ = w_.merge(bc[["rsid", "pos"]], on="rsid")
        gpos[s_] = w_.groupby("gene").pos.min()
    for cs in range(0, int(bc.pos.max()) + 1, CHUNK):
        ce = cs + CHUNK
        sel = (bc.pos >= cs - PAD) & (bc.pos < ce + PAD)
        if sel.sum() == 0: continue
        b = bc[sel].reset_index(drop=True)
        idx = b.idx.values
        P = {}
        for p in POPS:
            with open_bed(G.format(p) + ".bed", count_A1=True, num_threads=2) as bed:
                X = bed.read(index=np.s_[:, idx], dtype="int8")
            A, B = halves[p]
            P[p] = {"full": std_panel(X), "A": std_panel(X[A]), "B": std_panel(X[B])}
            del X
        rs2loc = dict(zip(b.rsid, range(len(b))))
        blocks = {p: np.array_split(np.arange(len(halves[p][0])), K) for p in POPS}
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
                if keep.sum() == 0: continue
                loc = loc[keep]; wa = gw.weight.values[keep] * sign[keep]
                pos = b.pos.values[loc]
                lo = np.searchsorted(b.pos.values, pos.min() - WIN); hi = np.searchsorted(b.pos.values, pos.max() + WIN, side="right")
                win = np.arange(lo, hi)
                Sl = loc - lo
                panels = {p: {k: (P[p][k][0][:, lo:hi], P[p][k][1][lo:hi]) for k in ("full", "A", "B")} for p in POPS}
                ZEf, sEf = panels["EUR"]["full"]; ZEa, sEa = panels["EUR"]["A"]
                vEf = wa * sEf[Sl]; vEa = wa * sEa[Sl]
                linEf, unbEf = stats(ZEf, Sl, vEf); linEa, unbEa = stats(ZEa, Sl, vEa)
                if linEf <= 1e-12 or linEa <= 1e-12: continue
                for p in POPS:
                    ZA, sA = panels[p]["A"]; ZB, sB = panels[p]["B"]
                    vA = wa * sA[Sl]
                    linA, unbA = stats(ZA, Sl, vA)
                    linB_vA, unbB_vA = stats(ZB, Sl, vA)
                    vE, linE, unbE = (vEa, linEa, unbEa) if p == "EUR" else (vEf, linEf, unbEf)
                    linB_vE, unbB_vE = stats(ZB, Sl, vE)
                    jk = [stats(np.delete(ZA, blk, axis=0), Sl, vA) for blk in blocks[p]]
                    Zf, sf = panels[p]["full"]; vf = wa * sf[Sl]
                    linF, unbF = stats(Zf, Sl, vf)
                    if min(linA, linF, linB_vA, linB_vE) <= 1e-12: continue
                    rows.append(dict(set=s, gene=gene, chr=ch, pop=p, n_w=len(loc), n_win=len(win),
                                     linA=linA, unbA=unbA, linB_vA=linB_vA, unbB_vA=unbB_vA,
                                     linE=linE, unbE=unbE, linB_vE=linB_vE, unbB_vE=unbB_vE,
                                     q_full=unbF / linF, qE=unbE / linE,
                                     jk_lin=";".join(f"{x[0]:.6g}" for x in jk), jk_unb=";".join(f"{x[1]:.6g}" for x in jk)))
        del P
    pd.DataFrame(rows).to_csv(fout, index=False)
    print(f"chr{ch} done, rows {len(rows)}, {time.time()-t0:.0f}s", flush=True)
import glob
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{out}.chr*.csv"))], ignore_index=True)
for s in SETS:
    df.loc[df.set == s, "pred_perf_pval"] = df.loc[df.set == s, "gene"].map(extra[s].set_index("gene")["pred.perf.pval"])
    df.loc[df.set == s, "pred_perf_R2"] = df.loc[df.set == s, "gene"].map(extra[s].set_index("gene")["pred.perf.R2"])
    df.loc[df.set == s, "genename"] = df.loc[df.set == s, "gene"].map(extra[s].set_index("gene")["genename"])
df.to_csv(out, index=False)
print("done", len(df), "rows", df.groupby("set").gene.nunique().to_dict())
