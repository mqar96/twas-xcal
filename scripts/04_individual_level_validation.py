import argparse, time, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from sklearn.linear_model import ElasticNetCV
from xcal.genotypes import simulate_region, POPS
from xcal.calibration import Panel
from xcal.core import satterthwaite_df
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser(description="Individual-level GWAS check of the summary-statistic null model, XCAL and XCAL-t")
ap.add_argument("out")
ap.add_argument("--regions", type=int, default=40)
ap.add_argument("--n-gwas", type=int, default=10000)
ap.add_argument("--reps", type=int, default=1000)
a = ap.parse_args()
N = a.n_gwas
CS = [0.0, 0.1, 0.4]
ALPHAS = [0.05, 0.01]


def train(X, rng, h2=0.15, k=2, n=500):
    Xt = X[:n].astype(float); frq = Xt.mean(0) / 2
    cand = np.where(np.minimum(frq, 1 - frq) >= 0.05)[0]
    b = np.zeros(X.shape[1]); b[rng.choice(cand, k, replace=False)] = rng.normal(0, 1, k)
    gv = (Xt - Xt.mean(0)) @ b; gv = gv / gv.std() * np.sqrt(h2)
    y = gv + rng.normal(0, np.sqrt(1 - h2), n)
    sd = Xt.std(0); ok = sd > 0
    Xs = np.zeros_like(Xt); Xs[:, ok] = (Xt[:, ok] - Xt[:, ok].mean(0)) / sd[ok]
    m = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=5, max_iter=2000, tol=1e-3, precompute=True, random_state=0).fit(Xs, y - y.mean())
    w = np.zeros(X.shape[1]); w[ok] = m.coef_[ok] / sd[ok]
    return w


rows = []
t0 = time.time()
import os as _os
_PART = a.out + ".parts"
_os.makedirs(_PART, exist_ok=True)
for r in range(a.regions):
    _pf = f"{_PART}/r{r:03d}.csv"
    if _os.path.exists(_pf):
        if _os.path.getsize(_pf) > 1:
            rows.extend(pd.read_csv(_pf).to_dict("records"))
        continue
    _start = len(rows)
    d = simulate_region(90_000 + r, splits={"train": 1000, "ldref": 500, "gwas": N})
    if int(d["n_snps"]) < 20:
        open(_pf, "w").write("")
        continue
    rng = np.random.default_rng(r)
    P = {p: Panel(d[f"{p}_ldref"]) for p in POPS}
    genes = []
    for tp in ("CEU", "YRI", "CEU", "YRI"):
        w = train(d[f"{tp}_train"], rng)
        if np.any(w):
            genes.append((tp, w))
    for g in POPS:
        Pg = Panel(d[f"{g}_gwas"])
        Xs = Pg.Z.astype(np.float32)
        pre = []
        for gi, (tp, w) in enumerate(genes):
            vT, vE, vG = P[g].v(w), P["CEU"].v(w), w * Pg.s
            mT, mE, mG = P[g].moments(vT), P["CEU"].moments(vE), P[g].moments(vG)
            Lb, Qb = P[g].jackknife(vG, K=20, seed=r)
            pre.append((gi, tp, w, vT, vE, vG, mT, mE, mG, Lb, Qb))
        for c in CS:
            B = rng.normal(0, np.sqrt(c / N), (Xs.shape[1], a.reps)).astype(np.float32)
            Gv = Xs @ B
            Y = Gv + rng.normal(0, 1, (N, a.reps)).astype(np.float32) * np.sqrt(np.clip(1 - Gv.var(0), 0.05, 1))
            Y = (Y - Y.mean(0)) / Y.std(0)
            Zg = (Xs.T @ Y) / np.sqrt(N)
            for gi, tp, w, vT, vE, vG, mT, mE, mG, Lb, Qb in pre:
                LT, LE, LG, QG = mT["L"], mE["L"], mG["L"], mG["Q"]
                qE = mE["Q"] / LE
                df = float(satterthwaite_df(Lb + c * Qb, LG + c * QG))
                numT = vT @ Zg; numE = vE @ Zg; numG = vG @ Zg
                RvT = Pg.R @ vT; RvE = Pg.R @ vE; RvG = Pg.R @ vG
                predT = float(vT @ RvT + c * RvT @ RvT); predE = float(vE @ RvE + c * RvE @ RvE); predG = float(vG @ RvG + c * RvG @ RvG)
                stats = {"STD_M": (numT / np.sqrt(LT), predT / LT), "STD_E": (numE / np.sqrt(LE), predE / LE),
                         "VC_E_E": (numE / np.sqrt(LE * (1 + c * qE)), predE / (LE * (1 + c * qE))),
                         "XCAL": (numG / np.sqrt(LG + c * QG), predG / (LG + c * QG))}
                row = dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, c=c, reps=a.reps, df=df)
                for k, (s, V) in stats.items():
                    row[f"emp_var_{k}"] = float(np.var(s)); row[f"pred_var_{k}"] = V
                    for al in ALPHAS:
                        row[f"rej_{k}_{al:g}"] = int(np.sum(np.abs(s) > norm.isf(al / 2)))
                        row[f"pred_rej_{k}_{al:g}"] = float(2 * norm.sf(norm.isf(al / 2) / np.sqrt(V)))
                s, V = stats["XCAL"]
                for al in ALPHAS:
                    thr = tdist.isf(al / 2, df)
                    row[f"rej_XCALT_{al:g}"] = int(np.sum(np.abs(s) > thr))
                    row[f"pred_rej_XCALT_{al:g}"] = float(2 * norm.sf(thr / np.sqrt(V)))
                rows.append(row)
    pd.DataFrame(rows[_start:]).to_csv(_pf + ".tmp", index=False)
    _os.replace(_pf + ".tmp", _pf)
    print(f"region {r} done {time.time() - t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
