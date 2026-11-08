import argparse, time, warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from xcal.genotypes import simulate_region, POPS
from xcal.calibration import Panel, true_var
from xcal.core import satterthwaite_df
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser(description="Effect of LD-panel size (500-5000) on XCAL / XCAL-t and EUR-Phi variance control")
ap.add_argument("out"); ap.add_argument("--regions", type=int, default=40)
a = ap.parse_args()
SIZES = [500, 1000, 2000, 5000]
CS = [0.1, 0.4, 0.8]


def train(X, rng, h2, k):
    Xt = X.astype(float); frq = Xt.mean(0) / 2
    cand = np.where(np.minimum(frq, 1 - frq) >= 0.05)[0]
    b = np.zeros(X.shape[1]); b[rng.choice(cand, k, replace=False)] = rng.normal(0, 1, k)
    gv = (Xt - Xt.mean(0)) @ b; gv = gv / gv.std() * np.sqrt(h2)
    y = gv + rng.normal(0, np.sqrt(1 - h2), len(gv))
    sd = Xt.std(0); ok = sd > 0
    Xs = np.zeros_like(Xt); Xs[:, ok] = (Xt[:, ok] - Xt[:, ok].mean(0)) / sd[ok]
    m = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=5, max_iter=2000, tol=1e-3, precompute=True, random_state=0).fit(Xs, y - y.mean())
    w = np.zeros(X.shape[1]); w[ok] = m.coef_[ok] / sd[ok]
    return w


rows = []; t = time.time()
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
    d = simulate_region(70_000 + r, splits={"train": 500, "ldref": 5000, "gwas": 20000})
    if int(d["n_snps"]) < 20:
        open(_pf, "w").write("")
        continue
    rng = np.random.default_rng(r)
    genes = [(tp, train(d[f"{tp}_train"], rng, float(rng.choice([0.05, 0.1, 0.2])), int(rng.choice([1, 2, 3])))) for tp in ("CEU", "YRI", "CEU", "YRI")]
    genes = [(tp, w) for tp, w in genes if np.any(w)]
    for g in POPS:
        PG = Panel(d[f"{g}_gwas"]); Rg = PG.R
        for m in SIZES:
            Pt = Panel(d[f"{g}_ldref"][:m]); Pe = Panel(d["CEU_ldref"][:m])
            for gi, (tp, w) in enumerate(genes):
                vt = w * PG.s; ve = Pe.v(w)
                mt = Pt.moments(vt); me = Pe.moments(ve)
                if mt["L"] <= 1e-12 or me["L"] <= 1e-12:
                    continue
                Lb, Qb = Pt.jackknife(vt, K=20, seed=r)
                for c in CS:
                    D = mt["L"] + c * mt["Q"]
                    rows.append(dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, panel_n=m, c=c,
                                     V_XCAL=true_var(vt, Rg, c) / D, df_XCAL=float(satterthwaite_df(Lb + c * Qb, D)),
                                     V_VCeur=true_var(ve, Rg, c) / (me["L"] * (1 + c * me["Q"] / me["L"])),
                                     floor_active=mt["floor_active"]))
    pd.DataFrame(rows[_start:]).to_csv(_pf + ".tmp", index=False)
    _os.replace(_pf + ".tmp", _pf)
    print(f"region {r} {time.time() - t:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
