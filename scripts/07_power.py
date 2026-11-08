import argparse, time, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from sklearn.linear_model import ElasticNetCV
from xcal.genotypes import simulate_region, POPS
from xcal.calibration import Panel, true_var, SIM_M_EFF
from xcal.core import satterthwaite_df
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser(description="Power at genome-wide alpha for oracle, XCAL, XCAL-t, EUR-Phi variance control and standard TWAS")
ap.add_argument("out"); ap.add_argument("--regions", type=int, default=60)
a = ap.parse_args()
SCEN = [(50_000, 0.3), (200_000, 0.5)]
Q2 = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3]
ALPHA = 2.5e-6


def pw(mean, var, thr):
    sd = np.sqrt(var)
    return norm.sf((thr - mean) / sd) + norm.cdf((-thr - mean) / sd)


rows = []; t0 = time.time()
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
    d = simulate_region(60_000 + r, splits={"train": 1000, "ldref": 1000, "gwas": 20000})
    if int(d["n_snps"]) < 20:
        open(_pf, "w").write("")
        continue
    rng = np.random.default_rng(r)
    Pg = {g: Panel(d[f"{g}_gwas"]) for g in POPS}
    for gi in range(4):
        tp = "CEU" if gi % 2 == 0 else "YRI"
        X = d[f"{tp}_train"][:500].astype(float); frq = X.mean(0) / 2
        cand = np.where(np.minimum(frq, 1 - frq) >= 0.05)[0]
        b = np.zeros(X.shape[1]); k = int(rng.choice([1, 2, 3])); b[rng.choice(cand, k, replace=False)] = rng.normal(0, 1, k)
        h2 = float(rng.choice([0.1, 0.2]))
        gv = (X - X.mean(0)) @ b; gv = gv / gv.std() * np.sqrt(h2); y = gv + rng.normal(0, np.sqrt(1 - h2), len(gv))
        sd = X.std(0); ok = sd > 0; Xs = np.zeros_like(X); Xs[:, ok] = (X[:, ok] - X[:, ok].mean(0)) / sd[ok]
        m = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=5, max_iter=2000, tol=1e-3, precompute=True, random_state=0).fit(Xs, y - y.mean())
        w = np.zeros(X.shape[1]); w[ok] = m.coef_[ok] / sd[ok]
        if not np.any(w):
            continue
        for g in POPS:
            P = Pg[g]; Pt = Panel(d[f"{g}_ldref"]); Pe = Panel(d["CEU_ldref"])
            bs = b * P.s; gvar = float(bs @ P.R @ bs)
            if gvar <= 0:
                continue
            mu_unit = P.R @ bs / np.sqrt(gvar)
            vG = w * P.s; vt = Pt.v(w); ve = Pe.v(w)
            mG = Pt.moments(vG); mt = Pt.moments(vt); me = Pe.moments(ve)
            if min(mG["L"], mt["L"], me["L"]) <= 1e-12:
                continue
            Lb, Qb = Pt.jackknife(vG, K=20, seed=r)
            qE = me["Q"] / me["L"]
            for N, h2t in SCEN:
                c = N * h2t / SIM_M_EFF
                varG = true_var(vG, P.R, c); varT = true_var(vt, P.R, c)
                Dx = mG["L"] + c * mG["Q"]
                df = float(satterthwaite_df(Lb + c * Qb, Dx))
                for q2 in Q2:
                    mG_ = np.sqrt(N * q2) * float(vG @ mu_unit); mT_ = np.sqrt(N * q2) * float(vt @ mu_unit)
                    zthr = norm.isf(ALPHA / 2); tthr = tdist.isf(ALPHA / 2, df)
                    spec = {"Oracle": (mG_, varG, varG), "XCAL": (mG_, varG, Dx),
                            "VCeur": (mT_, varT, mt["L"] * (1 + c * qE)), "Standard": (mT_, varT, mt["L"])}
                    res = {f"power_{k}": pw(mu / np.sqrt(den), var / den, zthr) for k, (mu, var, den) in spec.items()}
                    res["power_XCALt"] = pw(mG_ / np.sqrt(Dx), varG / Dx, tthr)
                    res.update({f"type1_{k}": pw(0.0, var / den, zthr) for k, (mu, var, den) in spec.items()})
                    res["type1_XCALt"] = pw(0.0, varG / Dx, tthr)
                    rows.append(dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, N=N, h2_trait=h2t, c=c, q2=q2, df=df, **res))
    pd.DataFrame(rows[_start:]).to_csv(_pf + ".tmp", index=False)
    _os.replace(_pf + ".tmp", _pf)
    print(f"region {r} {time.time() - t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
