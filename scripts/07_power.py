import sys, os, time, warnings, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.calibration import panel_stats, quad_R2_unbiased, jackknife_panels, LQ_from_panels, satterthwaite_df
from scipy.stats import norm, t as tdist
from sklearn.linear_model import ElasticNetCV
import stdpopsim
warnings.filterwarnings("ignore")
out = sys.argv[1]; n_regions = int(sys.argv[2])
POPS = ("YRI", "CEU", "CHB")
SPL = {"train": 1000, "ldref": 1000, "gwas": 3000}
M_EFF = 305726
SCEN = [(50_000, 0.3), (200_000, 0.5)]
Q2 = [2e-4, 5e-4, 1e-3, 2e-3]
ALPHA = 2.5e-6


def sim(seed):
    sp = stdpopsim.get_species("HomSap"); m = sp.get_demographic_model("OutOfAfrica_3G09")
    ts = stdpopsim.get_engine("msprime").simulate(m, sp.get_contig(length=1e6), {p: sum(SPL.values()) for p in POPS}, seed=seed)
    pop = np.array([ts.node(n).population for n in ts.samples()])
    names = [ts.population(i).metadata.get("name", str(i)) for i in range(ts.num_populations)]
    idx = {p: np.where(pop == names.index(p))[0] for p in POPS}
    keep = []
    for v in ts.variants():
        g = v.genotypes
        if g.max() > 1: continue
        f = np.array([g[idx[p]].mean() for p in POPS]); maf = np.minimum(f, 1 - f)
        if maf.max() >= 0.05 and maf.min() >= 0.01:
            keep.append(g.astype(np.int8))
    rng = np.random.default_rng(seed)
    sel = np.sort(rng.choice(len(keep), min(400, len(keep)), replace=False))
    H = np.vstack([keep[i] for i in sel]); res = {}
    for p in POPS:
        h = H[:, idx[p]]; dip = (h[:, 0::2] + h[:, 1::2]).T; o = 0
        for k, n in SPL.items():
            res[f"{p}_{k}"] = dip[o:o + n]; o += n
    return res


def pw(mean, var, thr):
    sd = np.sqrt(var)
    return norm.sf((thr - mean) / sd) + norm.cdf((-thr - mean) / sd)

rows = []; t0 = time.time()
for r in range(n_regions):
    d = sim(60_000 + r); rng = np.random.default_rng(r)
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
        if not np.any(w): continue
        for g in POPS:
            Pg = panel_stats(d[f"{g}_gwas"]); Pt = panel_stats(d[f"{g}_ldref"]); Pe = panel_stats(d["CEU_ldref"])
            JK = jackknife_panels(d[f"{g}_ldref"], K=20, seed=r)
            bs = b * Pg["s"]; gvar = bs @ Pg["R"] @ bs
            if gvar <= 0: continue
            mu_unit = Pg["R"] @ bs / np.sqrt(gvar)
            vt = w * Pt["s"]; Lf = vt @ Pt["R"] @ vt; Qf = quad_R2_unbiased(Pt, vt)
            ue = w * Pe["s"]; qE = quad_R2_unbiased(Pe, ue) / (ue @ Pe["R"] @ ue)
            Lk, Qk = LQ_from_panels(w, JK)
            for N, h2t in SCEN:
                c = N * h2t / M_EFF
                num_var = vt @ Pg["R"] @ vt + c * quad_R2_unbiased(Pg, vt)
                df = satterthwaite_df(Lk, Qk, c, Lf + c * Qf)
                for q2 in Q2:
                    mean_num = np.sqrt(N * q2) * (vt @ mu_unit)
                    zthr = norm.isf(ALPHA / 2); tthr = tdist.isf(ALPHA / 2, df)
                    D = {"Oracle": num_var, "XCAL": Lf + c * Qf, "VCeur": Lf * (1 + c * qE), "Standard": Lf}
                    res = {f"power_{k}": pw(mean_num / np.sqrt(v), num_var / v, zthr) for k, v in D.items()}
                    res["power_XCALt"] = pw(mean_num / np.sqrt(D["XCAL"]), num_var / D["XCAL"], tthr)
                    rows.append(dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, N=N, h2_trait=h2t, c=c, q2=q2, df=df, **res))
    print(f"region {r} {time.time()-t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(out, index=False)
