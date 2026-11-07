import sys, os, json, time, warnings, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.calibration import panel_stats, quad_R2_unbiased, poly_q
from scipy.stats import norm
from sklearn.linear_model import ElasticNetCV
from sklearn.exceptions import ConvergenceWarning
import stdpopsim
warnings.filterwarnings("ignore")

out = sys.argv[1]; n_regions = int(sys.argv[2]); N = int(sys.argv[3]); REPS = int(sys.argv[4])
POPS = ("YRI", "CEU", "CHB")
SPL = {"train": 1000, "ldref": 500, "gwas": N}
CS = [0.0, 0.1, 0.4]


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
    H = np.vstack([keep[i] for i in sel])
    res = {}
    for p in POPS:
        h = H[:, idx[p]]; dip = (h[:, 0::2] + h[:, 1::2]).T
        o = 0
        for k, n in SPL.items():
            res[f"{p}_{k}"] = dip[o:o + n]; o += n
    return res


def train(X, rng, h2=0.1, k=2, n=500):
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
t = time.time()
for r in range(n_regions):
    d = sim(90_000 + r)
    rng = np.random.default_rng(r)
    P = {f"{p}_ldref": panel_stats(d[f"{p}_ldref"]) for p in POPS}
    genes = []
    for tp in ("CEU", "YRI", "CEU", "YRI"):
        w = train(d[f"{tp}_train"], rng, h2=0.15)
        if np.any(w): genes.append((tp, w))
    for g in POPS:
        X = d[f"{g}_gwas"].astype(np.float32)
        s = X.std(0); ok = s > 0
        Xs = np.zeros_like(X); Xs[:, ok] = (X[:, ok] - X[:, ok].mean(0)) / s[ok]
        Pg = panel_stats(d[f"{g}_gwas"])
        for c in CS:
            B = rng.normal(0, np.sqrt(c / N), (X.shape[1], REPS)).astype(np.float32) if c > 0 else np.zeros((X.shape[1], REPS), np.float32)
            G = Xs @ B
            E = rng.normal(0, 1, (N, REPS)).astype(np.float32)
            Y = G + E * np.sqrt(np.clip(1 - G.var(0), 0.05, 1))
            Y = (Y - Y.mean(0)) / Y.std(0)
            Zg = (Xs.T @ Y) / np.sqrt(N)
            for gi, (tp, w) in enumerate(genes):
                for refname in (f"{g}_ldref", "CEU_ldref"):
                    Pr = P[refname]
                    v = w * Pr["s"]
                    Zt = (v @ Zg) / np.sqrt(v @ Pr["R"] @ v)
                    qE = poly_q(w, P["CEU_ldref"]); qT = poly_q(w, P[f"{g}_ldref"])
                    vt = w * P[f"{g}_ldref"]["s"]
                    Zx = (vt @ Zg) / np.sqrt(vt @ P[f"{g}_ldref"]["R"] @ vt + c * quad_R2_unbiased(P[f"{g}_ldref"], vt))
                    Rv = Pg["R"] @ v
                    V_pred = (v @ Rv + c * N / N * (Rv @ Rv)) / (v @ Pr["R"] @ v)
                    rows.append(dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, c=c, ref=refname.split("_")[0],
                                     emp_var_unadj=float(np.var(Zt)), pred_var_unadj=float(V_pred),
                                     t1_unadj=float(np.mean(np.abs(Zt) > 1.959964)),
                                     t1_VCeur=float(np.mean(np.abs(Zt / np.sqrt(1 + c * qE)) > 1.959964)),
                                     t1_VCtarget=float(np.mean(np.abs(Zt / np.sqrt(1 + c * qT)) > 1.959964)),
                                     t1_XCAL=float(np.mean(np.abs(Zx) > 1.959964)),
                                     emp_var_XCAL=float(np.var(Zx)), reps=REPS))
    print(f"region {r} done {time.time()-t:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(out, index=False)
