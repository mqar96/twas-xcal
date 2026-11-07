import sys, os, time, warnings, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.calibration import panel_stats, quad_R2_unbiased, jackknife_panels, LQ_from_panels, satterthwaite_df
from sklearn.linear_model import ElasticNetCV
import stdpopsim
warnings.filterwarnings("ignore")
out = sys.argv[1]; n_regions = int(sys.argv[2])
POPS = ("YRI", "CEU", "CHB")
SPL = {"train": 500, "ldref": 5000, "gwas": 5000}
SIZES = [500, 1000, 2000, 5000]
CS = [0.1, 0.4, 0.8]


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
        h = H[:, idx[p]]; dip = (h[:, 0::2] + h[:, 1::2]).T; o = 0
        for k, n in SPL.items():
            res[f"{p}_{k}"] = dip[o:o + n]; o += n
    return res


def train(X, rng, h2=0.1, k=2):
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
for r in range(n_regions):
    d = sim(70_000 + r); rng = np.random.default_rng(r)
    genes = [(tp, train(d[f"{tp}_train"], rng, h2=float(rng.choice([0.05, 0.1, 0.2])), k=int(rng.choice([1, 2, 3])))) for tp in ("CEU", "YRI", "CEU", "YRI")]
    genes = [(tp, w) for tp, w in genes if np.any(w)]
    for g in POPS:
        Pg = panel_stats(d[f"{g}_gwas"])
        for m in SIZES:
            Pt = panel_stats(d[f"{g}_ldref"][:m]); Pe = panel_stats(d["CEU_ldref"][:m])
            JK = jackknife_panels(d[f"{g}_ldref"][:m], K=20, seed=r)
            for gi, (tp, w) in enumerate(genes):
                vt = w * Pt["s"]; Lf = vt @ Pt["R"] @ vt; Qf = quad_R2_unbiased(Pt, vt)
                ue = w * Pe["s"]; qE = quad_R2_unbiased(Pe, ue) / (ue @ Pe["R"] @ ue)
                Lk, Qk = LQ_from_panels(w, JK)
                for c in CS:
                    num = vt @ Pg["R"] @ vt + c * quad_R2_unbiased(Pg, vt)
                    V_x = num / (Lf + c * Qf)
                    V_vce = (num / Lf) / (1 + c * qE)
                    rows.append(dict(region=r, gene=gi, train_pop=tp, gwas_pop=g, panel_n=m, c=c, V_XCAL=V_x,
                                     df_XCAL=satterthwaite_df(Lk, Qk, c, Lf + c * Qf), V_VCeur=V_vce))
    print(f"region {r} {time.time()-t:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(out, index=False)
