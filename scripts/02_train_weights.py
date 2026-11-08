import os, glob, time, json, numpy as np, warnings
from joblib import Parallel, delayed
from sklearn.linear_model import ElasticNetCV
from sklearn.exceptions import ConvergenceWarning

import argparse
ap = argparse.ArgumentParser(description="Train elastic-net expression models (CEU/YRI, n = 100-1000) for every simulated region")
ap.add_argument("genodir"); ap.add_argument("outdir")
ap.add_argument("--genes-per-region", type=int, default=8)
ap.add_argument("--jobs", type=int, default=2)
a = ap.parse_args()
gdir, out, genes_per_region = a.genodir, a.outdir, a.genes_per_region
os.makedirs(out, exist_ok=True)
POPS = ("YRI", "CEU", "CHB")


def train_region(f):
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    rid = int(os.path.basename(f)[7:11])
    of = f"{out}/weights_{rid:04d}.npz"
    if os.path.exists(of):
        return 0
    d = np.load(f)
    if int(d["n_snps"]) < 20:
        np.savez_compressed(of, recs=json.dumps([]))
        return 0
    rng = np.random.default_rng(50_000 + rid)
    recs = []
    for g in range(genes_per_region):
        train_pop = ["CEU", "YRI"][g % 2]
        n_eqtl = [100, 250, 500, 1000][(g // 2) % 4]
        h2 = float(rng.choice([0.05, 0.1, 0.2]))
        k = int(rng.choice([1, 2, 3, 5]))
        Xt = d[f"{train_pop}_train"][:n_eqtl].astype(float)
        p = Xt.shape[1]
        frq = d[f"{train_pop}_train"].mean(0) / 2
        cand = np.where(np.minimum(frq, 1 - frq) >= 0.05)[0]
        causal = rng.choice(cand, k, replace=False)
        beta = np.zeros(p)
        beta[causal] = rng.normal(0, 1, k)
        def expr(X, rng_):
            gv = (X - X.mean(0)) @ beta
            vg = gv.var() if gv.var() > 0 else 1.0
            gv = gv / np.sqrt(vg) * np.sqrt(h2)
            return gv + rng_.normal(0, np.sqrt(1 - h2), X.shape[0]), gv
        mu_ref = d[f"{train_pop}_train"].mean(0)
        sd_scale = np.sqrt(((d[f"{train_pop}_train"] - mu_ref) @ beta).var())
        y, _ = expr(Xt, rng)
        sd = Xt.std(0); ok = sd > 0
        Xs = np.zeros_like(Xt); Xs[:, ok] = (Xt[:, ok] - Xt[:, ok].mean(0)) / sd[ok]
        m = ElasticNetCV(l1_ratio=0.5, alphas=20, cv=5, max_iter=2000, tol=1e-3, precompute=True, random_state=0).fit(Xs, y - y.mean())
        w_std = m.coef_
        w = np.zeros(p); w[ok] = w_std[ok] / sd[ok]
        nz = np.flatnonzero(w)
        if len(nz) == 0:
            recs.append(dict(gene=f"r{rid}_g{g}", region=rid, train_pop=train_pop, n_eqtl=n_eqtl, h2=h2, k=k, n_w=0))
            continue
        r2_cv = float(max(0.0, 1 - m.mse_path_.mean(axis=-1).min() / y.var()))
        acc = {}
        for P in POPS:
            Xg = d[f"{P}_gwas_eval"].astype(float)
            gtrue = (Xg - Xg.mean(0)) @ beta / sd_scale * np.sqrt(h2)
            yP = gtrue + rng.normal(0, np.sqrt(1 - h2), len(gtrue))
            pred = (Xg - Xg.mean(0)) @ w
            acc[P] = float(np.corrcoef(pred, yP)[0, 1] ** 2) if pred.std() > 0 else 0.0
        recs.append(dict(gene=f"r{rid}_g{g}", region=rid, train_pop=train_pop, n_eqtl=n_eqtl, h2=h2, k=k,
                         n_w=int(len(nz)), r2_cv=r2_cv, idx=nz.tolist(), w=w[nz].tolist(),
                         causal=causal.tolist(), **{f"r2_{P}": acc[P] for P in POPS}))
    np.savez_compressed(of, recs=json.dumps(recs))
    return len(recs)

files = sorted(glob.glob(f"{gdir}/region_*.npz"))
t = time.time()
r = Parallel(n_jobs=a.jobs)(delayed(train_region)(f) for f in files)
print("trained", sum(r), "genes in", time.time() - t, "s")
