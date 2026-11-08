import argparse, os, time
import numpy as np
from joblib import Parallel, delayed
from xcal.genotypes import simulate_region, POPS
from xcal.core import standardize

ap = argparse.ArgumentParser(description="Simulate 1-Mb regions (YRI/CEU/CHB) with training, LD-reference and 20,000-person GWAS-LD samples")
ap.add_argument("outdir")
ap.add_argument("--n-regions", type=int, default=500)
ap.add_argument("--start", type=int, default=0)
ap.add_argument("--jobs", type=int, default=2)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)


def one(i):
    f = f"{a.outdir}/region_{i:04d}.npz"
    if os.path.exists(f):
        return -1
    r = simulate_region(seed=1000 + i)
    out = {"pos": r["pos"], "n_snps": r["n_snps"]}
    for p in POPS:
        out[f"{p}_train"] = r[f"{p}_train"]
        out[f"{p}_ldref"] = r[f"{p}_ldref"]
        out[f"{p}_gwas_eval"] = r[f"{p}_gwas"][:1100]
        Z, s, ok = standardize(r[f"{p}_gwas"])
        R = Z.T @ Z / Z.shape[0]
        out[f"{p}_R_gwas"] = R.astype(np.float32)
        out[f"{p}_sd_gwas"] = s.astype(np.float64)
        out[f"{p}_n_gwas"] = np.array(Z.shape[0])
    np.savez_compressed(f + ".tmp.npz", **out)
    os.replace(f + ".tmp.npz", f)
    return int(r["n_snps"])

t = time.time()
res = Parallel(n_jobs=a.jobs)(delayed(one)(i) for i in range(a.start, a.start + a.n_regions))
done = [x for x in res if x >= 0]
print("simulated", len(done), "regions in", round(time.time() - t), "s; snps/region", (min(done), int(np.median(done)), max(done)) if done else "-")
