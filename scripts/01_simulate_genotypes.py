import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.genotypes import simulate_region
from joblib import Parallel, delayed

out = sys.argv[1]; n = int(sys.argv[2]); start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
os.makedirs(out, exist_ok=True)

def one(i):
    f = f"{out}/region_{i:04d}.npz"
    if os.path.exists(f):
        return 0
    r = simulate_region(seed=1000 + i)
    np.savez_compressed(f, **r)
    return int(r["n_snps"])

t = time.time()
res = Parallel(n_jobs=2)(delayed(one)(i) for i in range(start, start + n))
print("done", n, "regions", time.time() - t, "s; snps/region min/med/max", min(res), np.median(res), max(res))
