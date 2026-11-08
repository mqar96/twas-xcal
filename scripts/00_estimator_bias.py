import argparse
import numpy as np
import pandas as pd
from xcal.core import standardize, panel_moments

ap = argparse.ArgumentParser(description="Bias of the quadratic estimator v'Sigma^2 v under the sample-correlation construction")
ap.add_argument("out")
ap.add_argument("--reps", type=int, default=1000)
a = ap.parse_args()
rng = np.random.default_rng(2026)
p = 200
idx = np.arange(p)


def weights(kind):
    w = np.zeros(p)
    if kind == "same_sign":
        w[[60, 62, 65]] = [1.0, 0.5, 0.7]
    elif kind == "opposite_sign":
        w[[60, 61]] = [1.0, -1.0]
    elif kind == "dense":
        w = np.random.default_rng(5).normal(size=p) * 0.1
    elif kind == "low_eigen":
        w = None
    return w


rows = []
for struct in ("AR1_0.5", "AR1_0.9", "AR1_0.98", "block_0.8"):
    if struct.startswith("AR1"):
        rho = float(struct.split("_")[1]); Sig = rho ** np.abs(np.subtract.outer(idx, idx))
    else:
        Sig = np.kron(np.eye(p // 20), np.full((20, 20), 0.8)) + 0.2 * np.eye(p)
    lam, U = np.linalg.eigh(Sig)
    Lc = np.linalg.cholesky(Sig + 1e-12 * np.eye(p))
    for kind in ("same_sign", "opposite_sign", "dense", "low_eigen"):
        w = U[:, 0] if kind == "low_eigen" else weights(kind)
        for n in (250, 500, 1000):
            Q, Qu, Qold, Qraw, T2, L, T1, fl, flold = [], [], [], [], [], [], [], [], []
            mu = rng.uniform(0, 2, p); sd = rng.uniform(0.3, 0.7, p)
            for _ in range(a.reps):
                X = (rng.normal(size=(n, p)) @ Lc.T) * sd + mu
                Z, s, ok = standardize(X)
                v = w * s
                m = panel_moments(Z, v, ok)
                Q.append(m["Q"]); Qu.append(m["Q_untruncated"]); Qraw.append(m["Q_raw"]); Qold.append(max(m["Q_untruncated"], m["L"]))
                L.append(m["L"]); fl.append(m["floor_active"]); flold.append(m["Q_untruncated"] < m["L"])
                T2.append(v @ Sig @ Sig @ v); T1.append(v @ Sig @ v)
            T2m = np.mean(T2)
            rows.append(dict(ld=struct, weights=kind, n=n, reps=a.reps, true_q=float(np.mean(T2) / np.mean(T1)),
                             relbias_L=np.mean(L) / np.mean(T1) - 1,
                             relbias_Q_v02=np.mean(Q) / T2m - 1, relbias_Q_untruncated=np.mean(Qu) / T2m - 1,
                             relbias_Q_v01_floor=np.mean(Qold) / T2m - 1, relbias_Q_uncorrected=np.mean(Qraw) / T2m - 1,
                             rel_rmse_Q_v02=np.sqrt(np.mean((np.array(Q) - np.array(T2)) ** 2)) / T2m,
                             floor_rate_v02=np.mean(fl), floor_rate_v01=np.mean(flold)))
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(a.out, index=False)
