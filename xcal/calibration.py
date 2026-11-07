import numpy as np
from scipy.stats import norm


def panel_stats(X):
    X = X.astype(np.float64)
    s = X.std(0)
    ok = s > 0
    Z = np.zeros_like(X)
    Z[:, ok] = (X[:, ok] - X[:, ok].mean(0)) / s[ok]
    n = X.shape[0]
    R = Z.T @ Z / n
    R[~ok, :] = 0
    R[:, ~ok] = 0
    R[~ok, ~ok] = 1.0
    return dict(R=R, s=s, n=n, ok=ok)


def quad_R2_unbiased(P, v):
    R, n = P["R"], P["n"]
    p = R.shape[0]
    Rv = R @ v
    raw = Rv @ Rv
    lin = v @ Rv
    return max((raw - (p / n) * lin) / (1.0 + 1.0 / n), lin)


def null_var_numerator(v, P_gwas, c):
    return v @ P_gwas["R"] @ v + c * quad_R2_unbiased(P_gwas, v)


def v_std(w, P_ref, P_gwas, c):
    v = w * P_ref["s"]
    return null_var_numerator(v, P_gwas, c) / (v @ P_ref["R"] @ v)


def poly_q(w, P):
    u = w * P["s"]
    return quad_R2_unbiased(P, u) / (u @ P["R"] @ u)


def v_xcal(w, P_target, P_gwas, c):
    vt = w * P_target["s"]
    den = vt @ P_target["R"] @ vt + c * quad_R2_unbiased(P_target, vt)
    return null_var_numerator(vt, P_gwas, c) / den


def type1(V, alpha):
    za = norm.isf(alpha / 2)
    return 2 * norm.sf(za / np.sqrt(np.asarray(V, float)))


def jackknife_LQ(w, X, K=20, seed=0):
    X = X.astype(np.float64)
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    blocks = np.array_split(rng.permutation(n), K)
    L, Q = [], []
    for b in blocks:
        keep = np.setdiff1d(np.arange(n), b)
        P = panel_stats(X[keep])
        v = w * P["s"]
        L.append(v @ P["R"] @ v)
        Q.append(quad_R2_unbiased(P, v))
    return np.array(L), np.array(Q)


def satterthwaite_df(L, Q, c, D_full):
    D = L + c * Q
    K = len(D)
    var = (K - 1) / K * np.sum((D - D.mean()) ** 2)
    tau2 = var / D_full ** 2
    return float(np.clip(2.0 / max(tau2, 1e-8), 4.0, 1e5))


def type1_t(V, df, alpha):
    from scipy.stats import t as tdist
    thr = tdist.isf(alpha / 2, np.asarray(df, float))
    return 2 * norm.sf(thr / np.sqrt(np.asarray(V, float)))


def jackknife_panels(X, K=20, seed=0):
    X = X.astype(np.float64)
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    blocks = np.array_split(rng.permutation(n), K)
    return [panel_stats(X[np.setdiff1d(np.arange(n), b)]) for b in blocks]


def LQ_from_panels(w, panels):
    L, Q = [], []
    for P in panels:
        v = w * P["s"]
        L.append(v @ P["R"] @ v)
        Q.append(quad_R2_unbiased(P, v))
    return np.array(L), np.array(Q)
