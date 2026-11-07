import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.core import ld_from_genotypes, unbiased_quad_R2, xcal_gene, run_xcal, estimate_c


def _ar1_genotypes(n, p, rho, rng):
    L = np.zeros((n, p)); L[:, 0] = rng.normal(size=n)
    for j in range(1, p):
        L[:, j] = rho * L[:, j - 1] + np.sqrt(1 - rho ** 2) * rng.normal(size=n)
    return (L > 0).astype(float) + (rng.normal(size=(n, p)) * 0.3 + L > 0.8)


def test_unbiased_quad_closer_to_truth():
    rng = np.random.default_rng(0)
    p = 200
    Sigma = 0.8 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
    v = np.zeros(p); v[[50, 51, 120]] = [1.0, -0.5, 0.7]
    truth = v @ Sigma @ Sigma @ v
    raw, unb = [], []
    for _ in range(50):
        X = rng.multivariate_normal(np.zeros(p), Sigma, size=400)
        R = np.corrcoef(X, rowvar=False)
        raw.append((R @ v) @ (R @ v)); unb.append(unbiased_quad_R2(R, 400, v))
    assert abs(np.mean(unb) - truth) < abs(np.mean(raw) - truth)
    assert abs(np.mean(unb) / truth - 1) < 0.05


def test_xcal_calibrated_under_polygenic_null():
    rng = np.random.default_rng(1)
    p, c = 150, 0.4
    Sigma = 0.9 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
    Xp = rng.multivariate_normal(np.zeros(p), Sigma, size=500)
    R, s, n = ld_from_genotypes(Xp)
    w = np.zeros(p); w[[40, 42, 90]] = [0.6, 0.3, -0.4]
    Lc = np.linalg.cholesky(Sigma + c * Sigma @ Sigma + 1e-9 * np.eye(p))
    zs = []
    for _ in range(4000):
        z = Lc @ rng.normal(size=p)
        zs.append(xcal_gene(w, z, R, s, n, c)[0])
    assert 0.85 < np.var(zs) < 1.15


def test_estimate_c_recovers_slope():
    rng = np.random.default_rng(2)
    l = rng.uniform(1, 60, 20000)
    chi2 = (1 + 0.05 * l) * rng.chisquare(1, size=l.size)
    est = estimate_c(np.sqrt(chi2), l)
    assert abs(est["c"] - 0.05) < 0.01 and abs(est["intercept"] - 1) < 0.1


def test_run_xcal_end_to_end():
    rng = np.random.default_rng(3)
    p = 60
    G = rng.binomial(2, 0.3, size=(300, p)).astype(float)
    snps = [f"rs{i}" for i in range(p)]
    pos = np.arange(p) * 10_000
    weights = pd.DataFrame(dict(gene="G1", snp=["rs5", "rs6"], weight=[0.2, -0.1]))
    gwas = pd.DataFrame(dict(snp=snps, z=rng.normal(size=p)))
    out = run_xcal(weights, gwas, G, snps, pos, c=0.1)
    assert len(out) == 1 and out.p_xcal.between(0, 1).all() and out.inflation_factor.iloc[0] >= 1
