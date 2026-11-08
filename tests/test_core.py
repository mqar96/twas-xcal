import numpy as np
import pandas as pd
import pytest

from xcal import (standardize, panel_moments, jackknife_moments, xcal_gene, run_xcal, estimate_c,
                  ld_scores, genomic_control, harmonize, sd_from_eaf, xcal_from_moments)
from scipy.stats import norm


def gaussian_panel(Sig, n, rng, mu=None, sd=None):
    p = Sig.shape[0]
    L = np.linalg.cholesky(Sig + 1e-12 * np.eye(p))
    X = rng.normal(size=(n, p)) @ L.T
    if sd is not None:
        X = X * sd
    if mu is not None:
        X = X + mu
    return X


def ar1(p, rho):
    i = np.arange(p)
    return rho ** np.abs(np.subtract.outer(i, i))


def test_audit_counterexample_q_below_one():
    rng = np.random.default_rng(0)
    Sig = np.array([[1.0, 0.8], [0.8, 1.0]])
    X = gaussian_panel(Sig, 200_000, rng)
    Z, s, ok = standardize(X)
    v = np.array([1.0, -1.0]) * s
    m = panel_moments(Z, v, ok)
    assert abs(m["L"] / (v @ v) * 2 - 0.4) < 0.01
    assert abs(m["Q"] / m["L"] - 0.2) < 0.01
    assert not m["floor_active"]


@pytest.mark.parametrize("rho,w_kind", [(0.5, "opp"), (0.9, "pos"), (0.5, "dense")])
def test_quadratic_estimator_nearly_unbiased(rho, w_kind):
    rng = np.random.default_rng(1)
    p, n = 120, 300
    Sig = ar1(p, rho)
    w = {"opp": np.r_[np.zeros(40), 1, -1, np.zeros(p - 42)],
         "pos": np.r_[np.zeros(40), 1, .5, .7, np.zeros(p - 43)],
         "dense": rng.normal(size=p) * 0.1}[w_kind]
    mu = rng.uniform(0, 2, p)
    sd = rng.uniform(0.3, 0.7, p)
    est, tru, raw = [], [], []
    for _ in range(150):
        Z, s, ok = standardize(gaussian_panel(Sig, n, rng, mu, sd))
        v = w * s
        m = panel_moments(Z, v, ok)
        est.append(m["Q"])
        raw.append(m["Q_raw"])
        tru.append(v @ Sig @ Sig @ v)
    rel = np.mean(est) / np.mean(tru) - 1
    assert abs(rel) < 0.04
    assert np.mean(raw) / np.mean(tru) - 1 > abs(rel)


def test_floor_is_cauchy_schwarz_bound():
    rng = np.random.default_rng(2)
    Z, s, ok = standardize(rng.normal(size=(30, 200)))
    v = rng.normal(size=200)
    m = panel_moments(Z, v, ok)
    assert m["Q"] >= m["L"] ** 2 / (v @ v) - 1e-12
    assert m["Q"] == max(m["Q_untruncated"], m["L"] ** 2 / (v @ v))


def test_jackknife_equals_recomputation_from_raw_genotypes_with_fixed_v():
    rng = np.random.default_rng(3)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, 60), size=(200, 60)).astype(float)
    G[:, 5] = 1.0
    v = rng.normal(size=60) * (rng.random(60) < 0.2)
    Z, s, ok = standardize(G)
    blocks = np.array_split(rng.permutation(200), 20)
    L, Q = jackknife_moments(Z, v, blocks=blocks)
    for k, b in enumerate(blocks):
        keep = np.setdiff1d(np.arange(200), b)
        Zk, sk, okk = standardize(G[keep])
        m = panel_moments(Zk, v, okk)
        assert np.isclose(m["L"], L[k], rtol=1e-10)
        assert np.isclose(m["Q"], Q[k], rtol=1e-10)


def test_reaudit_one_snp_counterexample_gives_normal_reference():
    rng = np.random.default_rng(31)
    for _ in range(20):
        G = rng.binomial(2, 0.3, size=(500, 1)).astype(float)
        r = xcal_gene(np.array([1.0]), np.array([3.0 * np.sqrt(1.4 * 0.42)]), G, c=0.4, snp_sd=sd_from_eaf([0.3]))
        assert r["df"] >= 1e5 - 1
        assert np.isclose(r["p_xcal_t"], r["p_xcal_normal"], rtol=1e-3)


def test_parameter_uncertainty_lowers_df():
    rng = np.random.default_rng(32)
    G = gaussian_panel(ar1(50, 0.8), 400, rng)
    w = np.zeros(50); w[[5, 20]] = [0.4, -0.3]
    z = rng.normal(size=50)
    r0 = xcal_gene(w, z, G, c=0.3)
    r1 = xcal_gene(w, z, G, c=0.3, se_c=0.05, se_a=0.02)
    assert r1["df"] < r0["df"] and r1["p_xcal_t"] >= r0["p_xcal_t"]


@pytest.mark.parametrize("kw", [dict(c=-1.0), dict(c=float("nan")), dict(c=0.1, intercept=-1.0),
                                dict(c=0.1, intercept=float("inf")), dict(c=0.1, K=1), dict(c=0.1, se_c=-0.1)])
def test_invalid_parameters_rejected(kw):
    rng = np.random.default_rng(33)
    G = rng.binomial(2, 0.3, size=(100, 5)).astype(float)
    with pytest.raises(ValueError):
        xcal_gene(np.ones(5), np.ones(5), G, **kw)


def test_invalid_vectors_rejected():
    rng = np.random.default_rng(34)
    G = rng.binomial(2, 0.3, size=(100, 5)).astype(float)
    with pytest.raises(ValueError):
        xcal_gene(np.ones((5, 1)), np.ones(5), G, c=0.1)
    with pytest.raises(ValueError):
        xcal_gene(np.r_[np.ones(4), np.nan], np.ones(5), G, c=0.1)
    with pytest.raises(ValueError):
        xcal_gene(np.ones(5), np.ones(5), G[:10], c=0.1)


def test_xcal_gene_outputs_and_t_is_more_conservative():
    rng = np.random.default_rng(4)
    G = gaussian_panel(ar1(80, 0.8), 400, rng)
    w = np.zeros(80); w[[10, 30]] = [0.5, -0.3]
    z = rng.normal(size=80) * 3
    r = xcal_gene(w, z, G, c=0.3)
    assert r["df"] >= 4 and r["p_xcal_t"] >= r["p_xcal_normal"]
    assert r["inflation_factor"] > 1


def test_xcal_t_calibrated_in_small_simulation():
    rng = np.random.default_rng(5)
    p, c = 60, 0.5
    Sig = ar1(p, 0.85)
    Lg = np.linalg.cholesky(Sig + c * Sig @ Sig)
    rej, n = 0, 0
    for g in range(40):
        G = gaussian_panel(Sig, 500, rng)
        w = np.zeros(p); idx = rng.choice(p, 3, replace=False); w[idx] = rng.normal(size=3)
        Z, s, ok = standardize(G)
        v = w * s
        m = panel_moments(Z, v, ok)
        Lb, Qb = jackknife_moments(Z, v)
        zs = (rng.normal(size=(2000, p)) @ Lg.T) @ v
        r = xcal_from_moments(zs, m["L"], m["Q"], c, 1.0, Lb, Qb)
        rej += np.sum(r["p_xcal_t"] < 0.01); n += len(zs)
    assert rej / n < 0.0125


def _toy_inputs(rng):
    snps = np.array([f"rs{i}" for i in range(50)])
    pos = np.arange(50) * 10_000
    G = rng.binomial(2, 0.3, size=(300, 50)).astype(float)
    weights = pd.DataFrame({"gene": ["A"] * 3 + ["B"] * 2, "snp": ["rs3", "rs4", "rs99", "rs20", "rs21"],
                            "weight": [0.2, -0.1, 0.4, 0.3, 0.1]})
    gwas = pd.DataFrame({"snp": snps, "z": rng.normal(size=50)})
    return weights, gwas, G, snps, pos


@pytest.mark.filterwarnings("ignore:gwas has no")
def test_run_xcal_reports_coverage_and_rejects_duplicates():
    rng = np.random.default_rng(6)
    weights, gwas, G, snps, pos = _toy_inputs(rng)
    res = run_xcal(weights, gwas, G, snps, pos, c=0.1, window_bp=100_000)
    a = res.set_index("gene").loc["A"]
    assert a.n_snps_used == 2 and a.n_not_in_panel == 1 and a.weight_coverage < 1
    assert set(res.status) == {"ok"}
    res2 = run_xcal(weights, gwas.assign(eaf=0.3), G, snps, pos, c=0.1, window_bp=100_000)
    assert set(res2.sd_source) == {"gwas"}
    with pytest.raises(ValueError):
        run_xcal(weights, gwas, G, snps, pos, c=-0.1)
    with pytest.raises(ValueError):
        run_xcal(pd.concat([weights, weights.iloc[:1]]), gwas, G, snps, pos, c=0.1)
    with pytest.raises(ValueError):
        run_xcal(weights, pd.concat([gwas, gwas.iloc[:1]]), G, snps, pos, c=0.1)
    G2 = G.copy(); G2[0, 0] = np.nan
    with pytest.raises(ValueError):
        run_xcal(weights, gwas, G2, snps, pos, c=0.1)


def test_estimate_c_recovers_slope_with_ld():
    rng = np.random.default_rng(7)
    c_true, blocks, p = 0.05, 300, 40
    zs, ls = [], []
    for b in range(blocks):
        Sig = ar1(p, rng.uniform(0.6, 0.95))
        Gp = gaussian_panel(Sig, 400, rng)
        ls.append(ld_scores(Gp, np.arange(p), window_bp=p))
        C = Sig + c_true * Sig @ Sig
        zs.append(np.linalg.cholesky(C) @ rng.normal(size=p))
    r = estimate_c(np.concatenate(zs), np.concatenate(ls))
    assert abs(r["c"] - c_true) < 3 * r["se_c"] + 0.01
    assert abs(r["intercept"] - 1) < 3 * r["se_intercept"] + 0.05
    assert r["se_c"] > 0


def test_genomic_control_and_harmonize():
    rng = np.random.default_rng(8)
    _, lam = genomic_control(rng.normal(size=200_000))
    assert abs(lam - 1) < 0.02
    df = pd.DataFrame({"snp": ["a", "b", "c", "d"], "effect_allele": ["A", "G", "A", "C"],
                       "other_allele": ["G", "A", "T", "T"], "weight": [1.0, 1.0, 1.0, 1.0]})
    panel = pd.DataFrame({"snp": ["a", "b", "c"], "a1": ["A", "A", "A"], "a2": ["G", "G", "T"]})
    out, rep = harmonize(df, panel)
    assert out.set_index("snp").weight.to_dict() == {"a": 1.0, "b": -1.0}
    assert rep == {"aligned": 1, "flipped": 1, "ambiguous": 1, "not_in_panel": 1}
    assert out.set_index("snp").loc["b", ["effect_allele", "other_allele"]].tolist() == ["A", "G"]
    assert out.set_index("snp").loc["b", "original_effect_allele"] == "G"
    again, rep2 = harmonize(out, panel)
    assert again.set_index("snp").weight.to_dict() == {"a": 1.0, "b": -1.0} and rep2 == {"aligned": 2}
    lower = panel.assign(a1=panel.a1.str.lower(), a2=panel.a2.str.lower())
    assert harmonize(df, lower)[0].weight.tolist() == out.weight.tolist()
    with pytest.raises(ValueError):
        harmonize(df, pd.concat([panel, panel.iloc[:1]]))


def test_estimate_c_rejects_degenerate_input():
    with pytest.raises(ValueError):
        estimate_c(np.ones(500), np.ones(500))
    with pytest.raises(ValueError):
        estimate_c(np.random.default_rng(0).normal(size=5000), np.full(5000, 3.0))


def test_snp_monomorphic_in_panel_is_dropped_from_numerator_with_gwas_sds():
    rng = np.random.default_rng(11)
    G = rng.binomial(2, 0.3, size=(300, 4)).astype(float)
    G[:, 2] = 0.0
    w = np.array([0.2, -0.1, 5.0, 0.3])
    z = np.array([1.0, -0.5, 40.0, 0.2])
    sd = np.full(4, 0.6)
    full = xcal_gene(w, z, G, c=0.1, snp_sd=sd)
    ref = xcal_gene(np.delete(w, 2), np.delete(z, 2), np.delete(G, 2, axis=1), c=0.1, snp_sd=np.delete(sd, 2))
    assert np.isclose(full["z_xcal"], ref["z_xcal"])
    assert np.isclose(full["df"], ref["df"])
