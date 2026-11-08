import numbers
import warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
DF_LO, DF_HI = 4.0, 1e5
MIN_PANEL = 20


def _check_params(c, intercept, se_c=0.0, se_a=0.0):
    for name, x in (("c", c), ("intercept", intercept), ("se_c", se_c), ("se_a", se_a)):
        if isinstance(x, (bool, np.bool_)) or not isinstance(x, numbers.Real) or not np.isfinite(x):
            raise ValueError(f"{name} must be a finite real number (got {x!r})")
    if c < 0:
        raise ValueError(f"c must be >= 0 (got {c})")
    if intercept <= 0:
        raise ValueError(f"intercept must be > 0 (got {intercept})")
    if se_c < 0 or se_a < 0:
        raise ValueError("standard errors must be >= 0")


def _check_K(K, n):
    if not isinstance(K, numbers.Integral) or K < 2:
        raise ValueError(f"K must be an integer >= 2 (got {K!r})")
    if n < MIN_PANEL or n < 2 * K:
        raise ValueError(f"panel too small for the jackknife: n = {n}, need n >= max({MIN_PANEL}, 2K)")


def sd_from_eaf(eaf):
    f = np.asarray(eaf, float)
    if np.any(~np.isfinite(f)) or np.any((f < 0) | (f > 1)):
        raise ValueError("allele frequencies must be finite and in [0, 1]")
    return np.sqrt(2.0 * f * (1.0 - f))


def standardize(G):
    X = np.asarray(G, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("genotype matrix must be 2-D (individuals x SNPs)")
    if not np.all(np.isfinite(X)):
        raise ValueError("genotype matrix contains missing or non-finite values; impute before calling")
    mu = X.mean(0)
    sd = X.std(0)
    ok = sd > 0
    Z = np.zeros_like(X)
    Z[:, ok] = (X[:, ok] - mu[ok]) / sd[ok]
    return Z, sd, ok


def panel_moments(Z, v, ok=None):
    """L = v'Rv and Q (estimate of v'Sigma^2 v) from a standardized panel Z (n x p) and a fixed vector v."""
    Z = np.asarray(Z, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    n = Z.shape[0]
    ok = (Z.std(0) > 0) if ok is None else np.asarray(ok, bool)
    v = np.where(ok, v, 0.0)
    p = int(np.count_nonzero(ok))
    u = Z @ v
    lin = float(u @ u) / n
    Rv = Z.T @ u / n
    raw = float(Rv @ Rv)
    q_unb = (raw - (p / n) * lin) / (1.0 + 1.0 / n)
    vv = float(v @ v)
    floor = lin * lin / vv if vv > 0 else 0.0
    Q = max(q_unb, floor)
    return dict(L=lin, Q=Q, Q_untruncated=q_unb, Q_raw=raw, floor_active=bool(q_unb < floor), n=n, p=p)


def jackknife_moments(Z, v, K=20, seed=0, blocks=None):
    """Delete-group jackknife replicates of (L, Q) with v held fixed.

    Each replicate re-centres and re-scales the retained individuals (i.e. recomputes the sample
    correlation matrix from their genotypes) while keeping v fixed, so the replicates describe the
    panel-sampling variability of D-hat = a L + c Q for the fixed target D = a v'Sigma v + c v'Sigma^2 v.
    """
    Z = np.asarray(Z, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    n = Z.shape[0]
    if blocks is None:
        _check_K(K, n)
        rng = np.random.default_rng(seed)
        blocks = np.array_split(rng.permutation(n), K)
    L, Q = [], []
    for b in blocks:
        keep = np.ones(n, bool)
        keep[b] = False
        Zk = Z[keep]
        m = Zk.mean(0)
        sk = Zk.std(0)
        okk = sk > 1e-12
        Zs = np.zeros_like(Zk)
        Zs[:, okk] = (Zk[:, okk] - m[okk]) / sk[okk]
        mom = panel_moments(Zs, v, okk)
        L.append(mom["L"])
        Q.append(mom["Q"])
    return np.array(L), np.array(Q)


def jackknife_var(x_blocks):
    x = np.asarray(x_blocks, float)
    K = x.shape[-1]
    return (K - 1) / K * np.sum((x - x.mean(-1, keepdims=True)) ** 2, axis=-1)


def satterthwaite_df(D_blocks, D_full, extra_var=0.0, lo=DF_LO, hi=DF_HI):
    """nu = 2 / tau^2 with tau^2 = [Var_jk(D-hat) + extra_var] / D-hat^2, clipped to [lo, hi]."""
    tau2 = (jackknife_var(D_blocks) + np.asarray(extra_var, float)) / np.asarray(D_full, float) ** 2
    return np.clip(2.0 / np.maximum(tau2, 1e-300), lo, hi)


def xcal_from_moments(num, L, Q, c, intercept=1.0, L_blocks=None, Q_blocks=None, se_c=0.0, se_a=0.0):
    """XCAL statistic and (if jackknife replicates are given) XCAL-t p-value.

    se_c / se_a: standard errors of c and of the intercept a (e.g. from estimate_c); their contribution
    Q^2 se_c^2 + L^2 se_a^2 is added to the denominator variance used for the degrees of freedom.
    """
    _check_params(c, intercept, se_c, se_a)
    L = np.asarray(L, float)
    Q = np.asarray(Q, float)
    D = intercept * L + c * Q
    if np.any(~(D > 0)):
        raise ValueError("non-positive denominator; check the panel and weights")
    num = np.asarray(num, float)
    if not np.all(np.isfinite(num)):
        raise ValueError("numerator must be finite")
    if L_blocks is not None and (np.asarray(L_blocks).shape[-1] < 2 or np.asarray(Q_blocks).shape != np.asarray(L_blocks).shape):
        raise ValueError("jackknife replicates must have at least two blocks and matching L/Q shapes")
    zx = num / np.sqrt(D)
    out = dict(z_xcal=zx, p_xcal_normal=2 * norm.sf(np.abs(zx)), D=D)
    if L_blocks is not None:
        Db = intercept * np.asarray(L_blocks, float) + c * np.asarray(Q_blocks, float)
        df = satterthwaite_df(Db, D, extra_var=(Q * se_c) ** 2 + (L * se_a) ** 2)
        out["df"] = df
        out["df_at_lower_clip"] = df <= DF_LO
        out["p_xcal_t"] = 2 * tdist.sf(np.abs(zx), df)
    return out


def xcal_gene(w, z, G, c, intercept=1.0, snp_sd=None, K=20, seed=0, se_c=0.0, se_a=0.0):
    """XCAL-t for one gene.

    w: weights on the allele-dosage scale; z: GWAS z-scores (same SNP order and allele);
    G: n x p dosage matrix of an LD panel of the GWAS ancestry;
    snp_sd: per-SNP genotype SD in the GWAS (e.g. sd_from_eaf(eaf)); recommended. If None, panel SDs are used
    (v then depends on the panel and the jackknife is conditional on the panel SDs).
    """
    _check_params(c, intercept, se_c, se_a)
    w = np.asarray(w, float)
    z = np.asarray(z, float)
    if w.ndim != 1 or z.ndim != 1:
        raise ValueError("w and z must be one-dimensional")
    if not (np.all(np.isfinite(w)) and np.all(np.isfinite(z))):
        raise ValueError("w and z must be finite")
    Z, s, ok = standardize(G)
    if w.shape != z.shape or w.shape[0] != Z.shape[1]:
        raise ValueError("w, z and the genotype window must have matching SNP dimension")
    _check_K(K, Z.shape[0])
    if snp_sd is None:
        sd, source = s, "panel"
    else:
        sd = np.asarray(snp_sd, float)
        if sd.shape != w.shape or np.any(~np.isfinite(sd)) or np.any(sd < 0):
            raise ValueError("snp_sd must be finite, non-negative and match w")
        source = "gwas"
    v = np.where(ok, w * sd, 0.0)
    mom = panel_moments(Z, v, ok)
    if mom["L"] <= 0:
        return None
    num = float(v @ z)
    Lb, Qb = jackknife_moments(Z, v, K=K, seed=seed)
    res = xcal_from_moments(num, mom["L"], mom["Q"], c, intercept, Lb, Qb, se_c, se_a)
    return dict(z_standard=num / np.sqrt(mom["L"]), L=mom["L"], Q=mom["Q"], q=mom["Q"] / mom["L"],
                floor_active=mom["floor_active"], inflation_factor=float(res["D"] / mom["L"]),
                z_xcal=float(res["z_xcal"]), df=float(res["df"]), df_at_lower_clip=bool(res["df_at_lower_clip"]),
                p_xcal_t=float(res["p_xcal_t"]), p_xcal_normal=float(res["p_xcal_normal"]),
                n_panel=mom["n"], sd_source=source)


def ld_scores(G, pos=None, window_bp=1_000_000):
    Z, s, ok = standardize(G)
    n, p = Z.shape
    if n < 3:
        raise ValueError("need at least 3 individuals for bias-corrected LD scores")
    pos = np.arange(p) if pos is None else np.asarray(pos)
    if np.any(np.diff(pos) < 0):
        raise ValueError("positions must be sorted")
    l = np.zeros(p)
    step = 2000
    for a in range(0, p, step):
        b = min(p, a + step)
        lo = np.searchsorted(pos, pos[a] - window_bp)
        hi = np.searchsorted(pos, pos[b - 1] + window_bp, side="right")
        r = Z[:, a:b].T @ Z[:, lo:hi] / n
        r2 = r ** 2
        r2 = r2 - (1 - r2) / (n - 2)
        dist = np.abs(pos[a:b][:, None] - pos[lo:hi][None, :])
        r2[dist > window_bp] = 0.0
        r2[:, ~ok[lo:hi]] = 0.0
        l[a:b] = r2.sum(1)
    l[~ok] = np.nan
    return l


def estimate_c(z, ld, n_blocks=200, fix_intercept=None):
    """LD-score regression of z^2 on LD scores: returns c (slope), intercept a, block-jackknife SEs and cov(a, c)."""
    z = np.asarray(z, float)
    ld = np.asarray(ld, float)
    if z.shape != ld.shape:
        raise ValueError("z and ld must have the same length")
    ok = np.isfinite(z) & np.isfinite(ld)
    chi2 = z[ok] ** 2
    l = ld[ok]
    m = len(chi2)
    if m < 1000:
        raise ValueError(f"too few SNPs for LD-score regression ({m}); need >= 1000")
    if np.std(l) < 1e-8 * max(1.0, abs(np.mean(l))):
        raise ValueError("LD scores have no variation; slope is not identified")
    if fix_intercept is not None:
        _check_params(0.0, fix_intercept)
    if n_blocks and m < 10 * n_blocks:
        n_blocks = max(10, m // 10)

    def fit(chi2, l):
        a, c = 1.0, max((chi2.mean() - 1) / max(l.mean(), 1e-8), 0.0)
        for _ in range(3):
            wt = 1.0 / (np.maximum(l, 1.0) * np.maximum(a + c * np.maximum(l, 0), 1e-3) ** 2)
            sw = np.sqrt(wt)
            if fix_intercept is None:
                X = np.column_stack([np.ones_like(l), l])
                coef, *_ = np.linalg.lstsq(X * sw[:, None], chi2 * sw, rcond=None)
                a, c = float(coef[0]), float(coef[1])
            else:
                a = float(fix_intercept)
                c = float(np.sum(wt * l * (chi2 - a)) / np.sum(wt * l * l))
        return a, c

    a, c = fit(chi2, l)
    if c < 0:
        warnings.warn(f"estimated slope is negative ({c:.3g}); returning c = 0")
    base = dict(c=max(c, 0.0), c_raw=c, intercept=a, n_snps=int(m), fixed_intercept=fix_intercept is not None)
    if n_blocks == 0:
        return dict(base, se_c=float("nan"), se_intercept=float("nan"), cov_intercept_c=float("nan"), n_blocks=0)
    edges = np.linspace(0, m, n_blocks + 1).astype(int)
    jk = []
    for i in range(n_blocks):
        keep = np.ones(m, bool)
        keep[edges[i]:edges[i + 1]] = False
        jk.append(fit(chi2[keep], l[keep]))
    jk = np.array(jk)
    dev = jk - jk.mean(0)
    cov = (n_blocks - 1) / n_blocks * dev.T @ dev
    return dict(base, se_c=float(np.sqrt(cov[1, 1])), se_intercept=float(np.sqrt(cov[0, 0])) if fix_intercept is None else 0.0,
                cov_intercept_c=float(cov[0, 1]) if fix_intercept is None else 0.0, n_blocks=int(n_blocks))


def genomic_control(z):
    z = np.asarray(z, float)
    lam = float(np.nanmedian(z ** 2) / 0.454936423119572)
    return z / np.sqrt(max(lam, 1.0)), lam


def harmonize(df, snp_alleles, effect_col="effect_allele", other_col="other_allele", value_col="weight"):
    """Align per-SNP values to the panel's counted allele a1.

    Returns a table whose effect/other allele columns are the panel's a1/a2 (original labels kept in
    original_<col>), with the value sign-flipped where needed. Strand-ambiguous (A/T, C/G), mismatched and
    absent SNPs are dropped. Repeated calls are stable.
    """
    pa = snp_alleles[["snp", "a1", "a2"]].copy()
    pa["a1"] = pa["a1"].astype(str).str.upper()
    pa["a2"] = pa["a2"].astype(str).str.upper()
    dup = pa.snp[pa.snp.duplicated()].unique()
    if len(dup):
        raise ValueError(f"{len(dup)} duplicated SNPs in snp_alleles (e.g. {list(dup[:3])})")
    m = df.drop(columns=[c for c in ("a1", "a2") if c in df]).merge(pa, on="snp", how="left")
    e = m[effect_col].astype(str).str.upper()
    o = m[other_col].astype(str).str.upper()
    ce, co = e.map(COMPLEMENT), o.map(COMPLEMENT)
    same = ((e == m.a1) & (o == m.a2)) | ((ce == m.a1) & (co == m.a2))
    flip = ((e == m.a2) & (o == m.a1)) | ((ce == m.a2) & (co == m.a1))
    ambiguous = ce == o
    status = np.select([m.a1.isna(), ambiguous, same, flip], ["not_in_panel", "ambiguous", "aligned", "flipped"], "mismatch")
    m["harmonize_status"] = status
    m[value_col] = np.where(status == "flipped", -m[value_col], m[value_col])
    for col in (effect_col, other_col):
        if f"original_{col}" not in m:
            m[f"original_{col}"] = m[col]
    keep = np.isin(status, ["aligned", "flipped"])
    m.loc[keep, effect_col] = m.loc[keep, "a1"]
    m.loc[keep, other_col] = m.loc[keep, "a2"]
    report = pd.Series(status).value_counts().to_dict()
    return m[keep].drop(columns=["a1", "a2"]).reset_index(drop=True), report


def _check_unique(series, what):
    d = series[series.duplicated()].unique()
    if len(d):
        raise ValueError(f"{len(d)} duplicated {what} (e.g. {list(d[:3])}); resolve duplicates before calling")


def run_xcal(weights, gwas, panel_genotypes, panel_snps, panel_pos, c, intercept=1.0, window_bp=1_000_000,
             K=20, seed=0, min_coverage=0.0, se_c=0.0, se_a=0.0):
    """XCAL-t for every gene on one chromosome.

    weights: DataFrame[gene, snp, weight] aligned to the panel's counted allele (see harmonize)
    gwas: DataFrame[snp, z] plus, recommended, 'eaf' (GWAS frequency of the counted allele) or 'sd'
    panel_genotypes: n x p dosage matrix (no missing values) from a panel of the GWAS ancestry
    c, intercept (and optionally se_c, se_a): from estimate_c on the same GWAS
    Returns one row per gene; p_xcal_t is the recommended p-value.
    """
    _check_params(c, intercept, se_c, se_a)
    if not (0.0 <= min_coverage <= 1.0):
        raise ValueError("min_coverage must be in [0, 1]")
    for col in ("gene", "snp", "weight"):
        if col not in weights:
            raise ValueError(f"weights is missing column '{col}'")
    for col in ("snp", "z"):
        if col not in gwas:
            raise ValueError(f"gwas is missing column '{col}'")
    _check_unique(weights["gene"].astype(str) + ":" + weights["snp"].astype(str), "gene:snp pairs in weights")
    _check_unique(gwas["snp"], "SNPs in gwas")
    panel_snps = np.asarray(panel_snps)
    panel_pos = np.asarray(panel_pos)
    _check_unique(pd.Series(panel_snps), "SNPs in panel")
    G = np.asarray(panel_genotypes)
    if G.ndim != 2 or G.shape[1] != len(panel_snps) or len(panel_pos) != len(panel_snps):
        raise ValueError("panel_genotypes must be n x p with p = len(panel_snps) = len(panel_pos)")
    _check_K(K, G.shape[0])
    if not np.all(np.isfinite(weights["weight"])) or not np.all(np.isfinite(gwas["z"])):
        raise ValueError("non-finite weights or z-scores")
    if "sd" in gwas:
        gsd = gwas["sd"].to_numpy(float)
    elif "eaf" in gwas:
        gsd = sd_from_eaf(gwas["eaf"].to_numpy(float))
    else:
        gsd = None
        warnings.warn("gwas has no 'eaf' or 'sd' column; using panel SDs (XCAL-t df then conditional on the panel SDs)")
    col = pd.Series(np.arange(len(panel_snps)), index=panel_snps)
    zmap = pd.Series(gwas["z"].values, index=gwas["snp"].values)
    sdmap = pd.Series(gsd, index=gwas["snp"].values) if gsd is not None else None
    out = []
    for gene, g in weights.groupby("gene", sort=False):
        total_w2 = float((g.weight ** 2).sum())
        in_panel = g.snp.isin(col.index)
        in_gwas = g.snp.isin(zmap.index)
        use = g[in_panel & in_gwas]
        row = dict(gene=gene, n_snps_model=len(g), n_not_in_panel=int((~in_panel).sum()),
                   n_not_in_gwas=int((in_panel & ~in_gwas).sum()), n_snps_used=len(use),
                   weight_coverage=float((use.weight ** 2).sum() / total_w2) if total_w2 > 0 else 0.0)
        if use.empty or row["weight_coverage"] < min_coverage or row["weight_coverage"] == 0:
            out.append(dict(row, status="insufficient_coverage"))
            continue
        widx = col[use.snp].values
        lo, hi = panel_pos[widx].min() - window_bp, panel_pos[widx].max() + window_bp
        win = np.where((panel_pos >= lo) & (panel_pos <= hi))[0]
        loc = np.searchsorted(win, widx)
        w = np.zeros(len(win)); z = np.zeros(len(win)); sdw = np.zeros(len(win))
        w[loc] = use.weight.values
        z[loc] = zmap[use.snp].values
        if sdmap is not None:
            sdw[loc] = sdmap[use.snp].values
        r = xcal_gene(w, z, G[:, win], c, intercept, sdw if sdmap is not None else None, K, seed, se_c, se_a)
        if r is None:
            out.append(dict(row, n_window_snps=len(win), status="zero_variance"))
            continue
        out.append(dict(row, n_window_snps=len(win), status="ok", **r))
    return pd.DataFrame(out)
