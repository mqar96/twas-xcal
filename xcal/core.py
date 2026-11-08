import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


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
    Z = np.asarray(Z, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    n = Z.shape[0]
    p = int(np.count_nonzero(ok)) if ok is not None else int(np.count_nonzero(Z.std(0) > 0))
    u = Z @ v
    lin = float(u @ u) / n
    Rv = Z.T @ u / n
    raw = float(Rv @ Rv)
    q_unb = (raw - (p / n) * lin) / (1.0 + 1.0 / n)
    vv = float(v @ v)
    floor = lin * lin / vv if vv > 0 else 0.0
    Q = max(q_unb, floor)
    return dict(L=lin, Q=Q, Q_untruncated=q_unb, Q_raw=raw, floor_active=bool(q_unb < floor), n=n, p=p)


def jackknife_moments(Z, s, w, K=20, seed=0, blocks=None):
    Z = np.asarray(Z, dtype=np.float64)
    n = Z.shape[0]
    if blocks is None:
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
        vk = np.where(okk, w * s * sk, 0.0)
        mom = panel_moments(Zs, vk, okk)
        L.append(mom["L"])
        Q.append(mom["Q"])
    return np.array(L), np.array(Q)


def satterthwaite_df(D_blocks, D_full, lo=4.0, hi=1e5):
    D_blocks = np.asarray(D_blocks, float)
    K = D_blocks.shape[-1]
    var = (K - 1) / K * np.sum((D_blocks - D_blocks.mean(-1, keepdims=True)) ** 2, axis=-1)
    tau2 = var / np.asarray(D_full, float) ** 2
    return np.clip(2.0 / np.maximum(tau2, 1e-12), lo, hi)


def xcal_from_moments(num, L, Q, c, intercept=1.0, L_blocks=None, Q_blocks=None):
    D = intercept * L + c * Q
    zx = num / np.sqrt(D)
    out = dict(z_xcal=zx, p_xcal_normal=2 * norm.sf(np.abs(zx)), D=D)
    if L_blocks is not None:
        df = satterthwaite_df(intercept * np.asarray(L_blocks) + c * np.asarray(Q_blocks), D)
        out["df"] = df
        out["p_xcal_t"] = 2 * tdist.sf(np.abs(zx), df)
    return out


def xcal_gene(w, z, G, c, intercept=1.0, K=20, seed=0):
    w = np.asarray(w, float)
    z = np.asarray(z, float)
    Z, s, ok = standardize(G)
    if w.shape != z.shape or w.shape[0] != Z.shape[1]:
        raise ValueError("w, z and the genotype window must have matching SNP dimension")
    v = np.where(ok, w * s, 0.0)
    mom = panel_moments(Z, v, ok)
    if mom["L"] <= 0:
        return None
    num = float(v @ z)
    Lb, Qb = jackknife_moments(Z, s, np.where(ok, w, 0.0), K=K, seed=seed)
    res = xcal_from_moments(num, mom["L"], mom["Q"], c, intercept, Lb, Qb)
    return dict(z_standard=num / np.sqrt(mom["L"]), L=mom["L"], Q=mom["Q"], q=mom["Q"] / mom["L"],
                floor_active=mom["floor_active"], inflation_factor=res["D"] / mom["L"],
                z_xcal=float(res["z_xcal"]), df=float(res["df"]), p_xcal_t=float(res["p_xcal_t"]),
                p_xcal_normal=float(res["p_xcal_normal"]), n_panel=mom["n"])


def ld_scores(G, pos=None, window_bp=1_000_000):
    Z, s, ok = standardize(G)
    n, p = Z.shape
    pos = np.arange(p) if pos is None else np.asarray(pos)
    order = np.argsort(pos, kind="stable")
    if not np.all(order == np.arange(p)):
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
    z = np.asarray(z, float)
    ld = np.asarray(ld, float)
    ok = np.isfinite(z) & np.isfinite(ld)
    chi2 = z[ok] ** 2
    l = ld[ok]
    m = len(chi2)
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
    if n_blocks == 0:
        return dict(c=max(c, 0.0), c_raw=c, se_c=float("nan"), intercept=a, se_intercept=float("nan"), n_snps=int(m), n_blocks=0)
    edges = np.linspace(0, m, n_blocks + 1).astype(int)
    jk = []
    for i in range(n_blocks):
        keep = np.ones(m, bool)
        keep[edges[i]:edges[i + 1]] = False
        jk.append(fit(chi2[keep], l[keep]))
    jk = np.array(jk)
    se = np.sqrt((n_blocks - 1) / n_blocks * ((jk - jk.mean(0)) ** 2).sum(0))
    return dict(c=max(c, 0.0), c_raw=c, se_c=float(se[1]), intercept=a,
                se_intercept=float(se[0]) if fix_intercept is None else 0.0, n_snps=int(m), n_blocks=int(n_blocks))


def genomic_control(z):
    z = np.asarray(z, float)
    lam = float(np.nanmedian(z ** 2) / 0.454936423119572)
    return z / np.sqrt(max(lam, 1.0)), lam


def harmonize(df, snp_alleles, effect_col="effect_allele", other_col="other_allele", value_col="weight"):
    """Align a per-SNP value to the panel's counted allele (a1). Drops strand-ambiguous and mismatched SNPs."""
    m = df.merge(snp_alleles[["snp", "a1", "a2"]], on="snp", how="left")
    e = m[effect_col].str.upper()
    o = m[other_col].str.upper()
    ce = e.map(COMPLEMENT)
    co = o.map(COMPLEMENT)
    same = ((e == m.a1) & (o == m.a2)) | ((ce == m.a1) & (co == m.a2))
    flip = ((e == m.a2) & (o == m.a1)) | ((ce == m.a2) & (co == m.a1))
    ambiguous = (e.map(COMPLEMENT) == o)
    status = np.select([m.a1.isna(), ambiguous, same, flip], ["not_in_panel", "ambiguous", "aligned", "flipped"], "mismatch")
    m["harmonize_status"] = status
    m[value_col] = np.where(status == "flipped", -m[value_col], m[value_col])
    keep = np.isin(status, ["aligned", "flipped"])
    return m[keep].drop(columns=["a1", "a2"]).reset_index(drop=True), m["harmonize_status"].value_counts().to_dict()


def _check_unique(series, what):
    d = series[series.duplicated()].unique()
    if len(d):
        raise ValueError(f"{len(d)} duplicated {what} (e.g. {list(d[:3])}); resolve duplicates before calling")


def run_xcal(weights, gwas, panel_genotypes, panel_snps, panel_pos, c, intercept=1.0, window_bp=1_000_000,
             K=20, seed=0, min_coverage=0.0):
    """XCAL / XCAL-t for every gene on one chromosome.

    weights: DataFrame[gene, snp, weight]; weights on the allele-dosage scale of the panel's counted allele
    gwas: DataFrame[snp, z] aligned to the same allele
    panel_genotypes: n x p dosage matrix (no missing values) from a panel of the GWAS ancestry
    panel_snps, panel_pos: SNP ids and positions of the panel columns
    c: polygenic scale N*h2/M (estimate_c); intercept: LDSC intercept (1.0 = no confounding)
    Returns one row per gene; p_xcal_t is the recommended p-value.
    """
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
    if not np.all(np.isfinite(weights["weight"])) or not np.all(np.isfinite(gwas["z"])):
        raise ValueError("non-finite weights or z-scores")
    col = pd.Series(np.arange(len(panel_snps)), index=panel_snps)
    zmap = pd.Series(gwas["z"].values, index=gwas["snp"].values)
    out = []
    for gene, g in weights.groupby("gene", sort=False):
        total_w2 = float((g.weight ** 2).sum())
        in_panel = g.snp.isin(col.index)
        in_gwas = g.snp.isin(zmap.index)
        use = g[in_panel & in_gwas]
        row = dict(gene=gene, n_snps_model=len(g), n_not_in_panel=int((~in_panel).sum()),
                   n_not_in_gwas=int((in_panel & ~in_gwas).sum()), n_snps_used=len(use),
                   weight_coverage=float((use.weight ** 2).sum() / total_w2) if total_w2 > 0 else 0.0)
        if use.empty or row["weight_coverage"] < min_coverage:
            out.append(dict(row, status="insufficient_coverage"))
            continue
        widx = col[use.snp].values
        lo, hi = panel_pos[widx].min() - window_bp, panel_pos[widx].max() + window_bp
        win = np.where((panel_pos >= lo) & (panel_pos <= hi))[0]
        loc = np.searchsorted(win, widx)
        w = np.zeros(len(win))
        z = np.zeros(len(win))
        w[loc] = use.weight.values
        z[loc] = zmap[use.snp].values
        r = xcal_gene(w, z, G[:, win], c, intercept, K, seed)
        if r is None:
            out.append(dict(row, n_window_snps=len(win), status="zero_variance"))
            continue
        out.append(dict(row, n_window_snps=len(win), status="ok", **r))
    return pd.DataFrame(out)
