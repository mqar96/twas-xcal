import numpy as np
import pandas as pd
from scipy.stats import norm


def ld_from_genotypes(G):
    G = np.asarray(G, dtype=np.float64)
    s = G.std(0)
    ok = s > 0
    Z = np.zeros_like(G)
    Z[:, ok] = (G[:, ok] - G[:, ok].mean(0)) / s[ok]
    R = Z.T @ Z / G.shape[0]
    R[~ok, :] = 0
    R[:, ~ok] = 0
    R[~ok, ~ok] = 1.0
    return R, s, G.shape[0]


def unbiased_quad_R2(R, n, v):
    p = R.shape[0]
    Rv = R @ v
    lin = float(v @ Rv)
    est = (float(Rv @ Rv) - (p / n) * lin) / (1.0 + 1.0 / n)
    return max(est, lin)


def xcal_gene(w, z, R, s, n, c):
    w = np.asarray(w, float)
    z = np.asarray(z, float)
    v = w * s
    lin = float(v @ R @ v)
    if lin <= 0:
        return np.nan, np.nan, np.nan, np.nan
    q = unbiased_quad_R2(R, n, v) / lin
    z_std = float(v @ z) / np.sqrt(lin)
    z_x = z_std / np.sqrt(1.0 + c * q)
    return z_x, 2 * norm.sf(abs(z_x)), z_std, q


def estimate_c(z, ld_scores, n_snps_total=None):
    z = np.asarray(z, float)
    l = np.asarray(ld_scores, float)
    ok = np.isfinite(z) & np.isfinite(l)
    chi2 = z[ok] ** 2
    l = l[ok]
    w = 1.0 / np.maximum(l, 1.0) ** 2
    X = np.column_stack([np.ones_like(l), l])
    W = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * W[:, None], chi2 * W, rcond=None)
    return dict(intercept=float(coef[0]), c=float(max(coef[1], 0.0)))


def run_xcal(weights, gwas, panel_genotypes, panel_snps, panel_pos, c, window_bp=1_000_000):
    """XCAL for every gene.

    weights: DataFrame[gene, snp, weight] (allele-dosage scale, same effect allele as the panel)
    gwas: DataFrame[snp, z] from the target-ancestry GWAS (aligned alleles)
    panel_genotypes: n x p dosage matrix of an LD panel of the SAME ancestry as the GWAS
    panel_snps, panel_pos: SNP ids and base-pair positions for the panel columns (one chromosome)
    c: polygenic parameter N*h2/M (use estimate_c with LD scores computed on the same SNP set)
    """
    panel_snps = np.asarray(panel_snps)
    panel_pos = np.asarray(panel_pos)
    col = {s_: i for i, s_ in enumerate(panel_snps)}
    zmap = dict(zip(gwas["snp"], gwas["z"]))
    out = []
    for gene, g in weights.groupby("gene"):
        g = g[g.snp.isin(col.keys()) & g.snp.isin(zmap.keys())]
        if g.empty:
            continue
        widx = np.array([col[s_] for s_ in g.snp])
        lo, hi = panel_pos[widx].min() - window_bp, panel_pos[widx].max() + window_bp
        win = np.where((panel_pos >= lo) & (panel_pos <= hi))[0]
        R, s, n = ld_from_genotypes(panel_genotypes[:, win])
        pos_in_win = {j: k for k, j in enumerate(win)}
        w = np.zeros(len(win))
        z = np.zeros(len(win))
        for j, wt, sn in zip(widx, g.weight.values, g.snp.values):
            w[pos_in_win[j]] = wt
            z[pos_in_win[j]] = zmap[sn]
        z_x, p, z_std, q = xcal_gene(w, z, R, s, n, c)
        out.append(dict(gene=gene, n_weight_snps=len(widx), n_window_snps=len(win), z_standard=z_std,
                        q_target=q, inflation_factor=1 + c * q, z_xcal=z_x, p_xcal=p))
    return pd.DataFrame(out)
