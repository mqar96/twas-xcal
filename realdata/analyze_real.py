import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import norm, t as tdist
inp, outdir = sys.argv[1], sys.argv[2]
d = pd.read_csv(inp)
M = 1836406
NH2 = [0, 10_000, 30_000, 60_000, 120_000, 250_000]
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
POPN = {"AFR": "African", "EUR": "European", "EAS": "East Asian", "SAS": "South Asian", "AMR": "Admixed American"}
d["sig"] = d.pred_perf_pval < 0.05
cov = d[d["pop"] == "EUR"].groupby("set").agg(genes=("gene", "nunique"), sig_genes=("sig", "sum"), median_snps=("n_w", "median"), median_window=("n_win", "median"))
cov.to_csv(f"{outdir}/real_coverage.csv")
d["phi_ratio"] = d.q_full / d.qE
pr = d[d.sig].groupby(["set", "pop"]).phi_ratio.describe(percentiles=[.05, .25, .5, .75, .95]).round(3)
pr.to_csv(f"{outdir}/real_phi_ratio.csv")
jl = np.array([np.array(x.split(";"), float) for x in d.jk_lin])
ju = np.array([np.array(x.split(";"), float) for x in d.jk_unb])
K = jl.shape[1]


def t1(V, a, df=None):
    V = np.asarray(V, float)
    thr = norm.isf(a / 2) if df is None else tdist.isf(a / 2, df)
    return 2 * norm.sf(thr / np.sqrt(V))

rows = []
for nh2 in NH2:
    c = nh2 / M
    numA = d.linB_vA + c * d.unbB_vA
    numE = d.linB_vE + c * d.unbB_vE
    V = {"Standard TWAS, matched LD": numA / d.linA,
         "Standard TWAS, EUR LD": numE / d.linE,
         "Variance control (EUR Phi), EUR LD": numE / d.linE / (1 + c * d.qE),
         "Variance control (EUR Phi), matched LD": numA / d.linA / (1 + c * d.qE),
         "XCAL": numA / (d.linA + c * d.unbA)}
    Dk = jl + c * ju; Df = (d.linA + c * d.unbA).values
    tau2 = (K - 1) / K * ((Dk - Dk.mean(1, keepdims=True)) ** 2).sum(1) / Df ** 2
    df_ = np.clip(2 / np.maximum(tau2, 1e-8), 4, 1e5)
    R_full = (1 + c * d.q_full) / (1 + c * d.qE)
    for (s, p), idx in d[d.sig].groupby(["set", "pop"]).groups.items():
        base = dict(set=s, pop=p, Nh2=nh2, c=c, n_genes=len(idx))
        for m, v in V.items():
            r = dict(base, method=m, mean_V=float(v[idx].mean()))
            for a in ALPHAS: r[f"ratio_{a:g}"] = float(t1(v[idx], a).mean() / a)
            rows.append(r)
        r = dict(base, method="XCAL-t", mean_V=float(V["XCAL"][idx].mean()))
        for a in ALPHAS: r[f"ratio_{a:g}"] = float(t1(V["XCAL"][idx], a, df_[d.index.get_indexer(idx)]).mean() / a)
        rows.append(r)
        r = dict(base, method="Residual of EUR-Phi variance control (full panels, ancestry bias only)", mean_V=float(R_full[idx].mean()))
        for a in ALPHAS: r[f"ratio_{a:g}"] = float(t1(R_full[idx], a).mean() / a)
        rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{outdir}/real_null_calibration.csv", index=False)
pool = T.groupby(["pop", "Nh2", "method"])[[c for c in T if c.startswith("ratio") or c == "mean_V"]].mean().reset_index()
pool.to_csv(f"{outdir}/real_null_calibration_pooled_sets.csv", index=False)

fig, ax = plt.subplots(figsize=(7.5, 4.2))
sub = d[d.sig]
order = ["AFR", "SAS", "EAS", "EUR", "AMR"]
data = [sub[sub["pop"] == p].phi_ratio.values for p in order]
ax.boxplot(data, tick_labels=[POPN[p] for p in order], showfliers=False)
ax.axhline(1, color="#888", lw=1)
ax.set_ylabel("Φ(target-ancestry LD) / Φ(European LD)")
n_models = int((sub["pop"] == "EUR").sum()); ax.set_title(f"Real data: {n_models:,} MESA gene models ({sub.gene.nunique():,} genes) x 2,490 1000 Genomes genomes", fontsize=10.5)
fig.tight_layout(); fig.savefig(f"{outdir}/fig_real_phi_ratio.png", dpi=200); plt.close(fig)

PLOT = [("Standard TWAS, EUR LD", "#eb6834", "--"), ("Variance control (EUR Phi), EUR LD", "#1baf7a", "-"),
        ("Variance control (EUR Phi), matched LD", "#eda100", ":"), ("XCAL", "#6250d6", "-"), ("XCAL-t", "#e34948", "-")]
fig, axs = plt.subplots(1, 5, figsize=(17, 4.4), sharey=True)
for ax, p in zip(axs, order):
    for m, col, st in PLOT:
        s_ = pool[(pool["pop"] == p) & (pool.method == m)].sort_values("Nh2")
        ax.plot(s_.Nh2 / 1000, s_["ratio_0.0001"], st, color=col, lw=2, marker="o", ms=4, label=m)
    ax.axhline(1, color="#888", lw=1); ax.set_yscale("log"); ax.set_title(f"{POPN[p]} GWAS", fontsize=10.5)
    ax.set_xlabel("N·h² (thousands)"); ax.grid(alpha=.25, lw=.5)
axs[0].set_ylabel("Type-I / nominal (α = 1e-4)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, fontsize=8.5, frameon=False)
fig.suptitle("Real 1000 Genomes genotypes + real MESA expression models: null calibration (split-half panels)", fontsize=11.5)
fig.tight_layout(rect=(0, 0.1, 1, 0.94)); fig.savefig(f"{outdir}/fig_real_null_calibration.png", dpi=200); plt.close(fig)
print(cov.to_string()); print(pr.to_string())
pd.set_option("display.width", 250)
print(pool[pool.Nh2.isin([0, 60000, 250000])].pivot_table(index=["pop", "Nh2"], columns="method", values="ratio_0.0001").round(2).to_string())
print(pool[pool.Nh2.isin([0, 60000, 250000])].pivot_table(index=["pop", "Nh2"], columns="method", values="ratio_2.5e-06").round(2).to_string())
