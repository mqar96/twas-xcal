import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, t as tdist
from xcal.core import satterthwaite_df

ap = argparse.ArgumentParser(description="Summarise the real-genotype benchmark (model-based predicted type-I error; independent half-panel LD as truth)")
ap.add_argument("inp"); ap.add_argument("outdir")
ap.add_argument("--n-snps", type=int, default=1_836_406, help="M: number of SNPs in the panel (c = N*h2/M)")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
d = pd.read_csv(a.inp)
M = a.n_snps
NH2 = [0, 10_000, 30_000, 60_000, 120_000, 250_000]
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
POPN = {"AFR": "African", "EUR": "European", "EAS": "East Asian", "SAS": "South Asian", "AMR": "Admixed American"}
d = d[d.pred_perf_pval < 0.05].reset_index(drop=True)

u = d.drop_duplicates(["set", "gene"])
cov = u.groupby("set").agg(distinct_models=("gene", "size"), median_snps=("n_w", "median"), median_window=("n_win", "median"))
cov.loc["all"] = [len(u), u.n_w.median(), u.n_win.median()]
for p_ in sorted(d["pop"].unique()):
    cov[f"evaluable_{p_}"] = d[d["pop"] == p_].groupby("set").size().reindex(cov.index)
    cov.loc["all", f"evaluable_{p_}"] = int((d["pop"] == p_).sum())
cov.to_csv(f"{a.outdir}/real_coverage.csv")
d["phi_ratio"] = d.q_full / d.qE
pr = d.groupby("pop").phi_ratio.describe(percentiles=[.05, .25, .5, .75, .95])
pr.to_csv(f"{a.outdir}/real_phi_ratio.csv")
d.groupby(["set", "pop"]).phi_ratio.describe(percentiles=[.05, .5, .95]).to_csv(f"{a.outdir}/real_phi_ratio_by_set.csv")
fl = d.groupby("pop")[["floorA", "floorF"]].mean()
fl.to_csv(f"{a.outdir}/real_floor_activation.csv")

jl = np.array([np.array(x.split(";"), float) for x in d.jk_lin])
ju = np.array([np.array(x.split(";"), float) for x in d.jk_unb])


def t1(V, al, df=None):
    thr = norm.isf(al / 2) if df is None else tdist.isf(al / 2, df)
    return 2 * norm.sf(thr / np.sqrt(np.asarray(V, float)))


rows = []
for nh2 in NH2:
    c = nh2 / M
    numA = d.linB_vA + c * d.unbB_vA
    numE = d.linB_vE + c * d.unbB_vE
    V = {"Standard TWAS, matched LD": numA / d.linA, "Standard TWAS, EUR LD": numE / d.linE,
         "Variance control (EUR Φ), EUR LD": numE / d.linE / (1 + c * d.qE),
         "Variance control (EUR Φ), matched LD": numA / d.linA / (1 + c * d.qE),
         "XCAL (normal reference)": numA / (d.linA + c * d.unbA)}
    df_ = satterthwaite_df(jl + c * ju, (d.linA + c * d.unbA).values)
    for grp_cols in (["set", "pop"], ["pop"]):
        for key, idx in d.groupby(grp_cols).groups.items():
            key = key if isinstance(key, tuple) else (key,)
            base = dict(zip(grp_cols, key), Nh2=nh2, c=c, n_models=len(idx), pooling="by_set" if len(grp_cols) == 2 else "all_models")
            for m, v in V.items():
                r = dict(base, method=m, mean_V=float(v[idx].mean()))
                for al in ALPHAS:
                    r[f"ratio_{al:g}"] = float(t1(v[idx], al).mean() / al)
                rows.append(r)
            r = dict(base, method="XCAL-t (t reference)", mean_V=float(V["XCAL (normal reference)"][idx].mean()), median_df=float(np.median(df_[idx])))
            for al in ALPHAS:
                r[f"ratio_{al:g}"] = float(t1(V["XCAL (normal reference)"][idx], al, df_[idx]).mean() / al)
            rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{a.outdir}/real_null_calibration.csv", index=False)
pool = T[T.pooling == "all_models"].drop(columns=["set"])
pool.to_csv(f"{a.outdir}/real_null_calibration_all_models.csv", index=False)
macro = T[T.pooling == "by_set"].groupby(["pop", "Nh2", "method"])[[c for c in T if c.startswith("ratio")]].mean().reset_index()
macro.to_csv(f"{a.outdir}/real_null_calibration_macro_over_sets.csv", index=False)

order = ["AFR", "SAS", "EAS", "EUR", "AMR"]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.boxplot([d[d["pop"] == p].phi_ratio.values for p in order], tick_labels=[POPN[p] for p in order], showfliers=False)
ax.axhline(1, color="#888", lw=1); ax.set_ylabel("Φ(target-ancestry LD) / Φ(European LD)")
ax.set_title(f"Real genotypes: {len(u):,} MESA models (P < 0.05) x 2,490 1000 Genomes genomes", fontsize=10.5)
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig_real_phi_ratio.png", dpi=300); plt.close(fig)

PLOT = [("Standard TWAS, EUR LD", "#eb6834", "--"), ("Variance control (EUR Φ), EUR LD", "#1baf7a", "-"),
        ("Variance control (EUR Φ), matched LD", "#eda100", ":"), ("XCAL (normal reference)", "#6250d6", "-"), ("XCAL-t (t reference)", "#e34948", "-")]
fig, axs = plt.subplots(2, 5, figsize=(17, 7.6), sharey="row")
for row_i, al in enumerate([1e-4, 2.5e-6]):
    for ax, p in zip(axs[row_i], order):
        for m, col, st in PLOT:
            s_ = pool[(pool["pop"] == p) & (pool.method == m)].sort_values("Nh2")
            ax.plot(s_.Nh2 / 1000, s_[f"ratio_{al:g}"], st, color=col, lw=2, marker="o", ms=4, label=m)
        ax.axhline(1, color="#888", lw=1); ax.set_yscale("log"); ax.grid(alpha=.25, lw=.5)
        if row_i == 0: ax.set_title(f"{POPN[p]} GWAS", fontsize=10.5)
        if row_i == 1: ax.set_xlabel("N·h² (thousands)")
    axs[row_i][0].set_ylabel(f"Predicted type-I / nominal (α = {al:g})")
h, l = axs[0][0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, fontsize=8.5, frameon=False)
fig.suptitle("Real 1000 Genomes genotypes x MESA models: model-based predicted null calibration (independent half-panel LD as truth)", fontsize=11.5)
fig.tight_layout(rect=(0, 0.06, 1, 0.95)); fig.savefig(f"{a.outdir}/fig_real_null_calibration.png", dpi=300); plt.close(fig)

pd.set_option("display.width", 250)
with open(f"{a.outdir}/summary_real.txt", "w") as fh:
    fh.write(cov.to_string() + "\n\n" + pr.round(3).to_string() + "\n\n" + fl.round(4).to_string() + "\n\n")
    for al in ALPHAS:
        fh.write(f"alpha={al:g}\n" + pool[pool.Nh2.isin([0, 60000, 250000])].pivot_table(index=["pop", "Nh2"], columns="method", values=f"ratio_{al:g}").round(3).to_string() + "\n\n")
print(open(f"{a.outdir}/summary_real.txt").read())
