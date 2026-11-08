import argparse, glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

ap = argparse.ArgumentParser(description="Summarise the null benchmark: counted type-I error with region-bootstrap CIs")
ap.add_argument("outdir")
ap.add_argument("--benchdir", help="directory with region_*.csv.gz (+ gc_region_*.csv.gz) from 03_null_benchmark.py")
ap.add_argument("--combined", help="combined per-gene raw file (results/raw/null_benchmark_per_gene.csv.gz) instead of --benchdir")
ap.add_argument("--lambda-only", action="store_true", help="only compute genomic-control lambdas from the pass-1 grids")
ap.add_argument("--boot", type=int, default=2000)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11.5, "legend.fontsize": 10})
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
GRID = np.logspace(-3, 3.5, 2000)

if a.lambda_only:
    grids = {}
    for f in sorted(glob.glob(f"{a.benchdir}/grid_*.npz")):
        z = np.load(f)
        for k in z.files:
            grids[k] = grids.get(k, 0) + z[k]
    d = pd.concat([pd.read_csv(f, usecols=["gwas_pop", "pi", "c", "reps"]) for f in sorted(glob.glob(f"{a.benchdir}/region_*.csv.gz")) if os.path.getsize(f) > 30])
    ntot = d.groupby(["gwas_pop", "pi", "c"]).reps.sum()
    rows = []
    for (g, pi, c), n in ntot.items():
        surv = grids[f"{g}|{pi}|{c}"] / n
        med = float(np.interp(0.5, surv[::-1], GRID[::-1]))
        rows.append(dict(gwas_pop=g, pi=pi, c=c, n_tests=int(n), median_chi2=med, lambda_gc=max(med / chi2.isf(0.5, 1), 1.0)))
    pd.DataFrame(rows).to_csv(f"{a.outdir}/table_gc_lambda.csv", index=False)
    print("lambdas written:", len(rows))
    raise SystemExit

if a.combined:
    d = pd.read_csv(a.combined)
else:
    d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{a.benchdir}/region_*.csv.gz")) if os.path.getsize(f) > 30], ignore_index=True)
    gcf = sorted(glob.glob(f"{a.benchdir}/gc_region_*.csv.gz"))
    if gcf:
        gc = pd.concat([pd.read_csv(f) for f in gcf if os.path.getsize(f) > 30], ignore_index=True)
        d = d.merge(gc.drop(columns=["reps"]), on=["gene", "region", "gwas_pop", "pi", "c"], how="left", validate="one_to_one")
        os.makedirs(f"{a.outdir}/raw", exist_ok=True)
        d.to_csv(f"{a.outdir}/raw/null_benchmark_per_gene.csv.gz", index=False)

METHODS = {"ORACLE": "Oracle (exact null variance)", "STD_M": "Standard TWAS, matched LD", "STD_E": "Standard TWAS, EUR LD",
           "VC_E_E": "Adapted variance control (EUR q), EUR LD", "VC_E_M": "Adapted variance control (EUR q), matched LD",
           "XCAL": "XCAL (normal reference)", "XCALT": "XCAL-t (GWAS SDs; recommended)", "XCALT_PSD": "XCAL-t (panel SDs)",
           "GC": "Genomic control (median chi2)"}
PHI_LIANG = 4.2e-5
qbar_eur = d[(d.gwas_pop == "CEU") & (d.c == 0) & (d.pi == 1)].qE.mean()
M_EFF = qbar_eur / PHI_LIANG
regions = np.sort(d.region.unique())
rng = np.random.default_rng(1)
W = rng.multinomial(len(regions), np.full(len(regions), 1 / len(regions)), size=a.boot)
rows = []
for (g, pi, c), s in d.groupby(["gwas_pop", "pi", "c"]):
    tests_r = s.groupby("region").reps.sum().reindex(regions, fill_value=0).values
    bt = W @ tests_r
    base = dict(gwas_pop=g, pi=pi, c=c, Nh2_equiv=int(round(c * M_EFF)), n_gene_models=len(s), n_tests=int(tests_r.sum()))
    for meth, lab in METHODS.items():
        if f"n_{meth}_0.05" not in s:
            continue
        r = dict(base, method=lab)
        for al in ALPHAS:
            cnt_r = s.groupby("region")[f"n_{meth}_{al:g}"].sum().reindex(regions, fill_value=0).values
            bc = W @ cnt_r / bt
            r.update({f"count_{al:g}": int(cnt_r.sum()), f"expected_{al:g}": tests_r.sum() * al, f"rate_{al:g}": cnt_r.sum() / tests_r.sum(),
                      f"ratio_{al:g}": cnt_r.sum() / tests_r.sum() / al, f"lo_{al:g}": np.quantile(bc, .025) / al,
                      f"hi_{al:g}": np.quantile(bc, .975) / al})
            if f"pred_{meth}_{al:g}" in s:
                r[f"pred_ratio_{al:g}"] = float((s[f"pred_{meth}_{al:g}"] * s.reps).sum() / s.reps.sum()) / al
        rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{a.outdir}/table_null_counts.csv", index=False)

orc = T[T.method == METHODS["ORACLE"]].set_index(["gwas_pop", "pi", "c"])
RO = T[["gwas_pop", "pi", "c", "Nh2_equiv", "method"]].copy()
for al in ALPHAS:
    RO[f"ratio_to_oracle_{al:g}"] = T[f"rate_{al:g}"].values / orc.loc[list(zip(T.gwas_pop, T.pi, T.c)), f"rate_{al:g}"].values
RO.to_csv(f"{a.outdir}/table_ratio_to_oracle.csv", index=False)

dense = d[(d.pi == 1) & (d.c == 0)].copy()
ceu_q = dense[dense.gwas_pop == "CEU"].set_index("gene").q_true
dense["phi_ratio_true"] = dense.q_true / dense.gene.map(ceu_q)
dense["phi_ratio_panel"] = dense.qT / dense.qE
phi = dense.groupby("gwas_pop")[["phi_ratio_true", "phi_ratio_panel"]].describe(percentiles=[.05, .25, .5, .75, .95]).T.reset_index()
phi.columns = ["statistic_source", "summary"] + list(phi.columns[2:])
phi.to_csv(f"{a.outdir}/table_phi_transfer.csv", index=False)

eb = dense.groupby("gwas_pop").apply(lambda x: pd.Series(dict(
    n_models=len(x), relbias_Q=x.Q_G.mean() / x.Q_true.mean() - 1, relbias_Q_uncorrected=x.Q_G_raw.mean() / x.Q_true.mean() - 1,
    relbias_Q_v01_floor=np.maximum(x.Q_G_untrunc, x.L_G).mean() / x.Q_true.mean() - 1, relbias_L=x.L_G.mean() / x.L_true.mean() - 1,
    median_abs_rel_err_Q=np.median(np.abs(x.Q_G / x.Q_true - 1)), floor_rate=x.floor_G.mean(),
    floor_rate_v01=(x.Q_G_untrunc < x.L_G).mean())), include_groups=False).reset_index()
eb.to_csv(f"{a.outdir}/table_estimator_genotype_panels.csv", index=False)

pu = dense.groupby(["gwas_pop", "train_pop", "n_eqtl"]).agg(n_models=("gene", "size"), mean_r2_target=("r2_target", "mean"),
                                                             tests=("reps", "sum"), rejections_005=("n_STD_M_0.05", "sum")).reset_index()
pu["type1_005"] = pu.rejections_005 / pu.tests
pu.to_csv(f"{a.outdir}/table_prediction_uncertainty_c0.csv", index=False)

dd = d[d.pi == 1]
DF = dd.groupby(["gwas_pop", "c"]).agg(median_df=("df", "median"), p05_df=("df", lambda x: x.quantile(.05)),
                                      frac_df_at_lower_clip=("df", lambda x: (x <= 4.0001).mean()),
                                      frac_df_at_upper_clip=("df", lambda x: (x >= 99999).mean()),
                                      median_df_panel_sd=("df_psd", "median")).reset_index()
DF.to_csv(f"{a.outdir}/table_xcalt_df.csv", index=False)

POPN = {"CEU": "European-ancestry GWAS", "CHB": "East Asian-ancestry GWAS", "YRI": "African-ancestry GWAS"}
PLOT = [("Standard TWAS, matched LD", "#2a78d6", "-"), ("Standard TWAS, EUR LD", "#eb6834", "--"),
        ("Adapted variance control (EUR q), EUR LD", "#1baf7a", "-"), ("Adapted variance control (EUR q), matched LD", "#eda100", ":"),
        ("Genomic control (median chi2)", "#8a8a8a", "-."), ("XCAL (normal reference)", "#6250d6", "-"),
        ("XCAL-t (GWAS SDs; recommended)", "#e34948", "-")]
for al, fname in [(1e-4, "fig_counted_type1_alpha1e-4.png"), (2.5e-6, "fig_counted_type1_genomewide.png")]:
    fig, axs = plt.subplots(1, 3, figsize=(14, 5.4), sharey=True)
    for ax, g in zip(axs, ["CEU", "CHB", "YRI"]):
        for lab, col, sty in PLOT:
            s_ = T[(T.gwas_pop == g) & (T.method == lab) & (T.pi == 1)].sort_values("c")
            y = np.maximum(s_[f"ratio_{al:g}"].values, 0.05)
            ax.errorbar(s_.Nh2_equiv / 1000, y, yerr=[y - np.maximum(s_[f"lo_{al:g}"], 0.05), np.maximum(s_[f"hi_{al:g}"] - y, 0)],
                        fmt=sty, color=col, lw=1.8, marker="o", ms=4, capsize=2, label=lab)
        ax.axhline(1, color="#555", lw=1); ax.set_yscale("log"); ax.set_title(POPN[g])
        ax.set_xlabel("Polygenic signal N·h² (thousands)"); ax.grid(alpha=.25, lw=.5)
    axs[0].set_ylabel(f"Counted type-I error / nominal (α = {al:g})")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.15, 1, 1)); fig.savefig(f"{a.outdir}/{fname}", dpi=300); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
for ax, c in zip(axs, [0.1, 0.4]):
    for g, col in zip(["CEU", "CHB", "YRI"], ["#2a78d6", "#eda100", "#1baf7a"]):
        for meth, sty in [("XCAL-t (GWAS SDs; recommended)", "-"), ("Oracle (exact null variance)", ":")]:
            s_ = T[(T.gwas_pop == g) & (T.method == meth) & (T.c == c)].sort_values("pi")
            ax.plot(s_.pi, np.maximum(s_["ratio_2.5e-06"], 0.05), sty, color=col, marker="o", ms=4,
                    label=f"{g} {'XCAL-t' if meth.startswith('XCAL') else 'oracle'}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.axhline(1, color="#555", lw=1)
    ax.set_xlabel("Fraction of SNPs with non-zero effect (π)"); ax.set_title(f"N·h² ≈ {int(round(c * M_EFF / 1000))}k")
    ax.grid(alpha=.25, lw=.5)
axs[0].set_ylabel("Counted type-I error / nominal (α = 2.5e-6)")
axs[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig_sparse_architectures.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.4))
ax.boxplot([dense[dense.gwas_pop == g].phi_ratio_true.values for g in ["CEU", "CHB", "YRI"]], tick_labels=["EUR", "EAS", "AFR"],
           whis=(5, 95), showfliers=False)
ax.axhline(1, color="#888", lw=1); ax.set_ylabel("q(GWAS-ancestry LD) / q(EUR LD)"); ax.set_xlabel("GWAS ancestry")
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig_phi_transfer.png", dpi=300); plt.close(fig)

with open(f"{a.outdir}/summary_null.txt", "w", encoding="utf-8") as fh:
    fh.write(f"M_eff={M_EFF:.0f} (qbar_EUR={qbar_eur:.3f}, Phi_Liang={PHI_LIANG})\n")
    fh.write(f"gene models per ancestry: {dense.groupby('gwas_pop').size().to_dict()}\n")
    cols = ["gwas_pop", "pi", "c", "Nh2_equiv", "method", "n_tests", "ratio_0.05", "ratio_0.0001", "lo_0.0001", "hi_0.0001",
            "ratio_2.5e-06", "lo_2.5e-06", "hi_2.5e-06"]
    fh.write(T[cols].round(3).to_string() + "\n\n" + phi.round(3).to_string() + "\n\n" + eb.round(4).to_string() + "\n\n" + DF.round(2).to_string() + "\n")
print(open(f"{a.outdir}/summary_null.txt", encoding="utf-8").read()[:3000])
