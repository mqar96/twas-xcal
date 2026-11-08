import argparse, glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2, norm

ap = argparse.ArgumentParser(description="Summarise the null benchmark: counted type-I error with region-bootstrap CIs")
ap.add_argument("benchdir"); ap.add_argument("outdir")
ap.add_argument("--boot", type=int, default=2000)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{a.benchdir}/region_*.csv.gz")) if os.path.getsize(f) > 30], ignore_index=True)
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
METHODS = {"ORACLE": "Oracle (true null variance)", "STD_M": "Standard TWAS, matched LD", "STD_E": "Standard TWAS, EUR LD",
           "VC_E_E": "Variance control (EUR Φ), EUR LD", "VC_E_M": "Variance control (EUR Φ), matched LD",
           "XCAL": "XCAL (normal reference)", "XCALT": "XCAL-t (t reference)", "GC": "Genomic control (median χ²)"}
PHI_LIANG = 4.2e-5
qbar_eur = d[(d.gwas_pop == "CEU") & (d.c == 0) & (d.pi == 1)].qE.mean()
M_EFF = qbar_eur / PHI_LIANG
GRID = np.logspace(-3, 3.5, 2000)
regions = np.sort(d.region.unique())
rng = np.random.default_rng(1)
W = rng.multinomial(len(regions), np.full(len(regions), 1 / len(regions)), size=a.boot)

grids = {}
for f in sorted(glob.glob(f"{a.benchdir}/grid_*.npz")):
    rid = int(os.path.basename(f)[5:9]); z = np.load(f)
    for k in z.files:
        grids.setdefault(k, {})[rid] = z[k]


def gc_rate(exceed, ntot, alpha):
    surv = exceed / ntot
    med = np.interp(0.5, surv[::-1], GRID[::-1])
    lam = max(med / chi2.isf(0.5, 1), 1.0)
    thr = lam * chi2.isf(alpha, 1)
    return float(np.exp(np.interp(np.log(thr), np.log(GRID), np.log(np.maximum(surv, 1e-300))))), lam


rows = []
for (g, pi, c), s in d.groupby(["gwas_pop", "pi", "c"]):
    tests_r = s.groupby("region").reps.sum().reindex(regions, fill_value=0).values
    base = dict(gwas_pop=g, pi=pi, c=c, Nh2_equiv=int(round(c * M_EFF)), n_gene_models=len(s), n_tests=int(tests_r.sum()))
    bt = W @ tests_r
    for meth, lab in METHODS.items():
        r = dict(base, method=lab)
        for al in ALPHAS:
            if meth == "GC":
                key = f"{g}|{pi}|{c}"
                ex_r = np.array([grids[key].get(rid, np.zeros(len(GRID))) for rid in regions])
                est, lam = gc_rate(ex_r.sum(0), tests_r.sum(), al)
                boots = [gc_rate(w @ ex_r, w @ tests_r, al)[0] for w in W[:300]]
                r.update({f"rate_{al:g}": est, f"ratio_{al:g}": est / al, f"lo_{al:g}": np.quantile(boots, .025) / al,
                          f"hi_{al:g}": np.quantile(boots, .975) / al, "lambda_gc": lam})
                continue
            cnt_r = s.groupby("region")[f"n_{meth}_{al:g}"].sum().reindex(regions, fill_value=0).values
            pred = float((s[f"pred_{meth}_{al:g}"] * s.reps).sum() / s.reps.sum())
            bc = W @ cnt_r / bt
            r.update({f"count_{al:g}": int(cnt_r.sum()), f"rate_{al:g}": cnt_r.sum() / tests_r.sum(),
                      f"ratio_{al:g}": cnt_r.sum() / tests_r.sum() / al, f"lo_{al:g}": np.quantile(bc, .025) / al,
                      f"hi_{al:g}": np.quantile(bc, .975) / al, f"pred_ratio_{al:g}": pred / al})
        rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{a.outdir}/table_null_counts.csv", index=False)
orc = T[T.method == METHODS["ORACLE"]].set_index(["gwas_pop", "pi", "c"])
RO = T.copy()
for al in ALPHAS:
    RO[f"ratio_to_oracle_{al:g}"] = RO[f"rate_{al:g}"].values / orc.loc[list(zip(RO.gwas_pop, RO.pi, RO.c)), f"rate_{al:g}"].values
RO[["gwas_pop", "pi", "c", "Nh2_equiv", "method"] + [f"ratio_to_oracle_{al:g}" for al in ALPHAS]].to_csv(f"{a.outdir}/table_ratio_to_oracle.csv", index=False)

dense = d[(d.pi == 1) & (d.c == 0)].copy()
ceu_q = dense[dense.gwas_pop == "CEU"].set_index("gene").q_true
dense["phi_ratio_panel"] = dense.qT / dense.qE
dense["phi_ratio_true"] = dense.q_true / dense.gene.map(ceu_q)
phi = dense.groupby("gwas_pop")[["phi_ratio_true", "phi_ratio_panel"]].describe(percentiles=[.05, .25, .5, .75, .95]).T
phi.to_csv(f"{a.outdir}/table_phi_transfer.csv")

est = dense.assign(Q_old=np.maximum(dense.Q_T_untrunc, dense.L_T))
eb = est.groupby("gwas_pop").apply(lambda x: pd.Series(dict(
    n_models=len(x), relbias_Q_v02=x.Q_T.mean() / x.Q_true.mean() - 1, relbias_Q_uncorrected=x.Q_T_raw.mean() / x.Q_true.mean() - 1,
    relbias_Q_v01_floor=x.Q_old.mean() / x.Q_true.mean() - 1, relbias_L=x.L_T.mean() / x.L_true.mean() - 1,
    median_abs_rel_err_Q=np.median(np.abs(x.Q_T / x.Q_true - 1)), floor_rate_v02=x.floor_T.mean(),
    floor_rate_v01=(x.Q_T_untrunc < x.L_T).mean())), include_groups=False)
eb.to_csv(f"{a.outdir}/table_estimator_genotype_panels.csv")

pu = dense.groupby(["gwas_pop", "train_pop", "n_eqtl"]).agg(
    n_models=("gene", "size"), mean_r2_target=("r2_target", "mean"),
    tests=("reps", "sum"), rejections_005=("n_STD_M_0.05", "sum")).reset_index()
pu["type1_005"] = pu.rejections_005 / pu.tests
pu.to_csv(f"{a.outdir}/table_prediction_uncertainty_c0.csv", index=False)

DF = d[d.pi == 1].groupby(["gwas_pop", "c"]).df.describe(percentiles=[.05, .5, .95]).reset_index()
DF["frac_df_at_floor4"] = d[d.pi == 1].groupby(["gwas_pop", "c"]).df.apply(lambda x: (x <= 4.0001).mean()).values
DF.to_csv(f"{a.outdir}/table_xcalt_df.csv", index=False)

POPN = {"CEU": "European-ancestry GWAS", "CHB": "East Asian-ancestry GWAS", "YRI": "African-ancestry GWAS"}
PLOT = [("Standard TWAS, matched LD", "#2a78d6", "-"), ("Standard TWAS, EUR LD", "#eb6834", "--"),
        ("Variance control (EUR Φ), EUR LD", "#1baf7a", "-"), ("Variance control (EUR Φ), matched LD", "#eda100", ":"),
        ("Genomic control (median χ²)", "#8a8a8a", "-."),
        ("XCAL (normal reference)", "#6250d6", "-"), ("XCAL-t (t reference)", "#e34948", "-")]
for al, fname in [(1e-4, "fig1_counted_type1_alpha1e-4.png"), (2.5e-6, "fig1b_counted_type1_genomewide.png")]:
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.8), sharey=True)
    for ax, g in zip(axs, ["CEU", "CHB", "YRI"]):
        for lab, col, sty in PLOT:
            s_ = T[(T.gwas_pop == g) & (T.method == lab) & (T.pi == 1)].sort_values("c")
            y = np.maximum(s_[f"ratio_{al:g}"].values, 0.05)
            ax.errorbar(s_.Nh2_equiv / 1000, y, yerr=[y - np.maximum(s_[f"lo_{al:g}"], 0.05), np.maximum(s_[f"hi_{al:g}"] - y, 0)],
                        fmt=sty, color=col, lw=1.8, marker="o", ms=3.5, capsize=2, label=lab)
        ax.axhline(1, color="#555", lw=1); ax.set_yscale("log"); ax.set_title(POPN[g], fontsize=11)
        ax.set_xlabel("Polygenic signal N·h² (thousands)"); ax.grid(alpha=.25, lw=.5)
    axs[0].set_ylabel(f"Counted type-I error / nominal (α = {al:g})")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8.5, frameon=False)
    fig.suptitle(f"Counted false positives under the infinitesimal null ({T.n_tests.max() / 1e6:.0f}M tests per point; 95% region-bootstrap CI)\n"
                 "In European GWAS, variance control with EUR Φ coincides exactly with XCAL", fontsize=11)
    fig.tight_layout(rect=(0, 0.13, 1, 0.93)); fig.savefig(f"{a.outdir}/{fname}", dpi=300); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for ax, c in zip(axs, [0.1, 0.4]):
    for g, col in zip(["CEU", "CHB", "YRI"], ["#2a78d6", "#eda100", "#1baf7a"]):
        for meth, sty in [("XCAL-t (t reference)", "-"), ("Oracle (true null variance)", ":")]:
            s_ = T[(T.gwas_pop == g) & (T.method == meth) & (T.c == c)].sort_values("pi")
            ax.plot(s_.pi, np.maximum(s_["ratio_2.5e-06"], 0.05), sty, color=col, marker="o", ms=4,
                    label=f"{g} {'XCAL-t' if meth.startswith('XCAL') else 'oracle'}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.axhline(1, color="#555", lw=1)
    ax.set_xlabel("Fraction of SNPs with non-zero effect (π)"); ax.set_title(f"N·h² ≈ {int(c * M_EFF / 1000)}k", fontsize=11)
    ax.grid(alpha=.25, lw=.5)
axs[0].set_ylabel("Counted type-I error / nominal (α = 2.5e-6)")
axs[1].legend(fontsize=8, frameon=False)
fig.suptitle("Sparse architectures: variance-based corrections cannot remove heavy tails", fontsize=11.5)
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig_sparse_architectures.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.8, 4))
data = [dense[dense.gwas_pop == g].phi_ratio_true.values for g in ["CEU", "CHB", "YRI"]]
ax.boxplot(data, tick_labels=["EUR", "EAS", "AFR"], showfliers=False)
ax.axhline(1, color="#888", lw=1); ax.set_ylabel("Φ(GWAS-ancestry LD) / Φ(EUR LD)"); ax.set_xlabel("GWAS ancestry")
ax.set_title("Polygenic inflation slope is ancestry-specific (20,000-person LD)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{a.outdir}/fig2_phi_transfer.png", dpi=300); plt.close(fig)

with open(f"{a.outdir}/summary_null.txt", "w") as fh:
    fh.write(f"M_eff={M_EFF:.0f} (qbar_EUR={qbar_eur:.3f}, Phi_Liang={PHI_LIANG})\n")
    fh.write(f"gene models per ancestry: {d[(d.pi == 1) & (d.c == 0)].groupby('gwas_pop').size().to_dict()}\n")
    cols = ["gwas_pop", "pi", "c", "Nh2_equiv", "method", "n_tests", "ratio_0.05", "ratio_0.0001", "lo_0.0001", "hi_0.0001",
            "ratio_2.5e-06", "lo_2.5e-06", "hi_2.5e-06", "pred_ratio_2.5e-06"]
    fh.write(T[cols].round(3).to_string() + "\n\n")
    fh.write(phi.round(3).to_string() + "\n\n" + eb.round(4).to_string() + "\n\n" + DF.round(2).to_string() + "\n")
print(open(f"{a.outdir}/summary_null.txt").read()[:6000])
