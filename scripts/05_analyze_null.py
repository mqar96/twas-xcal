import sys, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from xcal.calibration import type1, type1_t

inp, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
d = pd.read_csv(inp)
d = d[d.n_w > 0]
PHI_LIANG = 4.2e-5
qbar_eur = d[(d.gwas_pop == "CEU") & (d.c == 0)].qE.mean()
M_eff = qbar_eur / PHI_LIANG
d["Nh2_equiv"] = d.c * M_eff
METHODS = {
    "V_unadj_matchedLD": "Standard TWAS, matched-ancestry LD",
    "V_unadj_EURLD": "Standard TWAS, European LD reference",
    "V_VCeur_EURLD": "Variance control (EUR Phi), EUR LD",
    "V_VCeur_matchedLD": "Variance control (EUR Phi), matched LD",
    "V_VCtarget_matchedLD": "Variance control (target Phi), matched LD",
    "V_XCAL": "XCAL (proposed)",
}
ALPHAS = [0.05, 1e-3, 1e-4, 2.5e-6]
rows = []
for (g, c), s in d.groupby(["gwas_pop", "c"]):
    for m, lab in METHODS.items():
        V = s[m].values
        r = dict(gwas_pop=g, c=c, Nh2_equiv=round(c * M_eff), method=lab, n_genes=len(V),
                 mean_V=V.mean(), median_V=np.median(V), frac_V_gt_1_1=(V > 1.1).mean(), frac_V_lt_0_9=(V < 0.9).mean())
        for a in ALPHAS:
            r[f"type1_{a:g}"] = type1(V, a).mean()
            r[f"ratio_{a:g}"] = type1(V, a).mean() / a
        rows.append(r)
    V = s["V_XCAL"].values; df_ = s["df_XCAL"].values
    r = dict(gwas_pop=g, c=c, Nh2_equiv=round(c * M_eff), method="XCAL-t (proposed, panel-uncertainty aware)", n_genes=len(V),
             mean_V=V.mean(), median_V=np.median(V), frac_V_gt_1_1=(V > 1.1).mean(), frac_V_lt_0_9=(V < 0.9).mean(), median_df=np.median(df_))
    for a in ALPHAS:
        r[f"type1_{a:g}"] = type1_t(V, df_, a).mean()
        r[f"ratio_{a:g}"] = type1_t(V, df_, a).mean() / a
    rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(f"{outdir}/table_null_calibration.csv", index=False)

gc = []
for (g, c), s in d.groupby(["gwas_pop", "c"]):
    lam = np.median(s.V_unadj_matchedLD)
    V = s.V_unadj_matchedLD.values / lam
    gc.append(dict(gwas_pop=g, c=c, lambda_gc=lam, **{f"type1_{a:g}": type1(V, a).mean() for a in ALPHAS}))
pd.DataFrame(gc).to_csv(f"{outdir}/table_genomic_control.csv", index=False)

d["phi_ratio_target_over_EUR"] = d.qT / d.qE
pt = d[d.c == 0].groupby("gwas_pop").phi_ratio_target_over_EUR.describe(percentiles=[.05, .25, .5, .75, .95])
pt.to_csv(f"{outdir}/table_phi_transfer.csv")

u = d[d.c == 0].groupby(["gwas_pop", "train_pop", "n_eqtl"]).agg(
    n=("V_unadj_matchedLD", "size"), mean_V=("V_unadj_matchedLD", "mean"),
    type1_005=("V_unadj_matchedLD", lambda v: type1(v.values, 0.05).mean()),
    mean_r2_target=("r2_target", "mean")).reset_index()
u.to_csv(f"{outdir}/table_prediction_uncertainty_c0.csv", index=False)

PLOT = [("Standard TWAS, matched-ancestry LD", "#2a78d6", "-"),
        ("Standard TWAS, European LD reference", "#eb6834", "--"),
        ("Variance control (EUR Phi), EUR LD", "#1baf7a", "-"),
        ("Variance control (EUR Phi), matched LD", "#eda100", ":"),
        ("XCAL (proposed)", "#6250d6", "-"),
        ("XCAL-t (proposed, panel-uncertainty aware)", "#e34948", "-")]
for alpha_col, fname, ylab in [("ratio_0.0001", "fig1_type1_inflation_alpha1e-4.png", "Type-I error / nominal (α = 1e-4)"),
                               ("ratio_2.5e-06", "fig1b_type1_inflation_genomewide.png", "Type-I error / nominal (α = 2.5e-6)")]:
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.8), sharey=True)
    for ax, g in zip(axs, ["CEU", "CHB", "YRI"]):
        for lab, col, sty in PLOT:
            s_ = T[(T.gwas_pop == g) & (T.method == lab)].sort_values("c")
            ax.plot(s_.Nh2_equiv / 1000, s_[alpha_col], sty, color=col, lw=2, marker="o", ms=4, label=lab)
        ax.axhline(1, color="#888", lw=1)
        ax.set_yscale("log")
        ax.set_title({"CEU": "European-ancestry GWAS", "CHB": "East Asian-ancestry GWAS", "YRI": "African-ancestry GWAS"}[g], fontsize=11)
        ax.set_xlabel("Polygenic signal N·h² (thousands)")
        ax.grid(alpha=.25, lw=.5)
    axs[0].set_ylabel(ylab)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8.5, frameon=False)
    fig.suptitle("TWAS null calibration across ancestries: 3,538 gene models, 500 simulated regions", fontsize=12)
    fig.tight_layout(rect=(0, 0.12, 1, 0.95)); fig.savefig(f"{outdir}/{fname}", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.5, 4))
data = [d[(d.c == 0) & (d.gwas_pop == g)].phi_ratio_target_over_EUR.values for g in ["CEU", "CHB", "YRI"]]
ax.boxplot(data, tick_labels=["EUR", "EAS", "AFR"], showfliers=False)
ax.axhline(1, color="#888", lw=1)
ax.set_ylabel("Φ(target LD) / Φ(EUR LD)"); ax.set_xlabel("GWAS ancestry")
ax.set_title("Polygenic inflation slope is ancestry-specific", fontsize=11)
fig.tight_layout(); fig.savefig(f"{outdir}/fig2_phi_transfer.png", dpi=200); plt.close(fig)
print("M_eff", round(M_eff), "qbar_eur", round(qbar_eur, 2))
print(T[T.c.isin([0, 0.1, 0.4])][["gwas_pop", "c", "Nh2_equiv", "method", "mean_V", "type1_0.05", "ratio_0.0001", "ratio_2.5e-06"]].round(3).to_string())
print(pt.round(3))
print(u.round(3).to_string())
