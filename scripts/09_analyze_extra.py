import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from xcal.calibration import type1, type1_t

ap = argparse.ArgumentParser(description="Tables and figures for the individual-level check, panel size, power and estimated-c experiments")
ap.add_argument("--indiv"); ap.add_argument("--panel"); ap.add_argument("--power"); ap.add_argument("--estc")
ap.add_argument("outdir")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
rng = np.random.default_rng(3)
out = []

if a.indiv:
    d = pd.read_csv(a.indiv)
    d = d.drop_duplicates(["region", "gene", "gwas_pop", "c"])
    regions = d.region.unique()
    W = rng.multinomial(len(regions), np.full(len(regions), 1 / len(regions)), size=2000)
    rows = []
    for (g, c), s in d.groupby(["gwas_pop", "c"]):
        base = dict(gwas_pop=g, c=c, gene_models=len(s), phenotype_replicates=int(s.reps.iloc[0]),
                    gene_phenotype_evaluations=int(s.reps.sum()))
        for m in ["STD_M", "STD_E", "VC_E_E", "XCAL", "XCALT"]:
            r = dict(base, method=m)
            if m != "XCALT":
                ratio = s[f"emp_var_{m}"] / s[f"pred_var_{m}"]
                rr = ratio.groupby(s.region).mean().reindex(regions).values
                nn = s.groupby("region").size().reindex(regions).values
                bs = W @ (rr * nn) / (W @ nn)
                r.update(var_ratio_mean=float(ratio.mean()), var_ratio_lo=float(np.quantile(bs, .025)), var_ratio_hi=float(np.quantile(bs, .975)))
            for al in (0.05, 0.01):
                cnt = s.groupby("region")[f"rej_{m}_{al:g}"].sum().reindex(regions).values
                tests = s.groupby("region").reps.sum().reindex(regions).values
                bs = (W @ cnt) / (W @ tests)
                r.update({f"rate_{al:g}": cnt.sum() / tests.sum(), f"lo_{al:g}": np.quantile(bs, .025), f"hi_{al:g}": np.quantile(bs, .975),
                          f"pred_rate_{al:g}": float(s[f"pred_rej_{m}_{al:g}"].mean())})
            rows.append(r)
    T = pd.DataFrame(rows); T.to_csv(f"{a.outdir}/table_individual_level.csv", index=False)
    out.append(("individual-level", T.round(4)))

if a.panel:
    d = pd.read_csv(a.panel)
    rows = []
    for (g, n, c), s in d.groupby(["gwas_pop", "panel_n", "c"]):
        r = dict(gwas_pop=g, panel_n=n, c=c, n_models=len(s), floor_rate=s.floor_active.mean(), median_df=s.df_XCAL.median())
        for al in (1e-4, 2.5e-6):
            r[f"XCAL_{al:g}"] = type1(s.V_XCAL, al).mean() / al
            r[f"XCALt_{al:g}"] = type1_t(s.V_XCAL, s.df_XCAL, al).mean() / al
            r[f"VCeur_{al:g}"] = type1(s.V_VCeur, al).mean() / al
        rows.append(r)
    T = pd.DataFrame(rows); T.to_csv(f"{a.outdir}/table_panel_size.csv", index=False)
    out.append(("panel size (predicted ratios)", T.round(3)))

if a.power:
    d = pd.read_csv(a.power)
    meth = ["Oracle", "XCAL", "XCALt", "VCeur", "Standard"]
    T = d.groupby(["gwas_pop", "N", "q2"])[[f"power_{m}" for m in meth] + [f"type1_{m}" for m in meth]].mean().reset_index()
    for m in ["XCAL", "XCALt", "VCeur"]:
        T[f"rel_loss_{m}_vs_oracle"] = 1 - T[f"power_{m}"] / T.power_Oracle
    T["n_gene_models"] = d.groupby(["gwas_pop", "N", "q2"]).size().values
    T.to_csv(f"{a.outdir}/table_power.csv", index=False)
    out.append(("power", T.round(3)))
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, g in zip(axs, ["CEU", "CHB", "YRI"]):
        for m, col, st in [("Oracle", "#222", ":"), ("XCAL", "#6250d6", "-"), ("XCALt", "#e34948", "-"), ("VCeur", "#1baf7a", "--")]:
            for N, mk in [(50_000, "o"), (200_000, "s")]:
                s_ = T[(T.gwas_pop == g) & (T.N == N)].sort_values("q2")
                ax.plot(s_.q2 * 100, s_[f"power_{m}"], st, color=col, marker=mk, ms=4,
                        label=f"{ {'XCALt': 'XCAL-t', 'VCeur': 'VC (EUR Φ)'}.get(m, m)}, N={N // 1000}k")
        ax.set_xscale("log"); ax.set_xlabel("Trait variance explained by gene (%)"); ax.grid(alpha=.25, lw=.5)
        ax.set_title({"CEU": "European-ancestry GWAS", "CHB": "East Asian-ancestry GWAS", "YRI": "African-ancestry GWAS"}[g], fontsize=11)
    axs[0].set_ylabel("Power at α = 2.5e-6")
    h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.14, 1, 1)); fig.savefig(f"{a.outdir}/fig3_power.png", dpi=300); plt.close(fig)

if a.estc:
    T = pd.read_csv(a.estc); C = pd.read_csv(a.estc.replace(".csv", "_c_estimates.csv"))
    T.to_csv(f"{a.outdir}/table_estimated_c.csv", index=False); C.to_csv(f"{a.outdir}/table_estimated_c_values.csv", index=False)
    out.append(("estimated c: c-hat", C.round(4)))
    out.append(("estimated c: XCAL-t counted ratios", T.round(3)))

with open(f"{a.outdir}/summary_extra.txt", "w", encoding="utf-8") as fh:
    for name, T in out:
        fh.write(f"== {name}\n{T.to_string()}\n\n")
print(open(f"{a.outdir}/summary_extra.txt", encoding="utf-8").read()[:8000])
