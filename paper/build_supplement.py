import json
import os
import shutil
import sys

import pandas as pd

res, app, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(f"{out}/figures", exist_ok=True)
S = [
    ("S1", "Estimator bias, Gaussian settings", f"{res}/table_estimator_bias_gaussian.csv",
     "Relative bias of the moment-corrected estimator with the Cauchy-Schwarz floor (Q), the v0.1 floor at the linear form, and the uncorrected plug-in, over 48 Gaussian settings (1,000 panels each). floor_rate = share of panels where the floor binds."),
    ("S2", "Estimator on genotype panels", f"{res}/table_estimator_genotype_panels.csv",
     "500-person simulated genotype panels versus 20,000-person LD, pooled over 3,540 gene models per ancestry."),
    ("S3", "Individual-level GWAS check", f"{res}/table_individual_level.csv",
     "N = 10,000 per ancestry; 40 regions; 160 gene models x 1,000 phenotypes = 160,000 correlated gene-phenotype evaluations per cell. Counted rates with 95% region-bootstrap CIs (lo, hi); pred = analytic prediction. STD_M/STD_E = standard TWAS with matched/European LD; VC_E_E = adapted variance control (EUR q, EUR LD); XCALT = XCAL-t with GWAS SDs."),
    ("S4", "Counted null calibration", f"{res}/table_null_counts.csv",
     "Counted rejections over 1.06e8 tests per ancestry x architecture x c (3,540 gene models x 30,000 replicates). ratio = rate/alpha; lo/hi = 95% region-bootstrap CI (2,000 resamples); pred_ratio = analytic Gaussian prediction. pi = fraction of SNPs with non-zero effect (1 = infinitesimal). Genomic control counted exactly in a second pass."),
    ("S5", "Ratio to oracle", f"{res}/table_ratio_to_oracle.csv", "Counted error rate of each method divided by that of the oracle (exact null variance) in the same setting, at each alpha."),
    ("S6", "LD-panel size", f"{res}/table_panel_size.csv", "Predicted type-I ratios (exact truth from 20,000-person LD) by LD-panel size; 40 regions, 159 gene models. VCeur = adapted variance control (EUR q, EUR LD)."),
    ("S7", "Estimated a and c: XCAL-t", f"{res}/table_estimated_c.csv",
     "End-to-end benchmark with a and c re-estimated leave-one-chromosome-out (10 pseudo-chromosomes) by LD-score regression in each of 300 simulated genomes per setting (1.062e6 tests). Counted ratios with two-way (genome x region) bootstrap CIs; conditional_ratio = mean conditional rate given each genome's fitted parameters / alpha (infinitesimal settings only)."),
    ("S8", "Estimated a and c: estimates", f"{res}/table_estimated_c_values.csv", "Mean and SD of c-hat and a-hat over 300 genomes, free or fixed intercept, with the mean block-jackknife SE."),
    ("S9", "Power", f"{res}/table_power.csv", "Analytic power at alpha = 2.5e-6 (60 regions, 240 gene models per ancestry); q2 = trait variance explained by the gene; rel_loss = 1 - power/oracle power."),
    ("S10", "Inflation slope ratio (simulation)", f"{res}/table_phi_transfer.csv", "q(GWAS ancestry)/q(EUR); true = 20,000-person LD; panel = 500-person panels with the corrected estimator."),
    ("S11", "XCAL-t degrees of freedom", f"{res}/table_xcalt_df.csv", "Distribution of the Satterthwaite df (GWAS-SD scaling) and the share at the lower (4) and upper (1e5) clips; median_df_panel_sd = panel-SD fallback."),
    ("S12", "Prediction uncertainty (c = 0)", f"{res}/table_prediction_uncertainty_c0.csv", "Counted type-I error at 0.05 of standard TWAS with matched LD and no polygenic background, by training ancestry and eQTL sample size."),
    ("S13", "Genomic-control lambdas", f"{res}/table_gc_lambda.csv", "lambda = median(chi2)/0.455 over all genes and replicates of each setting (pass 1), used for the exact GC counts in S4."),
    ("S14", "Real genotypes: slope ratio", f"{res}/real/real_phi_ratio.csv", "1000 Genomes x MESA (P < 0.05): q(target)/q(EUR) from full panels."),
    ("S15", "Real genotypes: calibration", f"{res}/real/real_null_calibration_all_models.csv",
     "Model-based predicted type-I ratios pooled over all MESA models; truth = independent half-panel LD (estimated). No real phenotypes."),
    ("S16", "Real genotypes: by model set", f"{res}/real/real_null_calibration.csv", "As S15, by MESA model set (AFA, CAU, HIS) and pooled."),
    ("S17", "Real genotypes: coverage", f"{res}/real/real_coverage.csv", "Models with nominal prediction P < 0.05 and SNP coverage."),
    ("S18", "Biobank Japan height: summary", f"{app}/application_summary.csv", "Median-based lambda and Bonferroni-significant models (alpha = 0.05/12,195) per method."),
    ("S19", "Biobank Japan height: per model", f"{app}/application_per_model.csv.gz",
     "Per-model z-scores and p-values for every method; df_xcalt = XCAL-t degrees of freedom; c_used/a_used = leave-one-chromosome-out LDSC estimates; q_target/q_eur = slopes in East Asian/European LD."),
]
with pd.ExcelWriter(f"{out}/Supplementary_Tables_S1-S19.xlsx") as xw:
    readme = pd.DataFrame([(a, b, d) for a, b, c, d in S], columns=["table", "title", "description"])
    readme.to_excel(xw, sheet_name="README", index=False)
    for tag, title, f, desc in S:
        pd.read_csv(f, index_col=None).to_excel(xw, sheet_name=tag, index=False)
    j = json.load(open(f"{app}/application_ldsc.json"))
    rows = [dict(fit="genome-wide, free intercept", **{k: j["fit_free"][k] for k in ("c", "se_c", "intercept", "se_intercept")}),
            dict(fit="genome-wide, intercept fixed at 1", **{k: j["fit_fixed"][k] for k in ("c", "se_c", "intercept", "se_intercept")})]
    rows += [dict(fit=f"leave out chr{k}", **v) for k, v in j.get("loco", {}).items()]
    pd.DataFrame(rows).to_excel(xw, sheet_name="S19b_LDSC", index=False)
FIG = [(f"{res}/fig_counted_type1_genomewide.png", "Figure1_counted_type1_genomewide.png"),
       (f"{res}/fig_phi_transfer.png", "Figure2_inflation_slope_ratio_simulation.png"),
       (f"{res}/fig_sparse_architectures.png", "Figure3_sparse_architectures.png"),
       (f"{res}/fig3_power.png", "Figure4_power.png"),
       (f"{res}/real/fig_real_phi_ratio.png", "Figure5_real_inflation_slope_ratio.png"),
       (f"{res}/real/fig_real_null_calibration.png", "Figure6_real_genotype_calibration.png"),
       (f"{app}/fig_application_qq.png", "Figure7_bbj_height_qq.png"),
       (f"{res}/fig_counted_type1_alpha1e-4.png", "FigureS1_counted_type1_alpha1e-4.png")]
for a, b in FIG:
    shutil.copy(a, f"{out}/figures/{b}")
print("ok")
