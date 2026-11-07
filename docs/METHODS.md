# XCAL: theory

## Setup
- Weights `w` (allele-dosage scale) are trained in ancestry *A* and applied to GWAS summary z-scores from ancestry *G*.
- The standard summary-statistic TWAS (S-PrediXcan form), using LD reference panel *r* with SNP SDs `s_r` and correlation `R_r`:

  Z = v' z / sqrt(v' R_r v),  with v = w ∘ s_r

## Null distribution under a polygenic trait
Take the infinitesimal polygenic model in GWAS ancestry *G*: per-standardised-SNP effects β ~ N(0, h²/M), sample size N, LD Σ_G. Then

  z ~ N(0, Σ_G + c Σ_G²),  with c = N h² / M

So under H0 (the gene's genetically regulated expression has no effect beyond background),

  Var(Z) = [ v'Σ_G v + c · v'Σ_G² v ] / v' R_r v

This variance factors into:
- **LD-panel term:** v'Σ_G v / v'R_r v, which is 1 only if the reference matches the GWAS ancestry.
- **Polygenic term:** 1 + c · q_G, with q_G = v'Σ_G² v / v'Σ_G v. q_G is a weight-specific LD score in the GWAS ancestry. Liang et al.'s Φ equals q_G / M.

q_G is a property of the GWAS-ancestry LD, not of the training ancestry. A Φ computed from European LD therefore mis-scales the correction for African-ancestry GWAS (shorter LD gives smaller q) and East Asian-ancestry GWAS (longer LD gives larger q).

Prediction error in `w` changes power, not the null: under independent eQTL and GWAS samples, `w` is fixed with respect to z.

## XCAL statistic
  Z_XCAL = v_t' z / sqrt( v_t' R_t v_t + ĉ · Q̂_t(v_t) ),  with v_t = w ∘ s_t

- *t* is an LD panel **of the GWAS ancestry**, for example 1000 Genomes AFR/EAS/EUR.
- Q̂_t is a finite-sample-unbiased estimate of v'Σ²v from an n-sample panel with p SNPs (Wishart moments):

  Q̂_t = [ ||R_t v||² − (p/n) v'R_t v ] / (1 + 1/n)

  Without this correction, ||R_t v||² overstates v'Σ²v by about (p/n)·v'Σv. That is large for 1000G-sized panels (n ≈ 500, p ≈ 400 per Mb).
- ĉ = N h²/M is estimated genome-wide from the GWAS itself, as the LD-score-regression slope with ancestry-matched LD scores.

## Assumptions and limits (to be tested)
- Infinitesimal background. Sparse or heavy-tailed architectures are a phase-2 test.
- Independent eQTL and GWAS samples. Overlap combined with environmental correlation between expression and trait is a phase-2 test.
- Additive linear model, continuous trait. Case-control traits are phase 2.
- A homogeneous GWAS ancestry. Admixed GWAS are phase 2, where local-ancestry-specific LD is needed (see CADET).
