# TWAS-XCAL: research brief

**MUHAMMAD QUDDUSE AHMAD RAO** · Independent researcher · [ORCID](https://orcid.org/0009-0001-8731-7080)

## Research question
How much error is introduced when a polygenic TWAS variance correction derived using European linkage disequilibrium is transferred to other GWAS populations, and how well can a finite ancestry-matched LD panel address it?

## What the project provides
- Python research software implementing XCAL and the approximate XCAL-t reference.
- A counted-null benchmark with 3,540 retained gene models in 500 simulated regions and 106,200,000 correlated gene–replicate evaluations per setting.
- A finite-panel variance estimator and checks of estimator bias, jackknife behavior, genomic control and ancestry transfer.
- Model-based 1000 Genomes/MESA comparisons and a descriptive Biobank Japan height application using 12,195 models.
- Public source, result tables, computational records, installed-package tests and CI.

## Key empirical findings
In the tested nonzero infinitesimal settings, the adapted European-slope correction reached approximately 61 times the nominal false-positive rate in East Asian GWAS simulations at alpha 2.5e-6; genomic control reached approximately 26 times nominal. Using a European slope with ancestry-matched LD over-corrected some African settings.

XCAL-t had counted ratios of approximately 0.98–1.22 times nominal at that threshold with known parameters. It remains mildly liberal in some settings. Under the tested sparse point-normal architectures, even a Gaussian-reference oracle using the exact unconditional null variance failed in the extreme tail.

## Where it may be useful
Research teams developing or evaluating gene–trait association methods can use the benchmark to study ancestry transfer, LD-reference-panel uncertainty and false-positive behavior. Potential audiences include statistical genetics, biobank analysis, computational target prioritization and scientific-software teams. These are potential research uses; no industrial adoption or clinical validation is claimed.

## What the evidence does not establish
The estimated-parameter extreme-tail calculations are approximate and reuse pilot standard errors while omitting intercept–slope covariance. The genotype-only benchmark predicts rates using estimated truth. The real-GWAS application has no null ground truth. The work does not establish universal calibration or identify causal genes. See [the detailed qualifications](SUBMISSION_REVIEW.md).

## Verify and discuss
- [Scientific software/results v0.3.0](https://doi.org/10.5281/zenodo.24062489)
- [Source and reproduction instructions](https://github.com/mqar96/twas-xcal)
- [Maintenance-code CI: seven successful jobs](https://github.com/mqar96/twas-xcal/actions/runs/42559890080), including installed-package tests, minimum dependencies, Windows and a small end-to-end workflow. Subsequent discovery changes affect documentation/metadata only.
- [Author and research links](AUTHOR.md)
- Research contact: qudduse1234@gmail.com

The associated manuscript is prepared for peer review; no journal acceptance is claimed. This brief describes the existing public software and results, and is not a new scientific validation study.
