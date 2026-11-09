# Building the submission package

These scripts assemble the manuscript files and the supplementary workbook from the repository results. The manuscript source (`manuscript_source.md`, `cover_letter_source.md`) is distributed with the submission package, not in this repository.

```bash
python paper/build_supplement.py results results/application OUTDIR   # Supplementary_Tables_S1-S19.xlsx and figures/
cd OUTDIR && bash build_docs.sh                                        # .docx (Word equations) and .pdf (xelatex) with line and page numbers
```
Requirements: pandas, openpyxl, python-docx, pandoc ≥ 3, xelatex with the TeX Gyre Termes font and the `lineno` package.

Before running build_docs.sh in OUTDIR, copy that script and the manuscript_source.md and cover_letter_source.md files from the submission package into OUTDIR. The numerical figures are populated by build_supplement.py. The reviewed submission package also supplies a tested edit-replay script that preserves the existing Word mathematics.
