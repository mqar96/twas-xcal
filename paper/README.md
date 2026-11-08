# Building the submission package

These scripts assemble the manuscript files and the supplementary workbook from the repository results. The manuscript source (`manuscript.md`, `cover_letter.md`) is distributed with the submission package, not in this repository.

```bash
python paper/build_supplement.py results results/application OUTDIR   # Supplementary_Tables_S1-S19.xlsx and figures/
cd OUTDIR && bash build_docs.sh                                        # .docx (Word equations) and .pdf (xelatex) with line and page numbers
```
Requirements: pandas, openpyxl, python-docx, pandoc ≥ 3, xelatex with the TeX Gyre Termes font and the `lineno` package.
