#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
M=Rao_TWAS_XCAL_GeneticEpidemiology_manuscript
C=Rao_cover_letter_GeneticEpidemiology
python3 - <<'PY'
s = open("manuscript_source.md").read()
legends = s.split("## Figure legends")[1]
figs = ["Figure1_counted_type1_genomewide.png", "Figure2_inflation_slope_ratio_simulation.png", "Figure3_sparse_architectures.png",
        "Figure4_power.png", "Figure5_real_inflation_slope_ratio.png", "Figure6_real_genotype_calibration.png", "Figure7_bbj_height_qq.png"]
leg = [l.strip() for l in legends.strip().split("\n\n") if l.strip().startswith("**Figure")]
assert len(leg) == len(figs), (len(leg), len(figs))
block = "\n\n## Figures\n\n" + "\n\n".join(f"![]({'figures/' + f}){{width=16cm}}\n\n{l}\n\n" for f, l in zip(figs, leg))
import re
def sup(t):
    parts = re.split(r"(\$\$.*?\$\$|\$[^$]*\$)", t, flags=re.S)
    return "".join(p if p.startswith("$") else re.sub(r"(\d)\^(-?[0-9.]+)", r"\1^\2^", p) for p in parts)
open("manuscript_with_figures.md", "w").write(sup(s + block))
PY
pandoc manuscript_with_figures.md -o $M.docx --resource-path=.:figures
pandoc cover_letter_source.md -o $C.docx
python3 - "$M.docx" "$C.docx" <<'PY'
import sys
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt


def page_field(par):
    for t, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        r = par.add_run()
        if t:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), t)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt
        r._r.append(e)


for i, f in enumerate(sys.argv[1:]):
    d = Document(f)
    for st in ("Normal", "Body Text", "First Paragraph", "Compact"):
        if st in [s.name for s in d.styles]:
            d.styles[st].font.name = "Times New Roman"; d.styles[st].font.size = Pt(12 if i == 0 else 11)
            if st != "Compact":
                d.styles[st].paragraph_format.line_spacing = 2.0 if i == 0 else 1.15
    for sec in d.sections:
        if i == 1:
            from docx.shared import Cm
            sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Cm(2.0)
        if i == 0:
            ln = OxmlElement("w:lnNumType")
            ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "360")
            sp = sec._sectPr
            nxt = [c for c in sp if c.tag in {qn("w:pgNumType"), qn("w:cols"), qn("w:formProt"), qn("w:vAlign"), qn("w:noEndnote"), qn("w:titlePg"), qn("w:textDirection"), qn("w:bidi"), qn("w:rtlGutter"), qn("w:docGrid")}]
            nxt[0].addprevious(ln) if nxt else sp.append(ln)
        p = sec.footer.paragraphs[0]; p.alignment = 1
        e = OxmlElement("w:suppressLineNumbers"); p._p.get_or_add_pPr().append(e)
        page_field(p)
    d.save(f)
PY
cat > /tmp/xcal_header.tex <<'TEX'
\usepackage{lineno}
\linenumbers
\usepackage{newunicodechar}
\newunicodechar{α}{\ensuremath{\alpha}}
\newunicodechar{λ}{\ensuremath{\lambda}}
\newunicodechar{π}{\ensuremath{\pi}}
\newunicodechar{χ}{\ensuremath{\chi}}
\newunicodechar{ν}{\ensuremath{\nu}}
\newunicodechar{Φ}{\ensuremath{\Phi}}
\newunicodechar{τ}{\ensuremath{\tau}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{≥}{\ensuremath{\geq}}
\newunicodechar{≤}{\ensuremath{\leq}}
\newunicodechar{×}{\ensuremath{\times}}
\newunicodechar{−}{\ensuremath{-}}
\newunicodechar{·}{\ensuremath{\cdot}}
\newunicodechar{²}{\ensuremath{^{2}}}
\newunicodechar{⁻}{\ensuremath{^{-}}}
\usepackage{float}
\floatplacement{figure}{H}
TEX
pandoc manuscript_with_figures.md -o $M.pdf --pdf-engine=xelatex --resource-path=.:figures -H /tmp/xcal_header.tex \
  -V mainfont="TeX Gyre Termes" -V fontsize=12pt -V geometry:margin=2.5cm -V linestretch=2 -V colorlinks=true
soffice --headless --convert-to pdf $C.docx >/dev/null 2>&1
pdfinfo $M.pdf | grep Pages
pdfinfo $C.pdf | grep Pages
