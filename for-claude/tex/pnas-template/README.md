PNAS two-column typeset preview, for a visual sense of the page count against PNAS's six-page limit (main text 7 pages
including references; the SI Appendix follows in one column from page 8; 22 pages in all).

Source of truth is tex/pnas/pnas-main.tex. Everything here is generated: run `python3 build_typeset.py` in this
directory to regenerate pnas-typeset.tex, copy the figure, table, SI and .dat files from tex/pnas, and compile
(pdflatex, bibtex, pdflatex x2). Do not edit pnas-typeset.tex by hand. The generator replaces SI \ref's by "SI Text N"
names and \citet by "Author et al.~\cite", copies figure files with \textwidth -> \linewidth (Fig. 1 at column width
with \footnotesize labels; Figs. 3 and 4 as figure*), removes the journal footer, DOI and date (page numbers only), and
appends the SI block after the references.

Class: pnas-new.cls v1.44 (2018) and pnasresearcharticle.sty from a public mirror of the Overleaf PNAS template
(dereckmezquita/latex-template-pnas); watermark disabled. Built 2026-09-29.
