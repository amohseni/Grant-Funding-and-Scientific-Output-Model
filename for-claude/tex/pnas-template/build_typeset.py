"""Generate pnas-typeset.tex (PNAS two-column preview) from ../pnas/pnas-main.tex and compile it.
Run from tex/pnas-template. Source of truth is tex/pnas/pnas-main.tex; do not edit pnas-typeset.tex by hand."""
import re, os, shutil, subprocess, glob
SRC = '../pnas'
s = open(f'{SRC}/pnas-main.tex').read()
abstract = s[s.index('\\begin{abstract}') + len('\\begin{abstract}'):s.index('\\end{abstract}')].strip()
sig = s[s.index('\\noindent\\textbf{Significance.}') + len('\\noindent\\textbf{Significance.}'):s.index('The United States spent')].strip()
body = s[s.index('The United States spent'):s.index('\\section*{Materials and Methods}')]
methods = s[s.index('\\section*{Materials and Methods}') + len('\\section*{Materials and Methods}'):s.index('\\clearpage')].strip()
si = s[s.index('\\section*{Supporting Information Appendix}'):s.index('\\bibliographystyle{plainnat}')]
names = {'app:gap-rule': 'SI Text 1', 'app:lottery': 'SI Text 2', 'si:extended': 'SI Text 3', 'app:technology': 'SI Text 4',
         'si:sensitivity': 'SI Text 5', 'si:timing': 'SI Text 6', 'app:simulation': 'SI Materials and Methods', 'sec:model': 'SI Materials and Methods'}
def fix(t):
    for k, v in names.items(): t = t.replace('\\ref{%s}' % k, v)
    t = re.sub(r'Table~\\ref\{tab:(dynopt|strategies|settings|particles)\}', lambda m: {'dynopt': 'Table S4', 'strategies': 'Table S1', 'settings': 'Table S2', 'particles': 'Table S3'}[m.group(1)], t)
    t = re.sub(r'Fig\.~\\ref\{fig:(heuristics|targeting-value|overtrust|concentration)\}', lambda m: {'heuristics': 'Fig. S1', 'targeting-value': 'Fig. S2', 'overtrust': 'Fig. S3', 'concentration': 'Fig. S6'}[m.group(1)], t)
    t = t.replace('Corollary~\\ref{cor:targeting-vanishes}', 'Corollary 2').replace('Proposition~\\ref{prop:lottery}', 'Proposition 2')
    t = t.replace('\\citet{Carnehl2024}', 'Carnehl et al.~\\cite{Carnehl2024}').replace('\\citet{Sorzano2026}', 'Sorzano and Pueche-Granados~\\cite{Sorzano2026}')
    t = t.replace('\\citet{LindnerNakamura2015}', 'Lindner and Nakamura~\\cite{LindnerNakamura2015}').replace('\\citet{Sakamoto2023}', 'Sakamoto~\\cite{Sakamoto2023}')
    return t
body = fix(body); methods = fix(methods); abstract = fix(abstract); sig = fix(sig)
body = body.replace('The United States spent close to a trillion dollars', '\\dropcap{T}he United States spent close to a trillion dollars', 1)
methods = methods.replace('\\textbf{Model.}', '\\subsection*{Model}').replace('\\textbf{Information and strategies.}', '\\subsection*{Information and strategies}').replace('\\textbf{Parameters and scoring.}', '\\subsection*{Parameters and scoring}')
methods = re.sub(r'%.*', '', methods)
si = si.replace('\\section*{Supporting Information Appendix}', '\\section*{Supporting Information Appendix}\n\\noindent A.~Mohseni, S.~DeDeo, and K.\\,J.\\,S.~Zollman, A model of optimal science funding: targeting the capability-resource gap.\n')
doc = r'''\documentclass[9pt,twocolumn,twoside]{pnas-new}
\templatetype{pnasresearcharticle}
\usepackage{pgfplots}
\pgfplotsset{compat=1.17}
\usetikzlibrary{arrows.meta}
\usepackage{makecell}
\usepackage{amsthm,mathtools}
\usepackage{cleveref}
\usepackage{booktabs,subcaption}
\newtheorem{theorem}{Theorem}
\newtheorem{prop}{Proposition}
\newtheorem{corollary}{Corollary}
\newtheorem{lemma}{Lemma}
\newtheorem{example}{Example}
\newtheorem{remark}{Remark}
\newtheorem*{proprestate}{Proposition~1}
\newtheorem*{definition}{Definition}
\DeclareMathOperator*{\argmax}{arg \ max}
\makeatletter
\fancypagestyle{firststyle}{\fancyhf{}\fancyfoot[R]{\footerfont\textbf{\thepage}}\renewcommand{\headrulewidth}{0pt}}
\fancyhf{}\fancyfoot[RO]{\footerfont\textbf{\thepage}}\fancyfoot[LE]{\footerfont\textbf{\thepage}}
\makeatother
\title{A model of optimal science funding: targeting the capability-resource gap}
\author[a,1]{Aydin Mohseni}
\author[b]{Simon DeDeo}
\author[a]{Kevin J.\,S. Zollman}
\affil[a]{Department of Philosophy, Carnegie Mellon University, Pittsburgh, PA 15213}
\affil[b]{Department of Social and Decision Sciences, Carnegie Mellon University, Pittsburgh, PA 15213}
\leadauthor{Mohseni}
\significancestatement{''' + sig + r'''}
\authorcontributions{A.M., S.D., and K.J.S.Z. designed research; A.M. performed research and analyzed data; A.M., S.D., and K.J.S.Z. wrote the paper.}
\authordeclaration{The authors declare no competing interest.}
\correspondingauthor{\textsuperscript{1}To whom correspondence should be addressed. E-mail: amohseni@andrew.cmu.edu}
\keywords{science funding $|$ peer review $|$ lotteries $|$ Bayesian decision theory $|$ metascience}
\begin{abstract}
''' + abstract + r'''
\end{abstract}
\begin{document}
\maketitle
\thispagestyle{firststyle}
\ifthenelse{\boolean{shortarticle}}{\ifthenelse{\boolean{singlecolumn}}{\abscontentformatted}{\abscontent}}{}
''' + body + r'''
\matmethods{''' + methods + r'''}
\showmatmethods{}
\acknow{This work was supported by the John Templeton Foundation [grant number TODO].}
\showacknow{}
\bibliography{references}
\clearpage\onecolumn
\setboolean{shortarticle}{false}
''' + si + r'''
\end{document}
'''
open('pnas-typeset.tex', 'w').write(doc)
# supporting files
for f in glob.glob(f'{SRC}/*.dat') + [f'{SRC}/references.bib'] + glob.glob(f'{SRC}/si-*.tex') + glob.glob(f'{SRC}/figS*.tex') + glob.glob(f'{SRC}/tab-*.tex'):
    shutil.copy(f, '.')
for f in ['fig-p1', 'fig-p2', 'fig-p3', 'fig-p4', 'fig-regime']:
    t = open(f'{SRC}/{f}.tex').read().replace('\\textwidth', '\\linewidth').replace('\\ref{app:simulation}', 'SI Materials and Methods')
    if f == 'fig-p4':
        t = t.replace('\\begin{figure}[tp]', '\\begin{figure*}[t!]').replace('\\end{figure}', '\\end{figure*}')
        t = t.replace('width=0.76\\linewidth, height=5.4cm', 'width=0.5\\linewidth, height=4.2cm').replace('width=0.36\\linewidth, height=5.6cm', 'width=0.32\\linewidth, height=4.4cm')
    if f == 'fig-regime':
        t = t.replace('\\begin{figure}[t]', '\\begin{figure*}[t]').replace('\\end{figure}', '\\end{figure*}')
    if f == 'fig-p2':
        t = t.replace('height=6.0cm', 'height=4.6cm')
    if f == 'fig-p1':  # one-column figure: full column width, smaller label font
        t = t.replace('width=0.9\\linewidth, height=7.6cm,', 'width=\\linewidth, height=7.0cm,').replace('font=\\small', 'font=\\footnotesize')
    open(f'{f}.tex', 'w').write(t)
for cmd in ['pdflatex -interaction=nonstopmode pnas-typeset.tex', 'bibtex pnas-typeset', 'pdflatex -interaction=nonstopmode pnas-typeset.tex', 'pdflatex -interaction=nonstopmode pnas-typeset.tex']:
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
log = open('pnas-typeset.log').read()
print('errors', log.count('\n!'), 'undefined', log.lower().count('undefined'), 'pages', re.findall(r'\((\d+) pages', log))
