"""Generate the Appendix E / SI Text 5 tables and fig15-gap-gini.dat from the M = 1000 runs.
Writes tab-sens-*.tex into /home/claude/gf/tex and copies to /home/claude/gf/tex/pnas, and prints
the numbers quoted in prose."""
import csv, math, shutil
import numpy as np
D = '/home/claude/gf/tex/'
def load(f): return list(csv.DictReader(open(f)))
sens = load('sens_results.csv'); sens2 = load('sens2_results.csv'); sig = load('signals_results.csv')
def pick(rows, **kw):
    out = [r for r in rows if all(abs(float(r[k]) - v) < 1e-9 if isinstance(v, (int, float)) else r[k] == v for k, v in kw.items())]
    assert len(out) == 1, (kw, len(out)); return out[0]
f2 = lambda x: f"{float(x):.2f}"
def table(caption, label, header, body, colspec):
    L = ['\\begin{table}[t]', '\\centering', '\\small', '\\caption{%s}' % caption, '\\label{%s}' % label,
         '\\begin{tabular}{%s}' % colspec, '\\toprule'] + header + ['\\midrule'] + body + ['\\bottomrule', '\\end{tabular}', '\\end{table}', '']
    return '\n'.join(L)
KS = [1.3, 2.0, 3.5]
# ---- Table: mean capability (fixed scale 1 from sens/scale1; fixed mean from sens2/norm_fine)
body = []
for ks in KS:
    a = [pick(sens, label='scale1', k_shape=ks, tau_k=t) for t in (0.3, 1, 3)]
    b = [pick(sens, label='meanK2', k_shape=ks, tau_k=t) for t in (0.3, 1, 3)]
    body.append(f"{ks:g} & {float(a[0]['meanK']):.2f} & " + ' & '.join(f2(r['rev_share']) for r in a) + f" & {float(b[0]['meanK']):.2f} & " + ' & '.join(f2(r['rev_share']) for r in b) + ' \\\\')
t_mean = table("The value of review by field at a fixed Pareto scale and at a fixed mean capability. Share of the no-review shortfall recovered, at review noise $\\tau_K = 0.3$, $1$, and $3$; 100 simulated populations per entry; $M = 1{,}000$ particles. The mean capability shown is the realized mean over the populations.",
    'tab:sens-mean', ['& \\multicolumn{4}{c}{fixed scale 1} & \\multicolumn{4}{c}{fixed mean 2} \\\\', '$\\alpha_K$ & mean $K$ & $\\tau_K = 0.3$ & $1$ & $3$ & mean $K$ & $\\tau_K = 0.3$ & $1$ & $3$ \\\\'], body, '@{}lcccccccc@{}')
# ---- Table: matched informativeness (interpolate rev_share against sk_cor within norm_fine)
def matched(ks, target):
    rows = sorted([r for r in sens2 if r['label'] == 'norm_fine' and abs(float(r['k_shape']) - ks) < 1e-9], key=lambda r: float(r['sk_cor']) if 'sk_cor' in r else 0)
    return rows
# norm_fine in sens2 has no sk_cor column; use sens 'meanK2' label (same settings, M=1000 rerun of sens.R) which has sk_cor
def matched_val(ks, target):
    rows = sorted([r for r in sens if r['label'] == 'meanK2' and abs(float(r['k_shape']) - ks) < 1e-9], key=lambda r: float(r['sk_cor']))
    x = [float(r['sk_cor']) for r in rows]; y = [float(r['rev_share']) for r in rows]
    return float(np.interp(target, x, y))
body = [f"{ks:g} & " + ' & '.join(f2(matched_val(ks, c)) for c in (0.8, 0.55, 0.3)) + ' \\\\' for ks in KS]
t_matched = table("The value of review by field at equal informativeness. Share of the no-review shortfall recovered when the rank correlation between the review score and capability is held at the stated value; fixed mean capability; interpolated within the noise sweep of Table~\\ref{tab:sens-mean}.",
    'tab:sens-matched', ['$\\alpha_K$ & rank correlation 0.8 & 0.55 & 0.3 \\\\'], body, '@{}lccc@{}')
# ---- Table: rho, alpha_R, n (scale 1, tau_K = 1)
def cell(label, ks, **kw): return pick(sens, label=label, k_shape=ks, tau_k=1, **kw)
body = ['\\multicolumn{12}{@{}l}{\\emph{Value of targeting (share of uniform funding\'s gain)}} \\\\']
for ks in KS:
    v = [cell('rho', ks, rho=r)['targ_unif'] for r in (-0.5, 0, 0.5, 0.8)] + [cell('rshape', ks, r_shape=1.3)['targ_unif'], cell('rho', ks, rho=0)['targ_unif'], cell('rshape', ks, r_shape=3.5)['targ_unif']] + [cell('n', ks, n=25)['targ_unif'], cell('rho', ks, rho=0)['targ_unif'], cell('n', ks, n=100)['targ_unif'], cell('n', ks, n=200)['targ_unif']]
    body.append(f"{ks:g} & " + ' & '.join(f2(x) for x in v) + ' \\\\')
body.append('\\addlinespace'); body.append('\\multicolumn{12}{@{}l}{\\emph{Value of review (share of the no-review shortfall recovered)}} \\\\')
for ks in KS:
    v = [cell('rho', ks, rho=r)['rev_share'] for r in (-0.5, 0, 0.5, 0.8)] + [cell('rshape', ks, r_shape=1.3)['rev_share'], cell('rho', ks, rho=0)['rev_share'], cell('rshape', ks, r_shape=3.5)['rev_share']] + [cell('n', ks, n=25)['rev_share'], cell('rho', ks, rho=0)['rev_share'], cell('n', ks, n=100)['rev_share'], cell('n', ks, n=200)['rev_share']]
    body.append(f"{ks:g} & " + ' & '.join(f2(x) for x in v) + ' \\\\')
t_rest = table("The value of targeting (upper block) and of review (lower block) by field as the copula parameter $\\rho$ coupling capability and resources, the tail of the resource distribution, and the size of the applicant pool vary. Review noise $\\tau_K = 1$; Pareto scale 1; 50 simulated populations per entry; $M = 1{,}000$ particles. Default values: $\\rho = 0$, $\\alpha_R = 2$, $n = 50$.",
    'tab:sens-rest', ['& \\multicolumn{4}{c}{copula parameter $\\rho$} & \\multicolumn{3}{c}{resource tail $\\alpha_R$} & \\multicolumn{4}{c}{pool size $n$} \\\\', '$\\alpha_K$ & $-0.5$ & $0$ & $0.5$ & $0.8$ & $1.3$ & $2$ & $3.5$ & $25$ & $50$ & $100$ & $200$ \\\\'], body, '@{}lcccccccccccc@{}')
# ---- Table: contamination
body = []
for ks in KS:
    v = [pick(sens2, label='norm_fine', k_shape=ks, tau_k=0.3)['rev_share']] + [pick(sens2, label='contam', k_shape=ks, tau_k=0.3, beta=b)['rev_share'] for b in (0.25, 0.5, 1)] + \
        [pick(sens2, label='norm_fine', k_shape=ks, tau_k=1)['rev_share']] + [pick(sens2, label='contam', k_shape=ks, tau_k=1, beta=b)['rev_share'] for b in (0.25, 0.5, 1)]
    body.append(f"{ks:g} & " + ' & '.join(f2(x) for x in v) + ' \\\\')
t_contam = table("The value of review when the review score partly reflects resources. Share of the no-review shortfall recovered when the score is $K_i + \\beta R_{i0}$ plus noise, for $\\beta = 0$ (the main text's signal), $0.25$, $0.5$, and $1$; fixed mean capability; 50 simulated populations per entry (100 at $\\beta = 0$); $M = 1{,}000$ particles.",
    'tab:sens-contam', ['& \\multicolumn{4}{c}{$\\tau_K = 0.3$} & \\multicolumn{4}{c}{$\\tau_K = 1$} \\\\', '$\\alpha_K$ & $\\beta = 0$ & $0.25$ & $0.5$ & $1$ & $\\beta = 0$ & $0.25$ & $0.5$ & $1$ \\\\'], body, '@{}lcccccccc@{}')
# ---- Table: resource signal (tau_R sweep and none)
def sg(label, ks, par=None):
    rows = [r for r in sig if r['label'] == label and abs(float(r['k_shape']) - ks) < 1e-9 and (par is None or (r['par'] not in ('NA', '') and abs(float(r['par']) - par) < 1e-9))]
    assert len(rows) == 1, (label, ks, par, len(rows)); return rows[0]
body = ['\\multicolumn{7}{@{}l}{\\emph{No-review shortfall (share of the complete-information gain)}} \\\\']
for ks in KS:
    body.append(f"{ks:g} & " + ' & '.join(f2(sg('tauR', ks, t)['short4']) for t in (0.1, 0.3, 1, 3, 10)) + f" & {f2(sg('noRsignal', ks)['short4'])} \\\\")
body.append('\\addlinespace'); body.append('\\multicolumn{7}{@{}l}{\\emph{Value of review (share of the no-review shortfall recovered)}} \\\\')
for ks in KS:
    body.append(f"{ks:g} & " + ' & '.join(f2(sg('tauR', ks, t)['rev_share']) for t in (0.1, 0.3, 1, 3, 10)) + f" & {f2(sg('noRsignal', ks)['rev_share'])} \\\\")
t_tauR = table("The resource signal. The no-review shortfall (upper block) and the value of review (lower block) as the noise of the resource signal $\\tau_R$ varies from nearly exact to nearly uninformative, and with no resource signal; review noise $\\tau_K = 1$; fixed mean capability; 50 simulated populations per entry; $M = 1{,}000$ particles.",
    'tab:sens-tauR', ['& \\multicolumn{5}{c}{resource-signal noise $\\tau_R$} & none \\\\', '$\\alpha_K$ & $0.1$ & $0.3$ & $1$ & $3$ & $10$ & \\\\'], body, '@{}lcccccc@{}')
# ---- Table: multiplicative noise
body = ['\\multicolumn{6}{@{}l}{\\emph{Rank correlation between score and capability}} \\\\']
for ks in KS:
    body.append(f"{ks:g} & " + ' & '.join(f2(sg('mult', ks, s)['sk_cor']) for s in (0.1, 0.3, 0.6, 1, 1.5)) + ' \\\\')
body.append('\\addlinespace'); body.append('\\multicolumn{6}{@{}l}{\\emph{Value of review (share of the no-review shortfall recovered)}} \\\\')
for ks in KS:
    body.append(f"{ks:g} & " + ' & '.join(f2(sg('mult', ks, s)['rev_share']) for s in (0.1, 0.3, 0.6, 1, 1.5)) + ' \\\\')
t_mult = table("Multiplicative review noise. The review score is $K_i\\exp(e_i)$ with $e_i$ Gaussian of standard deviation $s$, and the funder knows the form; fixed mean capability; 50 simulated populations per entry; $M = 1{,}000$ particles. The upper block gives the informativeness of the score at each $s$.",
    'tab:sens-mult', ['& \\multicolumn{5}{c}{noise $s$} \\\\', '$\\alpha_K$ & $0.1$ & $0.3$ & $0.6$ & $1$ & $1.5$ \\\\'], body, '@{}lccccc@{}')
# ---- Table: fresh review each round (T = 2) with the once-observed comparison
body = []
for ks in KS:
    once = [pick(sens2, label='norm_fine', k_shape=ks, tau_k=t)['rev_share'] for t in (0.3, 1, 3)]
    fresh = [sg('fresh', ks, t)['rev_share'] for t in (0.3, 1, 3)]
    body.append(f"{ks:g} & " + ' & '.join(f2(x) for x in once) + ' & ' + ' & '.join(f2(x) for x in fresh) + ' \\\\')
t_fresh = table("Review observed once against review redrawn every round. Share of the no-review shortfall recovered over two rounds at review noise $\\tau_K = 0.3$, $1$, and $3$; fixed mean capability; 100 (once) and 50 (redrawn) simulated populations per entry; $M = 1{,}000$ particles.",
    'tab:sens-fresh', ['& \\multicolumn{3}{c}{observed once (main text)} & \\multicolumn{3}{c}{redrawn every round} \\\\', '$\\alpha_K$ & $\\tau_K = 0.3$ & $1$ & $3$ & $0.3$ & $1$ & $3$ \\\\'], body, '@{}lcccccc@{}')
for name, txt in [('tab-sens-mean', t_mean), ('tab-sens-matched', t_matched), ('tab-sens-rest', t_rest), ('tab-sens-contam', t_contam), ('tab-sens-tauR', t_tauR), ('tab-sens-mult', t_mult), ('tab-sens-fresh', t_fresh)]:
    open(D + name + '.tex', 'w').write(txt); shutil.copy(D + name + '.tex', D + 'pnas/' + name + '.tex')
# ---- fig15-gap-gini.dat and regressions over all sens settings
with open(D + 'fig15-gap-gini.dat', 'w') as f:
    f.write('gini targ share skcor a\n')
    for r in sens:
        f.write(f"{float(r['gap_gini']):.4f} {float(r['targ_unif']):.4f} {float(r['rev_share']):.4f} {float(r['sk_cor']):.3f} {float(r['k_shape']):g}\n")
shutil.copy(D + 'fig15-gap-gini.dat', D + 'pnas/fig15-gap-gini.dat')
g = np.array([float(r['gap_gini']) for r in sens]); tg = np.array([float(r['targ_unif']) for r in sens])
sh = np.array([float(r['rev_share']) for r in sens]); sk = np.array([float(r['sk_cor']) for r in sens])
mk = np.array([float(r['meanK']) for r in sens]); mr = np.array([float(r['meanR']) for r in sens]); a = np.array([float(r['k_shape']) for r in sens])
def r2(X, y):
    X = np.column_stack([np.ones(len(y))] + list(X)); b, *_ = np.linalg.lstsq(X, y, rcond=None); yh = X @ b
    return 1 - ((y - yh) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print('n settings', len(sens))
print('R2 log targ ~ gini', round(r2([g], np.log(tg)), 3), '+ log(meanK/meanR)', round(r2([g, np.log(mk / mr)], np.log(tg)), 3))
ok = (sh > 0.005) & (sh < 0.995) & (sk > 0.005) & (sk < 0.995)
lg = lambda p: np.log(p / (1 - p))
print('R2 logit share ~ gini + logit skcor', round(r2([g[ok], lg(sk[ok])], lg(sh[ok])), 3), '+ tail dummies', round(r2([g[ok], lg(sk[ok]), (a[ok] == 1.3) * 1.0, (a[ok] == 3.5) * 1.0], lg(sh[ok])), 3), 'n', ok.sum())
print('range gini', g.min().round(2), g.max().round(2))
