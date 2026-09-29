"""Table tab:dynopt from dynopt2_results.csv (T = 2 exact optimum) and dynopt5_results.csv
(T = 5 division of the budget), and the DYNLOSS / DYNCOR placeholders in the texts."""
import csv, shutil, re
D = '/home/claude/gf/tex/'
r2 = list(csv.DictReader(open('dynopt2_results.csv'))); r5 = list(csv.DictReader(open('dynopt5_results.csv')))
def find(rows, aK, eps, b):
    for r in rows:
        if abs(float(r['aK']) - aK) < 1e-9 and abs(float(r['eps']) - eps) < 1e-9 and abs(float(r['b']) - b) < 1e-9: return r
L = ['\\begin{table}[tp]', '\\centering', '\\footnotesize', '\\setlength{\\tabcolsep}{3pt}',
     '\\caption{The round-by-round gap rule against the sequential optimum under complete information. Two rounds: the output lost by applying the gap rule to half the budget in each round, as a percent of the research output that the exact two-round optimum adds over no funding (mean and maximum over 20 populations); the correlation between the two allocations\' first-round grants (mean and minimum); and the share of the budget the optimum spends in round one. Five rounds: the output lost by an even division of the budget across rounds, as a percent of what the best division adds, with the gap rule choosing recipients in every round, and the center of mass of the best division (0.5 is even; higher is later). Mean capability 2; $n = 50$.}',
     '\\label{tab:dynopt}',
     '\\begin{tabular}{@{}cccccccc@{}}', '\\toprule',
     '& & & \\multicolumn{3}{c}{two rounds, exact optimum} & \\multicolumn{2}{c}{five rounds, division only} \\\\',
     '$\\alpha_K$ & $\\epsilon$ & $b$ & \\makecell{loss, \\%\\\\(max)} & \\makecell{correlation\\\\(min)} & \\makecell{round-1\\\\share} & \\makecell{even division\\\\loses, \\%} & \\makecell{center\\\\of mass} \\\\', '\\midrule']
for aK in (1.3, 2.0, 3.5):
    for eps in (0.1, 0.5, 0.85):
        for b in (0.2, 1.0, 2.0):
            a = find(r2, aK, eps, b); f = find(r5, aK, eps, b)
            L.append(f"{aK:g} & {eps:g} & {b:g} & {100*float(a['loss_share']):.1f} ({100*float(a['loss_share_max']):.1f}) & {float(a['cor_g1']):.3f} ({float(a['cor_g1_min']):.3f}) & {float(a['share_r1']):.2f} & {100*float(f['loss_share_even']):.1f} & {float(f['com']):.2f} \\\\")
L += ['\\bottomrule', '\\end{tabular}', '\\end{table}', '']
open(D + 'tab-dynopt.tex', 'w').write('\n'.join(L)); shutil.copy(D + 'tab-dynopt.tex', D + 'pnas/tab-dynopt.tex')
maxloss = max(float(r['loss_share_max']) for r in r2); meanloss = max(float(r['loss_share']) for r in r2)
mincor = min(float(r['cor_g1_min']) for r in r2); shares = [float(r['share_r1']) for r in r2]
loss5 = max(float(r['loss_share_even']) for r in r5); loss5_default = max(float(r['loss_share_even']) for r in r5 if abs(float(r['eps']) - 0.1) < 1e-9)
print('max loss (any population)', round(100*maxloss, 2), 'max mean loss', round(100*meanloss, 2), 'min cor', round(mincor, 3), 'share range', round(min(shares), 2), round(max(shares), 2), 'T5 max loss', round(100*loss5, 1), 'T5 at eps 0.1', round(100*loss5_default, 2))
DYNLOSS = f"{100*maxloss:.1f}".rstrip('0').rstrip('.'); DYNCOR = f"{mincor:.2f}"
for p in ['main.tex', 'pnas/pnas-main.tex', 'appendix-gap-rule.tex']:
    s = open(D + p).read(); s = s.replace('DYNLOSS', DYNLOSS).replace('DYNCOR', DYNCOR); open(D + p, 'w').write(s)
print('filled', DYNLOSS, DYNCOR)
