"""Build the contour data for the regime map (Fig. 4 PNAS / Fig. 16 long version) from
regime_results.csv. Axes: x = log10 b, y = Gini of the capability distribution 1/(2 alpha - 1).
Writes fig16-targ.dat and fig16-rev.dat (one block per contour segment, blank-line separated,
with a 'level' column), plus fig16-grid.dat with the raw cell values."""
import csv, math, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

rows = list(csv.DictReader(open('regime_results.csv')))
bs = sorted(set(float(r['b']) for r in rows)); ks = sorted(set(float(r['k_shape']) for r in rows))
gini = lambda a: 1 / (2 * a - 1)
ys = [gini(a) for a in ks][::-1]          # increasing Gini (decreasing alpha)
ks_rev = ks[::-1]
def grid(col):
    G = np.zeros((len(ys), len(bs)))
    for r in rows:
        i = ks_rev.index(float(r['k_shape'])); j = bs.index(float(r['b']))
        G[i, j] = float(r[col])
    return G
X = np.log10(bs)
T = grid('targ_unif'); Rv = grid('rev_gain'); Rg = grid('rev_share'); S4 = grid('short4')
out = open('fig16-grid.dat', 'w'); out.write('logb b alpha gini targ rev revgain short4\n')
for i, a in enumerate(ks_rev):
    for j, b in enumerate(bs):
        out.write(f"{X[j]:.4f} {b:g} {a:g} {ys[i]:.4f} {T[i,j]:.4f} {Rv[i,j]:.4f} {Rg[i,j]:.4f} {S4[i,j]:.4f}\n")
out.close()
# fine grid, linear interpolation in (log b, gini); log-transform targeting for smoother contours
xf = np.linspace(X[0], X[-1], 241); yf = np.linspace(ys[0], ys[-1], 241)
XX, YY = np.meshgrid(xf, yf)
def interp(G, log=False):
    Z = np.log(np.maximum(G, 1e-6)) if log else G
    f = RegularGridInterpolator((ys, X), Z, method='linear')
    Zf = f(np.stack([YY.ravel(), XX.ravel()], -1)).reshape(XX.shape)
    return np.exp(Zf) if log else Zf
Tf = interp(T, log=True); Rf = interp(Rv)
def write_contours(Z, levels, fname):
    fig, ax = plt.subplots(); cs = ax.contour(XX, YY, Z, levels=levels)
    with open(fname, 'w') as f:
        f.write('x y level\n')
        for lev, segs in zip(cs.levels, cs.allsegs):
            for seg in segs:
                for x, y in seg: f.write(f"{x:.4f} {y:.4f} {lev:g}\n")
                f.write('\n')
    plt.close(fig)
    return {lev: [s for s in segs] for lev, segs in zip(cs.levels, cs.allsegs)}
lev_t = [float(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.1, 0.25, 0.5, 1, 2]
lev_r = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0.2, 0.4, 0.6, 0.8]
ct = write_contours(Tf, lev_t, 'fig16-targ.dat'); cr = write_contours(Rf, lev_r, 'fig16-rev.dat')
print('targ range', T.min().round(3), T.max().round(3), 'rev range', Rv.min().round(3), Rv.max().round(3))
for lev, segs in ct.items():
    for s in segs: print('targ', lev, 'ends', s[0].round(2), s[-1].round(2))
for lev, segs in cr.items():
    for s in segs: print('rev', lev, 'ends', s[0].round(2), s[-1].round(2))
