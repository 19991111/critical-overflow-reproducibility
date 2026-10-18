"""Draw the opening result figure from preserved endpoints and lookup nodes."""
from pathlib import Path
import csv, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

B = Path(__file__).resolve().parent
points = B / 'verification/plot_data/readout_endpoints.csv'
lookup = B / 'evidence/FROZEN_OVERFLOW_NODES.csv'
with points.open(encoding='utf-8') as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 4096
d, x, b, p = [np.array([float(r[k]) for r in rows])
              for k in ['D', 'observed', 'imbalance_only', 'profile']]
assert np.all(x >= b) and np.max(np.abs(b - np.abs(d) / 2)) < 1e-14
with lookup.open(encoding='utf-8') as f:
    table = list(csv.DictReader(f))
grid_d = np.array([float(r['d']) for r in table])
grid_s = np.array([float.fromhex(r['S0_hex']) for r in table])
assert len(table) == 16385 and grid_d[0] == 0
assert np.max(np.abs(np.interp(np.abs(d), grid_d, grid_s) - p)) < 1e-14
errors = [np.abs(x - b), np.abs(x - p)]
maes = [float(e.mean()) for e in errors]
blue, orange, green = '#0072B2', '#D55E00', '#00745A'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelsize': 10, 'axes.titlesize': 10.8, 'legend.fontsize': 8.3,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 1.4, 'xtick.major.width': 1.4,
    'ytick.major.width': 1.4, 'lines.linewidth': 1.8})
fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.7))
fig.subplots_adjust(left=.085, right=.985, bottom=.21, top=.87, wspace=.32)
ax = axs[0]
ax.scatter(d, x, s=5, alpha=.18, color='.42', edgecolors='none',
           rasterized=True, label='Saved endpoints')
xx = np.arange(-1536, 1537) / 512
ax.plot(xx, np.interp(np.abs(xx), grid_d, grid_s), color=blue,
        label=r'Internal profile $S_0(D)$')
ax.plot(xx, np.abs(xx) / 2, color=orange, ls='--', label=r'Totals only $|D|/2$')
ax.annotate('', xy=(0, grid_s[0]), xytext=(0, 0),
            arrowprops={'arrowstyle': '<->', 'color': green, 'lw': 1.5})
ax.text(.40, .10, r'$S_0(0)\approx0.6313$', color=green, fontsize=8.5)
ax.set(xlim=(-3, 3), ylim=(0, 1.68), xticks=[-3, -1.5, 0, 1.5, 3],
       xlabel=r'Observed imbalance $D=\Delta/q$', ylabel=r'Overflow $Z_c/q$',
       title='(a) Totals miss internal overflow')
ax.legend(frameon=False, loc='upper center', labelspacing=.2, handlelength=1.9)
ax.grid(axis='y', color='.92', lw=1.4)
ax = axs[1]
for e, color, style, name, mae in zip(errors, [orange, blue], ['--', '-'],
                                    ['Totals', 'Profile'], maes):
    vals, counts = np.unique(e, return_counts=True)
    ax.step(np.r_[0, vals], np.r_[0, np.cumsum(counts)/len(e)], where='post',
            color=color, ls=style, label=f'{name}: {mae:.4f} MAE')
ax.set(xlim=(0, .72), ylim=(0, 1.025), xticks=[0, .2, .4, .6],
       xlabel='Absolute overflow error', ylabel='Fraction of endpoints',
       title='(b) Error at the same endpoints')
ax.legend(frameon=False, loc='lower right', handlelength=1.9)
ax.grid(axis='y', color='.92', lw=1.4)
paths = []
for ext in ['pdf', 'svg', 'png']:
    out = B / f'manuscript/figures/readout_teaser.{ext}'
    fig.savefig(out, dpi=240, bbox_inches='tight', pad_inches=.04)
    paths.append(out)
plt.close(fig)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {'operation': 'replot_saved_readout_comparison', 'samples': len(rows),
    'MAE_totals_only': maes[0], 'MAE_profile': maes[1],
    'new_trajectories': 0, 'new_ODE_SDE_integrations': 0,
    'observed_imbalance_used_by_both_readouts': True,
    'new_confidence_claim': False,
    'input_sha256': {f.relative_to(B).as_posix(): sha(f) for f in [points, lookup]},
    'output_sha256': {f.relative_to(B).as_posix(): sha(f) for f in paths},
    'scope': 'Descriptive mechanism comparison; not a comparison with a published algorithm or an unconditional diffusion validation.'}
(B / 'verification/TEASER_FIGURE_CHECKS.json').write_text(
    json.dumps(receipt, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in receipt.items() if 'sha256' not in k}, indent=2))
