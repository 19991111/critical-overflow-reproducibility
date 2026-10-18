"""Replot saved profile atoms and frozen response nodes. No integration."""
from pathlib import Path
import csv, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

B = Path(__file__).resolve().parent
F = B / 'manuscript/figures'
E = B / 'verification/plot_data'
profile_file = B / 'evidence/raw_profiles.npz'
response_file = B / 'evidence/FROZEN_GBAR_NODES.csv'
blue, orange = '#0072B2', '#D55E00'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelsize': 10, 'axes.titlesize': 11,
    'legend.fontsize': 8.5, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 1.4, 'xtick.major.width': 1.4,
    'ytick.major.width': 1.4, 'lines.linewidth': 1.7,
})
records = []
with np.load(profile_file) as r:
    levels = r['levels']
    assert np.array_equal(levels, np.arange(-24, 25))
    for row, phase in [(0, 0), (128, 0.5)]:
        x, y = r['x_128'][row], r['y_128'][row]
        for side, tail in [('L', x), ('R', y)]:
            atoms = tail[:-1] - tail[1:]
            assert atoms.min() >= 0 and abs(atoms.sum() - 1) < 1e-12
            records.extend({'phase': phase, 'side': side, 'level': int(k),
                            'probability': float(a)} for k, a in zip(levels, atoms))
with response_file.open(encoding='utf-8') as f:
    rows = list(csv.DictReader(f))
d = np.array([float(r['d']) for r in rows])
g = np.array([float.fromhex(r['gbar_hex']) for r in rows])
assert len(d) == 4097 and d[0] == g[0] == 0 and d[-1] == 8
assert np.all(np.diff(g) < 0) and np.all(g >= -1)

fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.65))
fig.subplots_adjust(left=.075, right=.985, bottom=.20, top=.87, wspace=.32)
ax = axs[0]
for phase, style in [(0, '-'), (0.5, '--')]:
    for side, color, marker in [('L', blue, 'o'), ('R', orange, 's')]:
        rr = [r for r in records if r['phase'] == phase and r['side'] == side]
        xx = np.array([r['level'] for r in rr])
        yy = np.array([r['probability'] for r in rr])
        keep = (xx >= -4) & (xx <= 4)
        ax.plot(xx[keep], yy[keep], color=color, ls=style, marker=marker,
                ms=3.8, label=rf'${side},\ u={"0" if phase == 0 else "1/2"}$')
ax.set(xlim=(-3.2, 3.2), ylim=(-.005, .70), xticks=[-3, -2, -1, 0, 1, 2, 3],
       xlabel='Integer load k', ylabel='Atom probability',
       title='(a) Two-side lattice profile')
ax.legend(frameon=False, ncol=2, loc='upper center', columnspacing=.85,
          handlelength=1.8, handletextpad=.35, labelspacing=.25)
ax.grid(axis='y', color='.91', lw=1.4)
ax = axs[1]
keep = d <= 4
xx = np.r_[-d[keep][:0:-1], d[keep]]
yy = np.r_[-g[keep][:0:-1], g[keep]]
ax.axhline(0, color='.7', lw=1.4)
ax.axvline(0, color='.7', lw=1.4)
ax.plot(xx, yy, color=blue, lw=1.8)
ax.set(xlim=(-4, 4), ylim=(-1.08, 1.08), xticks=[-4, -2, 0, 2, 4],
       yticks=[-1, -.5, 0, .5, 1], xlabel='Community imbalance d',
       ylabel=r'Averaged response $\overline{g}(d)$',
       title='(b) Restoring response')
ax.grid(axis='y', color='.93', lw=1.4)
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(F / f'profile_mechanism.{ext}', bbox_inches='tight',
                pad_inches=.04, dpi=240)
plt.close(fig)
with (E / 'profile_atoms.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['phase', 'side', 'level', 'probability'])
    w.writeheader()
    w.writerows(records)
with (E / 'averaged_response.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['d', 'gbar'])
    w.writerows(zip(xx, yy))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {
    'operation': 'replot_preserved_nodes', 'new_trajectories': 0,
    'new_ODE_SDE_integrations': 0, 'profile_rows': len(records),
    'response_nodes': len(xx), 'negative_response': 'Odd reflection of frozen nodes',
    'profile_phase_indices': [0, 128], 'profile_array': ['x_128', 'y_128'],
    'source_sha256': {p.relative_to(B).as_posix(): sha(p)
                      for p in [profile_file, response_file]},
    'outputs': {p.relative_to(B).as_posix(): {'bytes': p.stat().st_size, 'sha256': sha(p)}
                for p in [*(F / f'profile_mechanism.{e}' for e in ['pdf', 'svg', 'png']),
                          E / 'profile_atoms.csv', E / 'averaged_response.csv']},
    'scope': 'Display of saved nodes; connecting segments do not certify intermediate values.',
}
(B / 'verification/MECHANISM_FIGURE_CHECKS.json').write_text(
    json.dumps(receipt, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in receipt.items() if k not in ['outputs', 'source_sha256']}, indent=2))
