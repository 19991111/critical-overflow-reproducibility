"""Illustrate two exact model kernels; no experimental data or simulation."""
from pathlib import Path
import hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

B = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
    'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig = plt.figure(figsize=(6.5, 2.35))
for rect, title, matrix in [
    ([.06, .34, .30, .55], '(a) Fixed matching', np.eye(4)/4),
    ([.64, .34, .30, .55], '(b) Product sampling', np.ones((4,4))/16),
]:
    ax = fig.add_axes(rect)
    assert np.all(matrix.sum(axis=0) == .25)
    assert np.all(matrix.sum(axis=1) == .25)
    ax.imshow(matrix, vmin=0, vmax=.25, cmap='Blues', interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which='minor', color='white', lw=1.5)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_title(title, fontsize=10.5, pad=6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.4); spine.set_edgecolor('.5')
fig.text(.5, .61, 'Same uniform\nendpoint marginals', ha='center', va='center',
         fontsize=10, linespacing=1.6)
box = FancyBboxPatch((.23,.055), .54,.15, transform=fig.transFigure,
                    boxstyle='round,pad=.012', fc='#EFF5F8', ec='#0072B2', lw=1.4)
fig.patches.append(box)
fig.text(.5,.13,'A common leading critical law',ha='center',va='center',fontsize=11)
for start, end in [((.21,.32),(.37,.22)),((.79,.32),(.63,.22))]:
    fig.patches.append(FancyArrowPatch(start,end, transform=fig.transFigure,
                       arrowstyle='-|>',mutation_scale=12,lw=1.4,color='#0072B2'))
paths=[]
for ext in ['pdf','svg','png']:
    out=B/f'manuscript/figures/kernel_universality.{ext}'
    fig.savefig(out,dpi=240,bbox_inches='tight',pad_inches=.04);paths.append(out)
plt.close(fig)
receipt={'operation':'exact_model_illustration','example_side_size':4,
    'conditional_pair_laws':['I/4','J/16'],'row_and_column_marginals':.25,
    'scope':'Illustration of two fixed admissible kernels. The common law also requires the main model and scaling assumptions.',
    'new_trajectories':0,'new_ODE_SDE_integrations':0,
    'outputs':{p.relative_to(B).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
(B/'verification/KERNEL_FIGURE_CHECKS.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print('Exact model illustration completed.')
