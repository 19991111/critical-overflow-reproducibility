"""Render saved measurements and frozen lookup values. No simulation or ODE."""
from pathlib import Path
import csv, hashlib, io, json, os, tarfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
import numpy as np

B=Path(__file__).resolve().parent
O=Path(os.environ.get('CRITICAL_SPILL_INPUT_ROOT',str(B.parent/'critical_spill_paper_20261018_v1')))
F=B/'manuscript/figures'
F.mkdir(exist_ok=True)
E=B/'verification/plot_data'
E.mkdir(exist_ok=True)
blue,orange,green='#0072B2','#D55E00','#009E73'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
 'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,
 'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,
 'axes.spines.right':False,'lines.linewidth':1.6})
def save(fig,name):
 # Scale line weights for the 5.125-inch SIAM text block. At final size the
 # thinnest illustration line is at least one point.
 from matplotlib.lines import Line2D
 from matplotlib.patches import Patch
 from matplotlib.text import Text
 for obj in fig.findobj():
  if isinstance(obj,(Line2D,Patch)):
   obj.set_linewidth(max(1.4,obj.get_linewidth()))
  if isinstance(obj,Text):
   obj.set_fontsize(obj.get_fontsize()*1.13)
 for ext in ['pdf','png','svg']:
  fig.savefig(F/(name+'.'+ext),dpi=240,bbox_inches='tight',pad_inches=.035)
 plt.close(fig)
def tarread(path,suffix):
 with tarfile.open(path) as t:
  members=[x for x in t.getmembers() if x.isfile() and x.name.endswith('/'+suffix)]
  assert len(members)==1,(path,suffix)
  return t.extractfile(members[0]).read()
def writecsv(name,rows):
 with (E/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
inputs=[O/'analysis/endpoint_comparison.csv',O/'inputs/sde_end_to_end_return.tar.gz',O/'inputs/paper_v4_release_return.tar.gz']
rows=list(csv.DictReader(inputs[0].open(encoding='utf-8')))
lines=tarread(inputs[1],'results/TABLES.hex').decode('ascii').splitlines()
ng,ns,density=map(int,lines[0].split());assert(ng,ns,density)==(4097,16385,512)
saved_s=np.array([float.fromhex(x) for x in lines[1+ng:]])

# Candidate graph diagram. Shapes denote entire sides, not literal bucket counts.
fig,axs=plt.subplots(1,2,figsize=(6.5,2.75),gridspec_kw={'width_ratios':[1.13,1]})
fig.subplots_adjust(wspace=.15)
ax=axs[0];ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax.text(.5,1.02,'(a) Candidate structure',ha='center',fontsize=10)
for y,label in [(.71,'0'),(.24,'1')]:
 ax.add_patch(FancyBboxPatch((.02,y-.155),.96,.35,boxstyle='round,pad=.01',fc='.97',ec='.68',lw=.8))
 ax.text(.5,y+.13,rf'Community $C_{label}$: $q$ buckets',ha='center',fontsize=8,bbox={'facecolor':'.97','edgecolor':'none','pad':.3})
 for x,side in [(.16,'L'),(.84,'R')]:
  for dx in [-.04,0,.04]:
   ax.add_patch(Circle((x+dx,y),.014,ec=blue,fc='white',lw=1))
  ax.text(x,y-.105,rf'${side}_{label}$ ($m$)',ha='center',fontsize=8)
 ax.annotate('',xy=(.73,y),xytext=(.27,y),arrowprops={'arrowstyle':'-','lw':1.5,'color':blue})
 ax.text(.5,y+.025,r'$(1-\varepsilon)/2$',ha='center',fontsize=8,color=blue)
for a,b in [((.20,.71),(.80,.24)),((.20,.24),(.80,.71))]:
 ax.annotate('',xy=b,xytext=a,zorder=2,arrowprops={'arrowstyle':'-','lw':1.3,'ls':'--','color':orange})
ax.text(.5,.465,r'$\varepsilon/2$ each',ha='center',va='center',fontsize=8,
 bbox={'facecolor':'white','edgecolor':'none','pad':1.2})
ax.text(.5,.0,'Choose the smaller load; left wins ties.',ha='center',fontsize=8)
ax=axs[1];ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
ax.text(.5,1.02,'(b) Critical description',ha='center',fontsize=10)
ax.text(.5,.87,r'$n_* = 2qc \sim \kappa q^2$',ha='center',fontsize=11)
ax.text(.5,.76,r'$\varepsilon q\to\lambda$; slow time $t=n/q^2$',ha='center',fontsize=8.5)
for y,text in [(.48,'Community imbalance\n'+r'$D_q(t)=\Delta_n/q\ \longrightarrow\ D_t$'),(.12,'Overflow near capacity\n'+r'$Z_c(n_*+\lfloor 2hq\rfloor)/q\ \longrightarrow\ S_h(D_\kappa)$')]:
 ax.add_patch(FancyBboxPatch((.025,y),.95,.24,boxstyle='round,pad=.012',fc='white',ec='.5',lw=.9))
 ax.text(.5,y+.12,text,ha='center',va='center',linespacing=1.7,fontsize=8.5)
ax.text(.5,.405,'phase averaged in drift',ha='center',fontsize=8,color=blue)
ax.text(.5,.015,'observation phase retained in readout',ha='center',fontsize=8,color=blue)
save(fig,'model_scales')

cell=[r for r in rows if r['batch']=='prospective' and r['q']=='512' and r['lambda']=='1' and r['kernel']=='matching']
assert len(cell)==4096
d=np.array([int(r['Dcrit_integer'])/512 for r in cell]);x=np.array([float(r['observed']) for r in cell])
p=np.array([float(r['profile_prediction']) for r in cell]);b=np.abs(d)/2
assert np.max(np.abs(saved_s[np.rint(np.abs(d)*512).astype(int)]-p))<1e-14
assert np.all(x>=b)
writecsv('readout_endpoints.csv',[{'D':di,'observed':xi,'imbalance_only':bi,'profile':pi} for di,xi,bi,pi in zip(d,x,b,p)])
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8));fig.subplots_adjust(wspace=.34,bottom=.19,top=.88)
ax=axs[0];ax.scatter(d,x,s=5,alpha=.20,color='.4',edgecolors='none',rasterized=True,label='Saved endpoints')
grid=np.linspace(-3,3,1537)
ax.plot(grid,np.interp(np.abs(grid),np.arange(ns)/density,saved_s),color=blue,label=r'Profile $S_0(D)$')
ax.plot(grid,np.abs(grid)/2,color=orange,ls='--',label=r'Imbalance $|D|/2$')
ax.set(xlim=(-3,3),ylim=(0,1.65),xlabel=r'Observed imbalance $D=\Delta/q$',ylabel=r'Overflow $Z_c/q$',title='(a) Overflow at critical fill')
ax.set_xticks([-3,-1.5,0,1.5,3]);ax.legend(loc='upper center',frameon=False,labelspacing=.2,handlelength=2.1)
ax.grid(axis='y',lw=.6,color='.9')
ax=axs[1];errors=[np.abs(x-b),np.abs(x-p)]
for e,color,style,label in zip(errors,[orange,blue],['--','-'],['Imbalance only','Full profile']):
 vals,counts=np.unique(e,return_counts=True)
 ax.step(np.r_[0,vals],np.r_[0,np.cumsum(counts)/len(e)],where='post',color=color,ls=style,label=label)
ax.set(xlim=(0,.72),ylim=(0,1.025),xlabel='Absolute overflow error',ylabel='Fraction of endpoints',title='(b) Readout errors')
ax.set_xticks([0,.2,.4,.6]);ax.grid(axis='y',color='.9',lw=.6)
ax.legend(loc='lower right',frameon=False)
save(fig,'internal_contribution')

original=list(csv.DictReader(io.StringIO(tarread(inputs[2],'results/validation/outcomes.csv').decode())))
assert len(original)==9216
summ=[]
for q in [64,128,256,512]:
 for lam in [1,2]:
  for kernel in [0,1]:
   rr=[r for r in original if int(r['q'])==q and int(r['lambda'])==lam and int(r['kernel'])==kernel]
   assert len(rr)=={64:128,128:128,256:256,512:1024}[q]
   summ.append({'q':q,'lambda':lam,'kernel':'matching' if kernel==0 else 'product','n':len(rr),
    'mean_feedback_RMS':np.mean([float(r['g_rmse']) for r in rr]),
    'mean_max_overflow_error':np.mean([float(r['max_spill_error']) for r in rr])})
writecsv('finite_scale_summary.csv',summ)
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8));fig.subplots_adjust(wspace=.33,bottom=.21,top=.86)
for ax,key,title,margin in zip(axs,['mean_feedback_RMS','mean_max_overflow_error'],['(a) Phase-specific feedback','(b) Critical-window overflow'],[.10,.05]):
 for lam,color,marker in [(1,blue,'o'),(2,orange,'s')]:
  for kern,style in [('matching','-'),('product','--')]:
   rr=[r for r in summ if r['lambda']==lam and r['kernel']==kern]
   ax.plot([r['q'] for r in rr],[r[key] for r in rr],color=color,marker=marker,ms=4,ls=style,label=rf'$\lambda={lam}$, {kern}')
 ax.axhline(margin,color='.4',ls=':',lw=1.2)
 ax.set(xscale='log',xticks=[64,128,256,512],xlabel='Buckets per community q',ylabel='Mean trajectory error',title=title)
 ax.set_xticklabels(['64','128','256','512']);ax.minorticks_off();ax.grid(axis='y',lw=.6,color='.9')
 ax.text(.97,.95,'Frozen margin',transform=ax.get_yaxis_transform(),ha='right',va='top',fontsize=7.5) if False else None
axs[0].legend(frameon=False,fontsize=7.4,loc='best')
save(fig,'finite_scale')

ref=list(csv.DictReader(io.StringIO(tarread(inputs[1],'results/reference_endpoints.csv').decode())))
assert len(ref)==65536
fig,axs=plt.subplots(2,2,figsize=(6.5,3.65),gridspec_kw={'height_ratios':[2,1]})
fig.subplots_adjust(hspace=.15,wspace=.32,left=.1,right=.985,bottom=.15,top=.9)
cdfrecords=[]
for j,lam in enumerate([1,2]):
 refvalues=np.array([float.fromhex(r['D_hex']) for r in ref if int(r['lambda'])==lam])
 assert len(refvalues)==32768
 refsorted=np.sort(refvalues)
 ax=axs[0,j];lo=axs[1,j]
 xx=np.linspace(-3.5,3.5,2501)
 yr=np.searchsorted(refsorted,xx,side='right')/len(refsorted)
 ax.plot(xx,yr,color='.2',lw=1.7,label='SDE reference',zorder=1)
 for kern,col,sty in [('matching',blue,'--'),('product',orange,':')]:
  micro=[r for r in rows if r['batch']=='prospective' and int(r['lambda'])==lam and r['kernel']==kern]
  v=np.sort(np.array([int(r['Dcrit_integer'])/512 for r in micro]));assert len(v)==4096
  ym=np.searchsorted(v,xx,side='right')/len(v)
  ax.plot(xx,ym,color=col,ls=sty,lw=1.3,label=kern.capitalize())
  # Difference is evaluated at all empirical jumps, with right-continuous CDFs.
  jumps=np.unique(np.r_[v,refsorted]);jm=np.searchsorted(v,jumps,side='right')/len(v)
  jr=np.searchsorted(refsorted,jumps,side='right')/len(refsorted)
  lo.step(jumps,jm-jr,where='post',color=col,ls=sty,lw=1)
  cdfrecords.extend({'lambda':lam,'kernel':kern,'D':di,'F_micro':fi,'F_reference':ri} for di,fi,ri in zip(xx,ym,yr))
 ax.set(xlim=(-3,3),ylim=(0,1.02),title=rf'$\lambda={lam}$',ylabel='Empirical CDF' if j==0 else '')
 ax.tick_params(labelbottom=False);ax.grid(axis='y',color='.92',lw=.5)
 lo.axhline(0,color='.5',lw=.8);lo.set(xlim=(-3,3),ylim=(-.035,.035),xlabel=r'Terminal imbalance $D$',ylabel=r'$F_{micro}-F_{ref}$' if j==0 else '')
 lo.set_yticks([-.02,0,.02]);lo.grid(axis='y',color='.92',lw=.5)
axs[0,0].legend(loc='upper left',frameon=False,fontsize=7.5)
writecsv('terminal_cdf_grid.csv',cdfrecords)
save(fig,'terminal_law')

receipt={'operation':'replot_preserved_outputs','new_trajectories':0,'new_ODE_SDE_integrations':0,
 'original_outcomes':len(original),'prospective_outcomes':16384,'reference_outcomes':len(ref),
 'readout_cell':{'n':len(cell),'MAE_imbalance':float(errors[0].mean()),'MAE_profile':float(errors[1].mean())},
 'source_hashes':{str(f.relative_to(O)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
 'outputs':{str(f.relative_to(B)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in list(F.iterdir())+list(E.iterdir())}}
(B/'verification/FIGURE_CHECKS.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','outputs']},indent=2))
