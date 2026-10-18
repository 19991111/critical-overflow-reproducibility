"""Prespecified fresh-micro intervals; reuse immutable old reference events."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from decimal import Context, Decimal, ROUND_CEILING
import csv,json,hashlib
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
REF=REPO/'research/paper_sde_end_to_end_20261017';GRID=100000;R=4096
def write(name,x):(ROOT/'results'/name).write_text(json.dumps(x,indent=2)+'\n')
def hx(s):
 m,e=s.lower().split('p');sign=-1 if m.startswith('-') else 1;m=m.lstrip('+-').removeprefix('0x');a,b=m.split('.');p=int(e)-4*len(b);return F(sign*int(a+b,16))*(2**p if p>=0 else F(1,2**(-p)))
def product(v):
 while len(v)>1:v=[v[j]*v[j+1] for j in range(0,len(v)-1,2)]+([v[-1]] if len(v)%2 else [])
 return v[0]
def interval(values):
 count=Counter(min(z,1536) for z in values);n=len(values)
 def reject(k,upper):
  den=(5 if upper else 1)*512*GRID;terms=[]
  for z,c in count.items():
   a=den+(2*(512*k-GRID*z) if upper else GRID*z-512*k);assert a>=0
   if a==0:return False
   terms.append(pow(a,c))
  return product(terms)>=1000*den**n
 l,u=0,GRID
 if not reject(0,False):L=0
 else:
  while u-l>1:
   k=(l+u)//2
   if reject(k,False):l=k
   else:u=k
  L=l
 l,u=GRID//2,3*GRID
 if reject(l,True):U=l
 elif not reject(u,True):U=u
 else:
  while u-l>1:
   k=(l+u)//2
   if reject(k,True):u=k
   else:l=k
  U=u
 return F(L,GRID),F(U,GRID)
def tail(U):
 gap=3-U
 if gap<=0:return None
 ctx=Context(prec=70,rounding=ROUND_CEILING);arg=-2*gap*gap;a=ctx.divide(Decimal(arg.numerator),Decimal(arg.denominator));e=ctx.next_plus(ctx.exp(a));rr=F(1,4)/gap;r=ctx.divide(Decimal(rr.numerator),Decimal(rr.denominator));v=F(ctx.multiply(r,e));return F((v.numerator*10**12+v.denominator-1)//v.denominator,10**12)
def directed(a,b):
 i=j=0;ab=ba=F(0);n,m=len(a),len(b)
 while i<n or j<m:
  x=b[j] if i==n else a[i] if j==m else min(a[i],b[j])
  while i<n and a[i]<=x:i+=1
  while j<m and b[j]<=x:j+=1
  d=F(i,n)-F(j,m);ab=max(ab,d);ba=max(ba,-d)
 return ab,ba
rows=list(csv.DictReader((ROOT/'results/micro/outcomes.csv').open()));assert len(rows)==16384
cells={}
for lam in [1,2]:
 for k in [0,1]:
  rr=sorted((x for x in rows if int(x['lambda'])==lam and int(x['kernel'])==k),key=lambda x:int(x['replicate']));assert len(rr)==R and [int(x['replicate']) for x in rr]==list(range(R));assert all(x['q']=='512' and x['n_completed']=='262144' for x in rr);cells[lam,k]=rr
refs={1:[],2:[]};readouts={1:[],2:[]};seen=set()
for x in csv.DictReader((REF/'results/reference_endpoints.csv').open()):
 lam,rid=int(x['lambda']),int(x['replicate']);assert (lam,rid) not in seen;seen.add((lam,rid));refs[lam].append(hx(x['D_hex']));readouts[lam].append(hx(x['S_hex']))
assert len(seen)==65536 and all((l,r) in seen for l in [1,2] for r in range(32768))
old=json.loads((REF/'logs/reference_means.stdout.log').read_text());assert old['lambda1']=={'lower_grid':69875,'upper_grid':70333} and old['lambda2']=={'lower_grid':67391,'upper_grid':67781}
quality=[];refci={}
for lam in [1,2]:
 r=F('.010') if lam==1 else F('.020');err=F('.0006') if lam==1 else F('.0010');L=F(old[f'lambda{lam}']['lower_grid'],GRID)-err;U=F(old[f'lambda{lam}']['upper_grid'],GRID)+err;refci[lam]=(L,U);refs[lam].sort();a=refs[lam];i=mx=0
 for j in range(len(a)):
  while a[j]-a[i]>2*r:i+=1
  mx=max(mx,j-i+1)
 width=F(mx,len(a))+F('.034');quality.append({'lambda':lam,'true_SDE_mean_interval':[str(L),str(U)],'CDF_envelope_width_upper':str(width),'mean_precision_passed':U-L<=F('.010'),'CDF_precision_passed':width<=F('.10')})
means=[];cdfs=[];micro=[]
for (lam,k),rr in cells.items():
 z=[int(x['Zcrit']) for x in rr];L,U=interval(z);mu=F(sum(z),512*R);coarse=mu+F('.030');beta=tail(coarse);upper=U+beta if beta is not None else F(512);gate=beta is not None and upper-L<=F('.04');rl,ru=refci[lam];lo=L-ru;hi=upper-rl;refgood=next(x['mean_precision_passed'] for x in quality if x['lambda']==lam);decision='supported' if lo>=-F('.05') and hi<=F('.05') else 'detected_mismatch' if lo>F('.05') or hi<-F('.05') else 'inconclusive'
 if not gate or not refgood:decision='inconclusive'
 micro.append({'lambda':lam,'kernel':k,'R':R,'observed_unconditional_mean':str(mu),'clipped_mean_lower':str(L),'clipped_mean_upper':str(U),'coarse_mean_upper':str(coarse),'tail_bias_upper':None if beta is None else str(beta),'unconditional_mean_interval':[str(L),str(upper)],'precision_gate_passed':gate,'clipped_observation_count':sum(v>1536 for v in z)})
 means.append({'lambda':lam,'kernel':k,'observed_micro_mean':str(mu),'true_mean_difference_interval':[str(lo),str(hi)],'absolute_margin':'1/20','decision':decision,'micro_precision_passed':gate,'reference_precision_passed':refgood})
 a=sorted(F(int(x['Dcrit_integer']),512) for x in rr);b=refs[lam];r=F('.010') if lam==1 else F('.020');ap,bp=directed(a,[x+r for x in b]);am,bm=directed(a,[x-r for x in b]);delta=F('.049');low=max(F(0),max(am,bp)-delta);up=min(F(1),max(ap,bm)+delta);ks=max(directed(a,b));refgood=next(x['CDF_precision_passed'] for x in quality if x['lambda']==lam);decision='supported' if up<=F('.20') else 'detected_mismatch' if low>F('.20') else 'inconclusive'
 if not refgood:decision='inconclusive'
 cdfs.append({'lambda':lam,'kernel':k,'empirical_KS_to_computed_reference':str(ks),'true_CDF_distance_interval':[str(low),str(up)],'margin':'1/5','horizontal_radius':str(r),'vertical_total':'49/1000','decision':decision,'reference_precision_passed':refgood})
write('MICRO_MEAN_INTERVALS.json',micro);write('REFERENCE_QUALITY_REUSED.json',quality);write('MEAN_COMPARISONS.json',means);write('CDF_COMPARISONS.json',cdfs)
summary={'status':'complete_fresh_micro_batch_analyzed','q':512,'kappa':1,'lambda_values':[1,2],'R_per_cell':R,'fresh_endpoint_rows':16384,'fresh_paired_runs':8192,'new_SDE_paths':0,'old_reference_rows_reused':65536,'mean_decisions':{d:sum(x['decision']==d for x in means) for d in ['supported','detected_mismatch','inconclusive']},'CDF_decisions':{d:sum(x['decision']==d for x in cdfs) for d in ['supported','detected_mismatch','inconclusive']},'all_precision_gates_passed':all(x['precision_gate_passed'] for x in micro) and all(x['mean_precision_passed'] and x['CDF_precision_passed'] for x in quality),'joint_coverage':'at least .95 for these eight comparisons, conditional on ideal within/across-stream iid words and audited reference arithmetic; unconditional over reused reference and fresh micro','reference_is_reused_not_independent_replication':True,'old_postselection_mean_claim_restored':False,'physical_CPU_models_added':0}
write('SUMMARY.json',summary);print(json.dumps(summary))
