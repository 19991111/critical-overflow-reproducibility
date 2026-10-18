"""Frozen new mean/CDF intervals, reusing old micro data without micro reruns."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from decimal import Context,Decimal,ROUND_CEILING
import csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];OLD=REPO/'research/paper_v4_release_20261017';GRID=100000
def write(name,x):(ROOT/'results'/name).write_text(json.dumps(x,indent=2)+'\n')
def product(v):
 while len(v)>1:v=[v[j]*v[j+1] for j in range(0,len(v)-1,2)]+([v[-1]] if len(v)%2 else [])
 return v[0]
def exact_interval(values):
 count=Counter(min(z,1536) for z in values);n=len(values)
 def reject(k,upper):
  den=(5 if upper else 1)*512*GRID
  terms=[]
  for z,c in count.items():
   a=den+(2*(512*k-GRID*z) if upper else GRID*z-512*k)
   assert a>=0
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
def tail_upper(U):
 gap=3-U
 if gap<=0:return None
 ctx=Context(prec=70,rounding=ROUND_CEILING);arg=-2*gap*gap;argup=ctx.divide(Decimal(arg.numerator),Decimal(arg.denominator));ex=ctx.next_plus(ctx.exp(argup));ratio=F(1,4)/gap;rr=ctx.divide(Decimal(ratio.numerator),Decimal(ratio.denominator));v=F(ctx.multiply(rr,ex));return F((v.numerator*10**12+v.denominator-1)//v.denominator,10**12)
rows=list(csv.DictReader((OLD/'results/validation/outcomes.csv').open()));assert len(rows)==9216
micro={}
for lam in [1,2]:
 for kernel in [0,1]:
  rr=sorted([r for r in rows if int(r['q'])==512 and int(r['lambda'])==lam and int(r['kernel'])==kernel],key=lambda r:int(r['replicate']));assert len(rr)==1024 and [int(r['replicate']) for r in rr]==list(range(1024));micro[lam,kernel]=rr
if sys.argv[1]=='--micro':
 result=[]
 for (lam,kernel),rr in micro.items():
  z=[int(r['Zcrit']) for r in rr];L,U=exact_interval(z);mean=F(sum(z),512*1024);coarse=mean+F('.059');tail=tail_upper(coarse);assert tail is not None;upper=U+tail
  result.append({'lambda':lam,'kernel':kernel,'R':1024,'observed_unconditional_mean':str(mean),'clipped_mean_lower':str(L),'clipped_mean_upper':str(U),'coarse_mean_upper':str(coarse),'tail_bias_upper':str(tail),'unconditional_mean_lower':str(L),'unconditional_mean_upper':str(upper),'precision_gate_width_0_04':upper-L<=F('.04')})
 write('MICRO_MEAN_INTERVALS.json',result);write('MICRO_PRECHECK.json',{'all_four_precision_gates_passed':all(r['precision_gate_width_0_04'] for r in result),'exact_integer_betting_inversions':True,'new_micro_trajectories':0,'source_sha256':hashlib.sha256((OLD/'results/validation/outcomes.csv').read_bytes()).hexdigest()});print(json.dumps(result));sys.exit(0)
assert sys.argv[1]=='--final'
def hx(s):
 mant,exp=s.lower().split('p');sg=-1 if mant.startswith('-') else 1;mant=mant.lstrip('+-').removeprefix('0x');a,b=mant.split('.');e=int(exp)-4*len(b);return F(sg*int(a+b,16))*(2**e if e>=0 else F(1,2**(-e)))
refs={1:[],2:[]};readouts={1:[],2:[]};seen=set()
for r in csv.DictReader((ROOT/'results/reference_endpoints.csv').open()):
 lam,rid=int(r['lambda']),int(r['replicate']);assert(lam,rid)not in seen;seen.add((lam,rid));refs[lam].append(hx(r['D_hex']));readouts[lam].append(hx(r['S_hex']))
assert len(seen)==65536 and all(len(v)==32768 for v in refs.values());assert all((lam,rid)in seen for lam in [1,2] for rid in range(32768))
raw=json.loads((ROOT/'logs/reference_means.stdout.log').read_text());microci=json.loads((ROOT/'results/MICRO_MEAN_INTERVALS.json').read_text());quality=[];referenceci={}
for lam in [1,2]:
 radius=F('.010') if lam==1 else F('.020');err=F('.0006') if lam==1 else F('.0010');ci=raw[f'lambda{lam}'];L=F(ci['lower_grid'],GRID)-err;U=F(ci['upper_grid'],GRID)+err;vals=sorted(refs[lam]);refs[lam]=vals;i=0;maxcount=0
 for j in range(len(vals)):
  while vals[j]-vals[i]>2*radius:i+=1
  maxcount=max(maxcount,j-i+1)
 envelope=F(maxcount,len(vals))+2*(F('.012')+F('.005'));goodmean=U-L<=F('.010');goodcdf=envelope<=F('.10');referenceci[lam]=(L,U)
 quality.append({'lambda':lam,'clipped_reference_mean_interval':[str(F(ci['lower_grid'],GRID)),str(F(ci['upper_grid'],GRID))],'true_SDE_unconditional_spill_mean_interval':[str(L),str(U)],'reference_approximation_allowance':str(err),'reference_CDF_envelope_width_upper':str(envelope),'mean_precision_passed':goodmean,'CDF_precision_passed':goodcdf})
write('REFERENCE_QUALITY.json',quality)
def directed(a,b):
 i=j=0;ab=ba=F(0);n,m=len(a),len(b)
 while i<n or j<m:
  x=b[j] if i==n else a[i] if j==m else min(a[i],b[j])
  while i<n and a[i]<=x:i+=1
  while j<m and b[j]<=x:j+=1
  d=F(i,n)-F(j,m);ab=max(ab,d);ba=max(ba,-d)
 return ab,ba
means=[];cdfs=[]
for (lam,kernel),rr in micro.items():
 mc=next(r for r in microci if r['lambda']==lam and r['kernel']==kernel);rl,ru=referenceci[lam];lo=F(mc['unconditional_mean_lower'])-ru;hi=F(mc['unconditional_mean_upper'])-rl;precision=next(x for x in quality if x['lambda']==lam)['mean_precision_passed']
 decision='supported' if lo>=-F('.05') and hi<=F('.05') else 'detected_mismatch' if lo>F('.05') or hi<-F('.05') else 'inconclusive'
 if not precision:decision='inconclusive'
 means.append({'lambda':lam,'kernel':kernel,'true_unconditional_mean_difference_interval':[str(lo),str(hi)],'absolute_margin':'1/20','decision':decision,'reference_precision_passed':precision})
 microd=sorted(F(int(r['Dcrit_integer']),512) for r in rr);ref=refs[lam];radius=F('.010') if lam==1 else F('.020');plus=[x+radius for x in ref];minus=[x-radius for x in ref];mplus,pplus=directed(microd,plus);mminus,pminus=directed(microd,minus);delta=F('.060')+F('.012')+F('.005');L=max(F(0),max(mminus,pplus)-delta);U=min(F(1),max(mplus,pminus)+delta);ks=max(directed(microd,ref));precision=next(x for x in quality if x['lambda']==lam)['CDF_precision_passed'];decision='supported' if U<=F('.20') else 'detected_mismatch' if L>F('.20') else 'inconclusive'
 if not precision:decision='inconclusive'
 cdfs.append({'lambda':lam,'kernel':kernel,'empirical_KS_to_computed_reference':str(ks),'true_CDF_distance_interval':[str(L),str(U)],'practical_margin':'1/5','horizontal_radius':str(radius),'micro_DKW':'3/50','reference_DKW':'3/250','vertical_coupling_allowance':'1/200','reference_precision_passed':precision,'decision':decision})
write('END_TO_END_MEAN_COMPARISONS.json',means);write('END_TO_END_D_CDF_COMPARISONS.json',cdfs)
summary={'status':'complete_followup_batch_analyzed','reference_endpoint_rows':65536,'new_micro_trajectories':0,'new_reference_replicates_per_lambda':32768,'reference_steps':16384,'new_primary_mean_decisions':{k:sum(r['decision']==k for r in means) for k in ['supported','detected_mismatch','inconclusive']},'new_primary_D_CDF_decisions':{k:sum(r['decision']==k for r in cdfs) for k in ['supported','detected_mismatch','inconclusive']},'all_reference_precision_gates_passed':all(x['mean_precision_passed'] and x['CDF_precision_passed'] for x in quality),'joint_confidence_at_least':'0.95 conditional on ideal independent replicate streams','scope':'q512,kappa1,lambda1/2,existing fixed/product micro kernels; finite practical tolerances only','old_primary_checks_or_manuscripts_rewritten':False}
write('SUMMARY.json',summary);print(json.dumps(summary));print(json.dumps(quality));print(json.dumps(means));print(json.dumps(cdfs))
