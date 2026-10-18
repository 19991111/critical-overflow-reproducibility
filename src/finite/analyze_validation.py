"""Frozen analysis: exact primary failure-rate intervals and DKW decisions."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial,sqrt
from decimal import Decimal,Context,ROUND_FLOOR,ROUND_CEILING
import json,csv,statistics,time,signal,struct
import numpy as np
ROOT=Path(__file__).resolve().parents[1];signal.alarm(115);started=time.monotonic();V=ROOT/'results/validation';plan=json.loads((ROOT/'FINITE_Q_VALIDATION_PLAN.v1.json').read_text())
rows=list(csv.DictReader((V/'outcomes.csv').open()));expected=sum(s['replicates'] for s in plan['scales'])*3*2;assert len(rows)==expected
groups={}
for r in rows:
 key=tuple(int(r[k]) for k in ['q','lambda','kernel']);groups.setdefault(key,[]).append(r);assert int(r['n_completed'])==key[0]**2+2*key[0] and r['feedback_count']=='16' and r['window_count']=='5'
for s in plan['scales']:
 for lam in plan['lambda']:
  for kernel in [0,1]:
   rr=groups[s['q'],lam,kernel];assert len(rr)==s['replicates'] and sorted(int(r['replicate']) for r in rr)==list(range(s['replicates']))
def binom_num(n,a,k,upper=False,den=1000):
 if a==0:return (den**n if (k<=0 if upper else k>=0) else 0)
 if a==den:return (den**n if (k<=n if upper else k>=n) else 0)
 b=den-a;t=b**n;z=0
 for j in range(n+1):
  if (j>=k if upper else j<=k):z+=t
  if j<n:t=t*(n-j)*a//((j+1)*b)
 return z
def cp(n,klo,khi):
 denominator=1000**n;delta=F(1,800)
 def lowertail(a):return F(binom_num(n,a,khi),denominator)
 def uppertail(a):return F(binom_num(n,a,klo,True),denominator)
 if klo==0:L=0
 else:
  l,u=0,1000
  while u-l>1:
   a=(l+u)//2
   if uppertail(a)<=delta:l=a
   else:u=a
  L=l
 if khi==n:U=1000
 else:
  l,u=0,1000
  while u-l>1:
   a=(l+u)//2
   if lowertail(a)<=delta:u=a
   else:l=a
  U=u
 return F(L,1000),F(U,1000)
descriptive=[];primary=[];dq=F('0.00004861');ds=F('0.00002440')
for key,rr in sorted(groups.items()):
 q,lam,kernel=key;R=len(rr);g=[float(r['g_rmse']) for r in rr];sp=[float(r['max_spill_error']) for r in rr];ys=[int(r['Zcrit'])/q for r in rr];dd=[int(r['Dcrit_integer'])/q for r in rr]
 descriptive.append({'q':q,'lambda':lam,'kernel':kernel,'R':R,'mean_spill_over_q':statistics.mean(ys),'spill_MCSE':statistics.stdev(ys)/sqrt(R),'mean_abs_D':statistics.mean(map(abs,dd)),'mean_D2':statistics.mean(x*x for x in dd),'mean_feedback_RMSE':statistics.mean(g),'mean_profile_response_RMSE':statistics.mean(float(r['profile_g_rmse']) for r in rr),'mean_pairing_response_RMSE':statistics.mean(float(r['pair_g_rmse']) for r in rr),'mean_window_max_error':statistics.mean(sp)})
 if q==512:
  metrics=[('spill_window','max_spill_error',F('0.05'),ds)]
  if lam:metrics.append(('feedback_RMSE','g_rmse',F('0.10'),dq))
  for label,field,margin,allowance in metrics:
   certain=sum(F(r[field])-allowance>margin for r in rr);possible=sum(F(r[field])+allowance>margin for r in rr);L,U=cp(R,certain,possible)
   decision='supported' if U<=F(1,10) else 'detected_mismatch' if L>F(1,10) else 'inconclusive'
   primary.append({'q':q,'lambda':lam,'kernel':kernel,'metric':label,'margin':str(margin),'numerical_allowance':str(allowance),'certain_failures':certain,'possible_failures':possible,'R':R,'failure_probability_lower':str(L),'failure_probability_upper':str(U),'decision':decision,'one_sided_alpha':'1/800'})
assert len(primary)==10
def ksnum(x,y):
 x=sorted(x);y=sorted(y);i=j=maximum=0
 for z in sorted(set(x+y)):
  while i<len(x) and x[i]<=z:i+=1
  while j<len(y) and y[j]<=z:j+=1
  maximum=max(maximum,abs(i-j))
 return maximum
x=2*1024*F(57,1000)**2;lower_exp=sum((x**k/factorial(k) for k in range(41)),F(0));assert lower_exp>640
universality=[];paired=[]
for q in [64,128,256,512]:
 for lam in [0,1,2]:
  a=sorted(groups[q,lam,0],key=lambda r:int(r['replicate']));b=sorted(groups[q,lam,1],key=lambda r:int(r['replicate']));diff=[(int(x['Zcrit'])-int(y['Zcrit']))/q for x,y in zip(a,b)]
  if lam==0:assert all(x['Dcrit_integer']==y['Dcrit_integer'] and x['Zcrit']==y['Zcrit'] and x['driver_hash']==y['driver_hash'] for x,y in zip(a,b))
  paired.append({'q':q,'lambda':lam,'R':len(a),'mean_fixed_minus_product_spill':statistics.mean(diff),'paired_MCSE':statistics.stdev(diff)/sqrt(len(diff)),'nominal_3MCSE_lower':statistics.mean(diff)-3*statistics.stdev(diff)/sqrt(len(diff)),'nominal_3MCSE_upper':statistics.mean(diff)+3*statistics.stdev(diff)/sqrt(len(diff)),'strict_coverage_claim':False})
  if q==512 and lam:
   for field in ['Dcrit_integer','Zcrit']:
    num=ksnum([int(r[field]) for r in a],[int(r[field]) for r in b]);obs=F(num,len(a));L=max(F(0),obs-F(114,1000));U=min(F(1),obs+F(114,1000));decision='supported' if U<=F(1,5) else 'detected_mismatch' if L>F(1,5) else 'inconclusive'
    universality.append({'lambda':lam,'quantity':'Dcrit/q' if field=='Dcrit_integer' else 'Zcrit/q','empirical_KS':str(obs),'true_CDF_distance_lower':str(L),'true_CDF_distance_upper':str(U),'margin':'1/5','decision':decision,'pair_band':'57/500','joint_family_alpha':'1/40'})
for name,obj in [('DESCRIPTIVE.json',descriptive),('PRIMARY_FAILURE_RATES.json',primary),('UNIVERSALITY_CDF.json',universality),('PAIRED_MEANS_NOMINAL.json',paired)]:
 (ROOT/'results'/name).write_text(json.dumps(obj,indent=2)+'\n')
# Check all labelled critical states against mass and spill identities independently.
count=0
with (V/'critical_configs.bin').open('rb') as f:
 assert f.read(7)==b'FQCFG1\0'
 while True:
  h=f.read(24)
  if not h:break
  assert len(h)==24;q,lam,rid,kernel,n=struct.unpack('<IIIIQ',h);raw=f.read(8*q);assert len(raw)==8*q;loads=struct.unpack('<'+'I'*(2*q),raw);assert sum(loads)==n==q*q
  rec=next(r for r in groups[q,lam,kernel] if int(r['replicate'])==rid);assert sum(max(a-q//2,0) for a in loads)==int(rec['Zcrit']);assert sum(loads[:q])-sum(loads[q:])==int(rec['Dcrit_integer']);count+=1
assert count==expected
# Deterministic lambda0 Brownian mean reference; no SDE Monte Carlo.
ctx=Context(prec=65);down=Context(prec=65,rounding=ROUND_FLOOR);up=Context(prec=65,rounding=ROUND_CEILING)
def dec(z,c):return c.divide(Decimal(z.numerator),Decimal(z.denominator))
def atan_bounds(x,N=80):
 s=sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(N)),F(0));return s,s+x**(2*N+1)/(2*N+1)
a0,a1=atan_bounds(F(1,5));b0,b1=atan_bounds(F(1,239));pilo=16*a0-4*b1;pihi=16*a1-4*b0;assert F(3)<pilo<pihi<F(22,7)
sqlo=ctx.next_minus(ctx.sqrt(down.multiply(Decimal(2),dec(pilo,down))));sqhi=ctx.next_plus(ctx.sqrt(up.multiply(Decimal(2),dec(pihi,up))));normlo=down.divide(Decimal(1),sqhi);normhi=up.divide(Decimal(1),sqlo)
table=[]
for line in (ROOT/'results/WAVE_TABLE.hex').read_text().splitlines()[1:]:table.append([float.fromhex(v) for v in line.split()])
table=np.array(table);levels=np.arange(-24,25);step=F(1,4096);sumlo=Decimal(0);sumhi=Decimal(0)
def sigma_vec(u):
 sh=np.floor(u).astype(int);ph=(u-sh)*256;j=ph.astype(int);t=ph-j;tails=(1-t[:,None])*table[j]+t[:,None]*table[j+1];p=-(np.diff(tails[:,:50],axis=1)+np.diff(tails[:,50:],axis=1));return np.sum(p*np.maximum(levels[None,:]+sh[:,None],0),axis=1)/2
for start in range(0,32768,512):
 dsamp=(np.arange(start,start+512)+.5)/4096;values=sigma_vec(dsamp/2)+sigma_vec(-dsamp/2)
 for d,v in zip(dsamp,values):
  dd=Decimal.from_float(float(d));arg=ctx.divide(ctx.multiply(dd,dd),Decimal(-2));ex=ctx.exp(arg);el=ctx.next_minus(ex);eh=ctx.next_plus(ex);sv=Decimal.from_float(float(v));assert sv>=0
  sumlo=down.add(sumlo,down.multiply(sv,down.multiply(normlo,el)));sumhi=up.add(sumhi,up.multiply(sv,up.multiply(normhi,eh)))
ilo=down.multiply(sumlo,dec(2*step,down));ihi=up.multiply(sumhi,dec(2*step,up));phi8_hi=up.multiply(normhi,ctx.next_plus(ctx.exp(Decimal(-32))));tail_hi=up.multiply(up.multiply(Decimal(2),phi8_hi),Decimal('1.079'))
err=step/2+F('8.632')*step+F('0.00002440')*(1+step);errhi=up.add(dec(err,up),tail_hi);reflo=down.subtract(ilo,errhi);refhi=up.add(ihi,errhi);assert errhi<Decimal('.005')
base=next(d for d in descriptive if d['q']==512 and d['lambda']==0 and d['kernel']==0);mean=base['mean_spill_over_q'];mc=base['spill_MCSE'];diaglo=mean-3*mc-float(refhi);diaghi=mean+3*mc-float(reflo)
reference={'status':'deterministic_Brownian_mean_reference_enclosed','lambda':0,'T':'1','cutoff':'8','midpoint_step':'1/4096','deterministic_error_upper':str(errhi),'mean_lower':str(reflo),'mean_upper':str(refhi),'wave_lookup_included':True,'normal_tail_upper':str(tail_hi),'new_SDE_samples':0,'largest_q_mean_spill':mean,'MCSE':mc,'mean_minus_reference_nominal_3MCSE_interval':[diaglo,diaghi],'nominal_interval_strict_coverage_claim':False,'eta':.05,'halfwidth_precision_gate':(diaghi-diaglo)/2<=.025,'nominal_interval_inside_eta':diaglo>=-.05 and diaghi<=.05}
(ROOT/'results/LAMBDA0_REFERENCE.json').write_text(json.dumps(reference,indent=2)+'\n')
summary={'status':'planned_batch_analyzed','outcome_rows':len(rows),'critical_configs_independently_checked':count,'primary_decisions':{k:sum(r['decision']==k for r in primary) for k in ['supported','detected_mismatch','inconclusive']},'universality_decisions':{k:sum(r['decision']==k for r in universality) for k in ['supported','detected_mismatch','inconclusive']},'lambda0_stepwise_identity_passed':True,'DKW_scalar_check_exact':True,'nonzero_lambda_SDE_reference_computed':False,'statistical_guarantees_conditional_on_ideal_independent_stream_reference':True,'elapsed_analysis_budget_seconds':time.monotonic()-started}
(ROOT/'results/VALIDATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));print(json.dumps(reference))
