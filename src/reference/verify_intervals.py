"""Read-only independent endpoint/CDF checks; no reference trajectory is generated."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Context,Decimal,ROUND_FLOOR,ROUND_CEILING
from bisect import bisect_right
import csv,hashlib,json,math,resource,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];OLD=REPO/'research/paper_v4_release_20261017';started=time.monotonic();cpu0=time.process_time();oscore=json.loads((ROOT/'protocol/ENVIRONMENT.json').read_text())['selected_single_core']
import os
os.sched_setaffinity(0,{oscore});resource.setrlimit(resource.RLIMIT_CPU,(180,181));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
def hx(s):
 m,e=s.lower().split('p');negative=m.startswith('-');m=m.lstrip('+-')[2:];a,b=m.split('.');power=int(e)-4*len(b);v=F((-1 if negative else 1)*int(a+b,16));return v*F(2**power) if power>=0 else v/F(2**(-power))
low=Context(prec=50,rounding=ROUND_FLOOR);high=Context(prec=50,rounding=ROUND_CEILING);near=Context(prec=50)
def logs(x):
 a=low.divide(Decimal(x.numerator),Decimal(x.denominator));b=high.divide(Decimal(x.numerator),Decimal(x.denominator));return low.next_minus(near.ln(a)),high.next_plus(near.ln(b))
threshold_lo,threshold_hi=logs(F(1000));refs={1:[],2:[]};Ds={1:[],2:[]}
for r in csv.DictReader((ROOT/'results/reference_endpoints.csv').open()):lam=int(r['lambda']);refs[lam].append(min(F(3),hx(r['S_hex'])));Ds[lam].append(hx(r['D_hex']))
intervals=json.loads((ROOT/'logs/reference_means.stdout.log').read_text());checks=[]
for lam in [1,2]:
 for upper,field,neighbor in [(False,'lower_grid',1),(True,'upper_grid',-1)]:
  k=intervals[f'lambda{lam}'][field]
  for kk in [k,k+neighbor]:
   m=F(kk,100000);lo=Decimal(0);hi=Decimal(0)
   for x in refs[lam]:
    factor=1+(m-x if upper else x-m)/8;a,b=logs(factor);lo=low.add(lo,a);hi=high.add(hi,b)
   endpoint=kk==k
   assert lo>threshold_hi if endpoint else hi<threshold_lo
   checks.append({'lambda':lam,'side':'upper' if upper else 'lower','candidate_grid':kk,'endpoint':endpoint,'log_capital_lower':str(lo),'log_capital_upper':str(hi),'strictly_correct_side_of_log1000':True})
rows=list(csv.DictReader((OLD/'results/validation/outcomes.csv').open()));micro=json.loads((ROOT/'results/MICRO_MEAN_INTERVALS.json').read_text());micro_checks=0
for ci in micro:
 rr=[r for r in rows if r['q']=='512' and int(r['lambda'])==ci['lambda'] and int(r['kernel'])==ci['kernel']];zz=[min(int(r['Zcrit']),1536) for r in rr];assert len(zz)==1024
 for upper,field,offset in [(False,'clipped_mean_lower',1),(True,'clipped_mean_upper',-1)]:
  k=int(F(ci[field])*100000)
  for kk in [k,k+offset]:
   d=(5 if upper else 1)*512*100000;terms=[d+(2*(512*kk-100000*z) if upper else 100000*z-512*kk) for z in zz];assert all(x>=0 for x in terms);left=math.prod(terms);right=1000*d**1024;assert (left>=right) if kk==k else (left<right);micro_checks+=1
def supremum(a,b):
 n,m=len(a),len(b);maximum=F(0)
 for x in sorted(set(a+b)):
  maximum=max(maximum,F(bisect_right(a,x),n)-F(bisect_right(b,x),m))
 return maximum
cdf=json.loads((ROOT/'results/END_TO_END_D_CDF_COMPARISONS.json').read_text());cdf_checks=0
for r in cdf:
 lam,k=r['lambda'],r['kernel'];a=sorted(F(int(x['Dcrit_integer']),512) for x in rows if x['q']=='512' and int(x['lambda'])==lam and int(x['kernel'])==k);b=sorted(Ds[lam]);rad=F(r['horizontal_radius']);plus=[x+rad for x in b];minus=[x-rad for x in b];allow=F('.077');L=max(F(0),max(supremum(a,minus),supremum(plus,a))-allow);U=min(F(1),max(supremum(a,plus),supremum(minus,a))+allow);assert [str(L),str(U)]==r['true_CDF_distance_interval'];cdf_checks+=1
report={'status':'passed','reference_betting_endpoint_and_neighbor_directed_Decimal_checks':checks,'micro_betting_endpoint_neighbor_exact_product_checks':micro_checks,'CDF_global_shift_suprema_independently_recomputed_by_bisection':cdf_checks,'new_reference_or_micro_trajectories':0,'CPU_seconds':time.process_time()-cpu0,'wall_seconds':time.monotonic()-started,'source_and_reference_SHA256':{'code':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference_endpoints':hashlib.sha256((ROOT/'results/reference_endpoints.csv').read_bytes()).hexdigest()}}
(ROOT/'results/INDEPENDENT_INTERVAL_CHECK.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='reference_betting_endpoint_and_neighbor_directed_Decimal_checks'}))
