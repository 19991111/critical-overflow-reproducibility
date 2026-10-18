"""Independent raw-configuration, exact-capital and CDF audit; no trajectories."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Context,Decimal,ROUND_CEILING
from bisect import bisect_right
import csv,hashlib,json,math,struct,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];REF=REPO/'research/paper_sde_end_to_end_20261017';started=time.monotonic()
rows=list(csv.DictReader((ROOT/'results/micro/outcomes.csv').open()));indexed={(int(x['lambda']),int(x['replicate']),int(x['kernel'])):x for x in rows};assert len(indexed)==len(rows)==16384
seen=set();record=struct.Struct('<IIIIQ1024I')
files=sorted((ROOT/'results/micro').glob('terminal_configs_*.bin'));assert len(files)==8
for p in files:
 with p.open('rb') as f:
  assert f.read(7)==b'PMCFG1\0';nrecords=0
  while b:=f.read(record.size):
   assert len(b)==record.size;q,lam,rid,k,n,*a=record.unpack(b);key=(lam,rid,k);assert key in indexed and key not in seen;seen.add(key);r=indexed[key]
   assert q==512 and n==262144 and sum(a)==n and sum(a[:512])-sum(a[512:])==int(r['Dcrit_integer']) and sum(max(x-256,0) for x in a)==int(r['Zcrit']);nrecords+=1
  assert nrecords==2048
assert len(seen)==16384
progress=list(csv.DictReader((ROOT/'results/micro/progress.csv').open()));assert len(progress)==8192
for rid in range(4096):
 hashes={indexed[lam,rid,k]['driver_hash'] for lam in [1,2] for k in [0,1]};assert len(hashes)==1
for x in progress:
 assert int(x['arrivals'])==262144 and x['driver_hash']==indexed[int(x['lambda']),int(x['replicate']),0]['driver_hash']
ci=json.loads((ROOT/'results/MICRO_MEAN_INTERVALS.json').read_text());endpoint_checks=0
for x in ci:
 rr=[indexed[x['lambda'],j,x['kernel']] for j in range(4096)];zz=[min(int(r['Zcrit']),1536) for r in rr]
 for upper,field,offset in [(False,'clipped_mean_lower',1),(True,'clipped_mean_upper',-1)]:
  k=int(F(x[field])*100000)
  assert 0<k<100000 if not upper else 50000<k<300000
  for kk in [k,k+offset]:
   d=(5 if upper else 1)*512*100000;terms=[d+(2*(512*kk-100000*z) if upper else 100000*z-512*kk) for z in zz];assert all(t>=0 for t in terms);left=math.prod(terms);right=1000*d**4096;assert (left>=right) if kk==k else (left<right);endpoint_checks+=1
 mean=F(sum(int(r['Zcrit']) for r in rr),512*4096);assert str(mean)==x['observed_unconditional_mean'];coarse=mean+F('.030');assert str(coarse)==x['coarse_mean_upper'];gap=3-coarse;ctx=Context(prec=100,rounding=ROUND_CEILING);d=ctx.divide(Decimal(gap.numerator),Decimal(gap.denominator));arg=-2*gap*gap;v=ctx.multiply(ctx.divide(Decimal(1),ctx.multiply(Decimal(4),d)),ctx.next_plus(ctx.exp(ctx.divide(Decimal(arg.numerator),Decimal(arg.denominator)))));assert F(x['tail_bias_upper'])>=F(v)
 assert [x['clipped_mean_lower'],str(F(x['clipped_mean_upper'])+F(x['tail_bias_upper']))]==x['unconditional_mean_interval']
def parse(s):
 mant,expo=s.split('p');sign=-1 if mant.startswith('-') else 1;mant=mant.lstrip('+-')[2:];left,right=mant.split('.');p=int(expo)-4*len(right);value=F(sign*int(left+right,16));return value*2**p if p>=0 else value/2**(-p)
ref={1:[],2:[]}
for r in csv.DictReader((REF/'results/reference_endpoints.csv').open()):ref[int(r['lambda'])].append(parse(r['D_hex']))
for a in ref.values():assert len(a)==32768;a.sort()
def sup(a,b):
 result=F(0)
 for v in sorted(set(a+b)):result=max(result,F(bisect_right(a,v),len(a))-F(bisect_right(b,v),len(b)))
 return result
cdf_checks=0
for x in json.loads((ROOT/'results/CDF_COMPARISONS.json').read_text()):
 lam,k=x['lambda'],x['kernel'];a=sorted(F(int(indexed[lam,j,k]['Dcrit_integer']),512) for j in range(4096));b=ref[lam];r=F(x['horizontal_radius']);plus=[v+r for v in b];minus=[v-r for v in b];L=max(F(0),max(sup(a,minus),sup(plus,a))-F('.049'));U=min(F(1),max(sup(a,plus),sup(minus,a))+F('.049'));assert [str(L),str(U)]==x['true_CDF_distance_interval'];assert str(max(sup(a,b),sup(b,a)))==x['empirical_KS_to_computed_reference'];cdf_checks+=1
means=json.loads((ROOT/'results/MEAN_COMPARISONS.json').read_text());cdf=json.loads((ROOT/'results/CDF_COMPARISONS.json').read_text())
for x in means:
 m=next(y for y in ci if (y['lambda'],y['kernel'])==(x['lambda'],x['kernel']));rl,ru=([F('.69815'),F('.70393')] if x['lambda']==1 else [F('.67291'),F('.67881')]);L,U=F(m['unconditional_mean_interval'][0])-ru,F(m['unconditional_mean_interval'][1])-rl;assert [str(L),str(U)]==x['true_mean_difference_interval'];dec='supported' if -F('.05')<=L and U<=F('.05') else 'detected_mismatch' if L>F('.05') or U<-F('.05') else 'inconclusive'
 if not m['precision_gate_passed'] or not x['reference_precision_passed']:dec='inconclusive'
 assert dec==x['decision']
for x in cdf:
 L,U=map(F,x['true_CDF_distance_interval']);dec='supported' if U<=F('.20') else 'detected_mismatch' if L>F('.20') else 'inconclusive'
 if not x['reference_precision_passed']:dec='inconclusive'
 assert dec==x['decision']
result={'status':'passed','terminal_configurations_independently_checked':len(seen),'progress_records_checked':len(progress),'paired_and_cross_lambda_driver_hashes_checked':4096,'exact_endpoint_and_neighbor_products_checked':endpoint_checks,'CDF_comparisons_independently_checked_by_bisection':cdf_checks,'clipping_tail_bounds_checked_at_100_digits':4,'mean_difference_and_all_eight_decisions_checked':True,'new_random_words_or_trajectories':0,'wall_seconds':time.monotonic()-started}
(ROOT/'results/INDEPENDENT_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
