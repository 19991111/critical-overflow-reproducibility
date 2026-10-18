"""Independent read-only reconstruction of saved checkpoint identities and metrics."""
from pathlib import Path
from collections import defaultdict
import csv,json,math,time
ROOT=Path(__file__).resolve().parents[1];V=ROOT/'results/validation';t0=time.monotonic()
def hist(s):return {int(k):int(v) for k,v in (x.split(':') for x in s.split(';'))}
def pairs(s):return {(int(a),int(b)):int(n) for a,b,n in (x.split(':') for x in s.split(';'))}
acc=defaultdict(lambda:{'g2':0.,'p2':0.,'r2':0.,'maxspill':0.,'f':0,'w':0})
seen=set();count=0;identities=0;zero=None
for file in sorted(V.glob('observations_part*.csv')):
 for r in csv.DictReader(file.open()):
  q,lam,rep,kernel,n=[int(r[k]) for k in ['q','lambda','replicate','kernel','n']];key=(q,lam,rep,kernel);full=key+(n,)
  assert full not in seen;seen.add(full)
  hh=[hist(r[k]) for k in ['L0','R0','L1','R1']];assert all(sum(h.values())==q//2 and all(v>0 for v in h.values()) for h in hh)
  masses=[sum(k*v for h in hh[i:i+2] for k,v in h.items()) for i in [0,2]]
  assert masses==[int(r['M0']),int(r['M1'])] and sum(masses)==n
  spill=sum(max(k-q//2,0)*v for h in hh for k,v in h.items());assert spill==int(r['Z'])
  def H(l,rr):return sum(nl*nr for a,nl in l.items() for b,nr in rr.items() if a<=b)
  prod=H(hh[0],hh[3])-H(hh[2],hh[1]);assert prod==int(r['g_product_num']) and int(r['g_product_den'])==(q//2)**2
  if kernel==0:
   joint=[pairs(r['pairs01']),pairs(r['pairs10'])]
   for table,li,ri in [(joint[0],0,3),(joint[1],2,1)]:
    left=defaultdict(int);right=defaultdict(int)
    for (a,b),nn in table.items():left[a]+=nn;right[b]+=nn
    assert dict(left)==hh[li] and dict(right)==hh[ri]
   gnum=sum(v for (a,b),v in joint[0].items() if a<=b)-sum(v for (a,b),v in joint[1].items() if a<=b)
   assert gnum==int(r['g_kernel_num']) and int(r['g_kernel_den'])==q//2
  else:assert int(r['g_kernel_num'])==prod and int(r['g_kernel_den'])==(q//2)**2
  gw=float(r['g_wave']);gk=int(r['g_kernel_num'])/int(r['g_kernel_den']);gp=prod/(q//2)**2;a=acc[key]
  if int(r['feedback_id'])>=0:a['g2']+=(gk-gw)**2;a['p2']+=(gp-gw)**2;a['r2']+=(gk-gp)**2;a['f']+=1
  if int(r['window_id'])>=0:a['maxspill']=max(a['maxspill'],abs(spill/q-float(r['S_wave'])));a['w']+=1
  if lam==0:
   if kernel==0:zero=(q,rep,n,hh,spill)
   else:assert zero==(q,rep,n,hh,spill);identities+=1
  count+=1
maxdiff=0
for r in csv.DictReader((V/'outcomes.csv').open()):
 key=tuple(int(r[k]) for k in ['q','lambda','replicate','kernel']);a=acc.pop(key);assert a['f']==16 and a['w']==5
 for field,truth in [('g_rmse',math.sqrt(a['g2']/16)),('profile_g_rmse',math.sqrt(a['p2']/16)),('pair_g_rmse',math.sqrt(a['r2']/16)),('max_spill_error',a['maxspill'])]:
  diff=abs(float(r[field])-truth);maxdiff=max(maxdiff,diff);assert diff<5e-14,(key,field,diff)
assert not acc
report={'status':'passed','saved_observations':count,'unique_checkpoint_keys':len(seen),'zero_lambda_paired_checkpoint_identities':identities,'outcome_metrics_reconstructed':9216,'all_histogram_masses_spills_product_feedback_and_fixed_joint_tables_checked':True,'largest_metric_roundoff_difference':maxdiff,'tolerance':5e-14,'new_trajectories':0,'elapsed_seconds':time.monotonic()-t0}
(ROOT/'results/INDEPENDENT_OBSERVATION_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
