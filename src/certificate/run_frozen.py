"""Frozen stage sequence, directed arithmetic, fixed budget and unchanged acceptance."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,subprocess,sys,os,time,csv
from bounds import ROOT,OLD,S,inputs,wave_tail,candidate_tail,check_anchors,ceil_scaled
from math import factorial
for p,h in json.loads((ROOT/'protocol/SOURCE_FREEZE.json').read_text())['files'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
plan=json.loads((ROOT/'protocol/PLAN_V1.json').read_text());anchor=check_anchors();data=inputs();best=None;completed=[];skipped=[];started=time.monotonic();env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
env['PYTHONDONTWRITEBYTECODE']='1'
(ROOT/'results/ANALYTIC_ANCHORS.json').write_text(json.dumps(anchor,indent=2)+'\n')
budget=[]
old=json.loads((OLD/'results/posterior_attempt_v2/CERTIFICATE_STATUS.json').read_text());oldr=F(old['candidate_one_period_residual_upper']);oldbeta=F(old['beta_lower_rational']);oldtv=F(old['candidate_tail_tv_exact']);oldtail=F(old['orbit_tail_term_upper']);oldtotal=oldr/oldbeta+2*oldtv+oldtail
(ROOT/'results/OLD_BOUND_DECOMPOSITION.json').write_text(json.dumps({'residual_over_beta':str(oldr/oldbeta),'two_candidate_tail':str(2*oldtv),'old_wave_tail':str(oldtail),'total':str(oldtotal),'wave_tail_fraction':str(oldtail/oldtotal),'wave_tail_fraction_float_diagnostic':float(oldtail/oldtotal)},indent=2)+'\n')
oldb=F(ceil_scaled(1/sum((F(2,5)**j/factorial(j) for j in range(25)),F(0))),S)
for K in [4,6,8,10,12,24,29,30,32]:
 tau,parts=wave_tail(K);tv=candidate_tail(data['laws'],K)
 for B in [8,16,24,32,40,48]:
  remain=F(1,1000)-2*tv-2*tau;req=B*data['r']/remain if remain>0 else None
  old_tail=100*oldb**K;oldremain=F(1,1000)-2*tv-old_tail;oldreq=B*data['r']/oldremain if oldremain>0 else None
  budget.append({'K':K,'B':B,'two_wave_tail_upper':str(2*tau),'two_candidate_tail':str(2*tv),'beta_needed_with_old_r':str(req) if req else None,'beta_needed_float_diagnostic':float(req) if req else None,'old_100exp_tail_upper':str(old_tail),'beta_needed_with_old_tail':str(oldreq) if oldreq else None})
(ROOT/'results/ANALYTIC_BUDGET.json').write_text(json.dumps(budget,indent=2)+'\n')
with (ROOT/'logs/ATTEMPTS.jsonl').open('x') as ledger:
 for stage in plan['stages']:
  if best is not None and F(best,S)<=F(1,10000):skipped.append({'stage':stage['id'],'reason':'predeclared_optional_stricter_stop_met'});continue
  remaining=300-(time.monotonic()-started)
  if remaining<=0:skipped.append({'stage':stage['id'],'reason':'predeclared_total_budget_exhausted'});break
  cmd=[sys.executable,'-B',str(ROOT/'src/run_stage.py'),stage['id'],str(best) if best is not None else 'none'];t=time.monotonic()
  with (ROOT/'logs'/f"{stage['id']}.stdout.log").open('x') as out,(ROOT/'logs'/f"{stage['id']}.stderr.log").open('x') as err:
   try:p=subprocess.run(cmd,stdout=out,stderr=err,env=env,timeout=min(remaining,90));code=p.returncode;timed=False
   except subprocess.TimeoutExpired:code=None;timed=True
  record={'stage':stage,'command':cmd,'prior_counts':str(best) if best is not None else None,'exit_code':code,'timed_out':timed,'elapsed_budget_seconds':time.monotonic()-t};ledger.write(json.dumps(record)+'\n');ledger.flush();print(json.dumps(record),flush=True)
  if code!=0:break
  result=json.loads((ROOT/'results'/stage['id']/'CERTIFICATE.json').read_text());completed.append(result)
  if result['posterior_W1_upper_counts'] is not None:
   val=int(result['posterior_W1_upper_counts']);best=val if best is None else min(best,val)
  assert sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file())<=64*1024**2
summary={'status':'acceptance_met' if best is not None and F(best,S)<=F(1,1000) else 'acceptance_not_met','completed_stages':len(completed),'skipped':skipped,'best_W1_upper_rational':str(F(best,S)) if best is not None else None,'best_W1_upper_float_diagnostic':float(F(best,S)) if best is not None else None,'acceptance_threshold':'1/1000','acceptance_met':best is not None and F(best,S)<=F(1,1000),'elapsed_matrix_budget_seconds':time.monotonic()-started,'new_ODE_runs':0,'new_random_experiments':0}
(ROOT/'results/SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
