"""Exact saved-path proof attempt. No integration, fitting or random calls."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import json,hashlib,signal,resource
import numpy as np
ROOT=Path(__file__).resolve().parents[1];signal.alarm(85);resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
for p,h in json.loads((ROOT/'protocol/POSTERIOR_ATTEMPT_V2_FREEZE.json').read_text())['files'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
theta=F(2,5);exphi=sum((theta**j/factorial(j) for j in range(6)),F(0))+theta**6/F(factorial(6))/(1-theta/7);assert exphi<F(373,250) and 3-2*exphi>F(2,125)
elo=sum((F(1,factorial(j)) for j in range(7)),F(0));assert elo>F(125,46)
norm=json.loads((ROOT/'results/MEAN_ZERO_CANDIDATE.json').read_text());correction=F(norm['exact_W1_sum_change']);assert norm['exact_candidate_mean']=='0'
laws=[{int(k):F(v) for k,v in p.items()} for p in norm['laws']];assert all(sum(p.values())==1 and min(p.values())>=0 for p in laws) and sum(sum(k*v for k,v in p.items()) for p in laws)==0
with np.load(ROOT/'results/tail_h512_K24/raw_profiles.npz') as raw:
 levels=list(map(int,raw['levels']));XX=[[F(float(z)) for z in row] for row in raw['x_128']];YY=[[F(float(z)) for z in row] for row in raw['y_128']]
lo=levels[0];hi=levels[-1];n=len(levels);Q=256;h=F(1,Q);BITS=80;S=1<<BITS
def mass(xs):return [xs[k]-xs[k+1] for k in range(n)]
def tails(ps):
 out=[F(0)]*(n+1)
 for k in range(n-1,-1,-1):out[k]=out[k+1]+ps[k]
 return out
def cap_field(x,y):
 dx=[F(0)]+[2*(x[k-1]-x[k])*y[k-1] for k in range(1,n)]+[F(0)]
 dy=[F(0)]+[2*(y[k-1]-y[k])*x[k] for k in range(1,n)]+[F(0)]
 return dx,dy
def ctrl(p0,p1,d0,d1):
 c1=[p0[k]+h*d0[k]/3 for k in range(n)];back=[p1[k]-h*d1[k]/3 for k in range(n)];neg=-sum(min(F(0),x) for x in back);total=sum(max(F(0),x) for x in back);c2=[max(F(0),x)/total for x in back]
 cs=[p0,c1,c2,p1]
 assert all(sum(p)==1 and min(p)>=0 for p in cs)
 return cs,neg
prodweights=[[(i,k-i,F(comb(3,i)*comb(3,k-i),comb(6,k))) for i in range(max(0,k-3),min(3,k)+1)] for k in range(7)]
elev=[[(i,F(comb(2,i)*comb(4,k-i),comb(6,k))) for i in range(max(0,k-4),min(2,k)+1)] for k in range(7)]
def product3(a,b):return [sum((w*a[i]*b[j] for i,j,w in row),F(0)) for row in prodweights]
def defect(cs,ds):
 total=F(0)
 for side in range(2):
  for k in range(1,n):
   c=cs if side==0 else ds;other=ds if side==0 else cs
   v=[c[j][k] for j in range(4)];d2=[3*(v[j+1]-v[j])/h for j in range(3)]
   derivative=[sum((w*d2[i] for i,w in row),F(0)) for row in elev]
   diff=[c[j][k-1]-c[j][k] for j in range(4)];othercoords=[other[j][k-1 if side==0 else k] for j in range(4)];flow=product3(diff,othercoords)
   total+=sum(abs(derivative[j]-2*flow[j]) for j in range(7))
 boundary=product3([cs[j][n-1] for j in range(4)],[ds[j][n-1] for j in range(4)]);total+=2*sum(boundary)
 return h*total/7
error=F(0);rows=[];maxproj=F(0)
for cell in range(Q):
 x0,y0=XX[cell],YY[cell];x1,y1=XX[cell+1],YY[cell+1]
 for xs in [x0,y0,x1,y1]:assert xs[0]==1 and xs[-1]==0 and min(mass(xs))>=0
 f0,g0=cap_field(x0,y0);f1,g1=cap_field(x1,y1)
 ac,pa=ctrl(mass(x0),mass(x1),mass(f0),mass(f1));bc,pb=ctrl(mass(y0),mass(y1),mass(g0),mass(g1));maxproj=max(maxproj,pa,pb)
 xc=[tails(a) for a in ac];yc=[tails(b) for b in bc];local=defect(xc,yc);local=F((local.numerator*S+local.denominator-1)//local.denominator,S);error+=local;enclosure=error+correction
 def lower(vals):
  z=max(F(0),min(vals)-enclosure);return z.numerator*S//z.denominator
 columns={}
 for name,controls in [('x',xc),('y',yc),('a',ac),('b',bc)]:columns[name]=[lower([row[k-lo] for row in controls]) for k in range(-5,6)]
 rows.append({'cell':cell,'defect_integral_upper':str(local),'cumulative_defect_upper':str(error),'normalized_trajectory_W1_error_upper':str(enclosure),'lower_counts':columns})
def tail_at(xs,k):return F(1) if k<=lo else F(0) if k>hi else xs[k-lo]
return_distance=sum((abs(tail_at(end,k+1)-tail_at(start,k)) for end,start in [(XX[-1],XX[0]),(YY[-1],YY[0])] for k in range(lo-1,hi+2)),F(0));residual=return_distance+error+2*correction
FLO=-4;FHI=4;states=[(s,k) for s in [0,1] for k in range(FLO,FHI+1)];index={v:i for i,v in enumerate(states)};dim=len(states);I=np.eye(dim,dtype=object)*S;M=I.copy()
for row in rows:
 rates=row['lower_counts'];B=np.zeros((dim,dim),dtype=object)
 def rate(key,k):return rates[key][k+5]
 for s,k in states:
  destinations=[((0,k+1),rate('y',k)),((1,k),rate('b',k-1))] if s==0 else [((1,k+1),rate('x',k+1)),((0,k+1),rate('a',k))]
  for dest,value in destinations:
   if dest in index:B[index[dest],index[s,k]]+=value
 assert all(sum(B[:,j])<=S for j in range(dim))
 C=B//Q;term=I.copy();E=I.copy()
 for degree in range(1,5):term=(C@term)//(S*degree);E+=term
 E=(E*127)//128;M=(E@M)//S
 assert min(M.flat)>=0 and all(sum(M[:,j])<=S for j in range(dim))
core=[index[s,k] for s in [0,1] for k in range(-2,3)];common=[min(M[row,j] for j in core) for row in range(dim)];beta=F(sum(common),S)
tv=sum((sum((max(k-2,0)+max(-k-3,0))*v for k,v in law.items()) for law in laws),F(0));tail_upper=F(100)/sum((F(4,5)**j/factorial(j) for j in range(19)),F(0));bound=residual/beta+2*tv+tail_upper if beta>0 else None
out=ROOT/'results/posterior_attempt_v2';out.mkdir()
(out/'TRAJECTORY_LOWER_ENCLOSURE.json').write_text(json.dumps({'levels':list(range(-5,6)),'denominator':str(S),'rows':rows},indent=2)+'\n')
(out/'LOWER_MATRIX.json').write_text(json.dumps({'states':states,'denominator':str(S),'counts':[[str(v) for v in row] for row in M.tolist()],'core_columns':core,'common_row_counts':list(map(str,common))},indent=2)+'\n')
result={'status':'validated_but_precision_target_failed' if beta>0 else 'no_positive_minorization_at_fixed_resolution','candidate_exact_mean':0,'rational_strong_moment_scalar_checks_passed':True,'strong_orbit_cosh_bound':25,'strong_theta':'2/5','exact_probability_Bezier_controls_verified':True,'max_backward_control_negative_mass_repaired':str(maxproj),'raw_period_endpoint_distance':str(return_distance),'integrated_full_ODE_defect_upper':str(error),'exact_normalization_W1_correction':str(correction),'candidate_one_period_residual_upper':str(residual),'candidate_one_period_residual_float_diagnostic':float(residual),'block_B':1,'core_K':2,'chain_radius':4,'beta_lower_rational':str(beta),'beta_lower_float_diagnostic':float(beta),'candidate_tail_tv_exact':str(tv),'orbit_tail_term_upper':str(tail_upper),'posterior_W1_upper_rational':str(bound) if bound is not None else None,'posterior_W1_upper_float_diagnostic':float(bound) if bound is not None else None,'target_W1':'1/1000','precision_target_met':bound is not None and bound<=F(1,1000),'new_ODE_runs':0,'random_calls':0,'true_orbit_digits_certified':False}
assert not result['precision_target_met']
(out/'CERTIFICATE_STATUS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
