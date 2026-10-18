"""Independent backward terminal-indicator certificate, own bootstrap and64-bit floor."""
from fractions import Fraction as F
from math import factorial
import json,hashlib,signal,resource,time
import numpy as np
from bounds import ROOT,S,inputs,AnalyticFront,wave_tail,candidate_tail,exp_minus_lower,ceil_scaled
signal.alarm(85);resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));started=time.monotonic()
for p,h in json.loads((ROOT/'protocol/VERIFICATION_FREEZE.json').read_text())['files'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
plan=json.loads((ROOT/'protocol/PLAN_V1.json').read_text());selected=json.loads((ROOT/'protocol/VERIFICATION_PLAN.json').read_text())['stages'];data=inputs();prior=None;results=[];Q=1<<64;ratio=S//Q
def down(x):return x//ratio
def up(x):return (x+ratio-1)//ratio
for stage in plan['stages']:
 if stage['id'] not in selected:continue
 K=stage['K'];B=stage['B'];R=stage['radius'];J=stage['J'];n=2*R+1;ext=list(range(-R-1,R+2));levels=list(range(-R,R+1));offset=-data['levels'][0];front=AnalyticFront(J);ret80=exp_minus_lower(F(2,J))
 retfrac=sum(((-F(2,J))**j/factorial(j) for j in range(10)),F(0));ret64=retfrac.numerator*Q//retfrac.denominator
 fields=[]
 for cell in range(J):
  start=cell*(256//J);end=(cell+1)*(256//J);cv=[];pv=[]
  for side in [0,1]:
   lows=np.array([data['low'][side][start,j+offset] for j in ext],dtype=object);highs=np.array([data['high'][side][end,j+offset] for j in ext],dtype=object);mass=np.array([data['mass_low'][side][start,j+offset] for j in ext],dtype=object)
   al=[];au=[];am=[]
   for idx,j in enumerate(ext):
    l=front.lower(j,side,cell);u=front.upper(j,side,cell+1);a=ret80*max(0,front.lower(j,side,cell)-front.upper(j+1,side,cell))//S
    if prior is not None:l=max(l,max(0,int(lows[idx])-prior));u=min(u,min(S,int(highs[idx])+prior));a=max(a,ret80*max(0,int(mass[idx])-prior)//S)
    assert 0<=l<=u<=S;al.append(l);au.append(u);am.append(a)
   cv.append((lows,highs,mass));pv.append((np.array(al,dtype=object),np.array(au,dtype=object),np.array(am,dtype=object)))
  fields.append((cv,pv))
 V=np.eye(2*n,dtype=object).reshape(2,n,2*n)*Q
 for cycle in range(B-1,-1,-1):
  V[:,1:]=V[:,:-1].copy();V[:,0]=0
  dr=cycle*data['r_counts']
  for cell in range(J-1,-1,-1):
   cv,pv=fields[cell];lo=[np.maximum(0,a[0]-dr) for a in cv];hi=[np.minimum(S,a[1]+dr) for a in cv];mass=[np.maximum(0,a[2]-dr)*ret80//S for a in cv]
   fL=lo[1][1:-1]+pv[1][0][1:-1];cL=mass[1][:-2]+pv[1][2][:-2];fR=lo[0][2:]+pv[0][0][2:];cR=mass[0][1:-1]+pv[0][2][1:-1]
   eL=np.minimum(2*S,hi[1][:-2]+pv[1][1][:-2]);eR=np.minimum(2*S,hi[0][1:-1]+pv[0][1][1:-1]);assert all(fL+cL<=eL) and all(fR+cR<=eR)
   el=np.array([up(int(x)) for x in eL],dtype=object);er=np.array([up(int(x)) for x in eR],dtype=object)
   diag=np.array([(2*Q-el)//J,(2*Q-er)//J],dtype=object);fl=np.array([down(int(x))//J for x in fL],dtype=object);cl=np.array([down(int(x))//J for x in cL],dtype=object);fr=np.array([down(int(x))//J for x in fR],dtype=object);cr=np.array([down(int(x))//J for x in cR],dtype=object)
   def action(T):
    U=diag[:,:,None]*T;U[0,:-1]+=fl[:-1,None]*T[0,1:];U[0]+=cl[:,None]*T[1];U[1,:-1]+=fr[:-1,None]*T[1,1:];U[1,:-1]+=cr[:-1,None]*T[0,1:];return U
   term=V.copy();out=V.copy()
   for degree in range(1,5):term=action(term)//(Q*degree);out+=term
   V=(out*ret64)//Q
 assert min(V.flat)>=0 and all(V.sum(axis=2).reshape(-1)<=Q)
 full=V.reshape(2*n,2*n);core=[side*n+j+R for side in [0,1] for j in range(-K,K+1)];counts=sum(min(full[i,h] for i in core) for h in range(2*n));beta=F(counts,Q);assert beta>0
 tv=candidate_tail(data['laws'],K);tau,_=wave_tail(K);bound=B*data['r']/beta+2*tv+2*tau;current=ceil_scaled(bound);prior=current if prior is None else min(prior,current)
 original=json.loads((ROOT/'results'/stage['id']/'LOWER_MATRIX.json').read_text())['counts'];orig=np.array([[[int(v)/S for v in row] for row in side] for side in original]).reshape(2*n,len(core));check=np.array([[int(full[i,h])/Q for i in core] for h in range(2*n)])
 rec={'stage':stage['id'],'beta_lower_rational':str(beta),'beta_lower_float_diagnostic':float(beta),'W1_upper_rational':str(bound),'W1_upper_float_diagnostic':float(bound),'own_prior_upper_counts_after_stage':str(prior),'max_forward_backward_entry_difference_diagnostic':float(np.max(np.abs(orig-check))),'acceptance_met':bound<=F(1,1000)};results.append(rec)
 (ROOT/'results'/f"BACKWARD_{stage['id']}.json").write_text(json.dumps({'denominator':str(Q),'rows_are_initial_states_columns_terminal':True,'state_levels':levels,'counts':[[str(v) for v in row] for row in full.tolist()]},indent=2)+'\n');print(json.dumps(rec),flush=True)
summary={'status':'passed_independent_backward_certificate','bits':64,'direction':'backward all terminal indicators','stages':results,'best_W1_upper_rational':str(F(prior,S)),'best_W1_upper_float_diagnostic':float(F(prior,S)),'acceptance_threshold':'1/1000','acceptance_met':F(prior,S)<=F(1,1000),'elapsed_budget_seconds':time.monotonic()-started,'new_ODE_runs':0,'random_calls':0}
(ROOT/'results/INDEPENDENT_BACKWARD_CHECK.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
