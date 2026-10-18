"""Sparse, downward integer propagation of every core column; one frozen stage."""
import json,sys,signal,resource,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from bounds import ROOT,S,inputs,wave_tail,candidate_tail,AnalyticFront,exp_minus_lower,ceil_scaled
signal.alarm(85);resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));started=time.monotonic()
plan=json.loads((ROOT/'protocol/PLAN_V1.json').read_text());stage=next(x for x in plan['stages'] if x['id']==sys.argv[1]);prior=None if sys.argv[2]=='none' else int(sys.argv[2]);data=inputs()
K=stage['K'];B=stage['B'];R=stage['radius'];J=stage['J'];assert 256%J==0
lev=list(range(-R,R+1));ext=list(range(-R-1,R+2));n=len(lev);nc=2*(2*K+1);ret=exp_minus_lower(F(2,J));front=AnalyticFront(J);offset=-data['levels'][0]
first=[];orbit=[]
for cell in range(J):
 start=cell*(256//J);end=(cell+1)*(256//J);cl=[];cu=[];cm=[];pl=[];pu=[];pm=[]
 for side in [0,1]:
  l=np.array([data['low'][side][start,j+offset] for j in ext],dtype=object);u=np.array([data['high'][side][end,j+offset] for j in ext],dtype=object);a=np.array([data['mass_low'][side][start,j+offset] for j in ext],dtype=object)
  low=[];up=[];atom=[]
  for j in ext:
   al=front.lower(j,side,cell);au=front.upper(j,side,cell+1);am=max(0,front.lower(j,side,cell)-front.upper(j+1,side,cell));am=(am*ret)//S
   if prior is not None:
    al=max(al,max(0,int(l[j-ext[0]])-prior));au=min(au,min(S,int(u[j-ext[0]])+prior));am=max(am,(ret*max(0,int(a[j-ext[0]])-prior))//S)
   assert 0<=al<=au<=S and 0<=am<=S
   low.append(al);up.append(au);atom.append(am)
  cl.append(l);cu.append(u);cm.append(a);pl.append(np.array(low,dtype=object));pu.append(np.array(up,dtype=object));pm.append(np.array(atom,dtype=object))
 first.append((cl,cu,cm));orbit.append((pl,pu,pm))
M=np.zeros((2,n,nc),dtype=object);col=0
for side in [0,1]:
 for j in range(-K,K+1):M[side,j+R,col]=S;col+=1
trace=[]
def common(m):return sum(min(row) for row in m.reshape(2*n,nc))
for cycle in range(B):
 dr=cycle*data['r_counts']
 for cell in range(J):
  cl,cu,cm=first[cell];pl,pu,pm=orbit[cell]
  vl=[np.maximum(0,a-dr) for a in cl];vu=[np.minimum(S,a+dr) for a in cu];vm=[(np.maximum(0,a-dr)*ret)//S for a in cm]
  fL=vl[1][1:-1]+pl[1][1:-1];cL=vm[1][:-2]+pm[1][:-2];fR=vl[0][2:]+pl[0][2:];cR=vm[0][1:-1]+pm[0][1:-1]
  eL=np.minimum(2*S,vu[1][:-2]+pu[1][:-2]);eR=np.minimum(2*S,vu[0][1:-1]+pu[0][1:-1]);assert all(fL+cL<=eL) and all(fR+cR<=eR)
  diag=np.array([(2*S-eL)//J,(2*S-eR)//J],dtype=object);fl=fL//J;cr=cR//J;fr=fR//J;clow=cL//J
  def apply(T):
   U=diag[:,:,None]*T
   U[0,1:]+=fl[:-1,None]*T[0,:-1];U[1]+=clow[:,None]*T[0];U[1,1:]+=fr[:-1,None]*T[1,:-1];U[0,1:]+=cr[:-1,None]*T[1,:-1]
   return U
  term=M.copy();out=M.copy()
  for degree in range(1,5):term=apply(term)//(S*degree);out+=term
  M=(out*ret)//S
 M[:,:-1]=M[:,1:].copy();M[:,-1]=0
 assert min(M.flat)>=0 and all(M.sum(axis=(0,1))<=S)
 beta_counts=common(M);trace.append({'cycle_completed':cycle+1,'beta_lower_counts':str(beta_counts),'minimum_retained_column_mass_counts':str(min(M.sum(axis=(0,1))))})
 print(json.dumps({'stage':stage['id'],'period':cycle+1,'beta_lower_diagnostic':beta_counts/S}),flush=True)
beta=F(common(M),S);tv=candidate_tail(data['laws'],K);tau,parts=wave_tail(K);bound=B*data['r']/beta+2*tv+2*tau if beta else None;bound_counts=ceil_scaled(bound) if bound is not None else None
outdir=ROOT/'results'/stage['id'];outdir.mkdir()
(outdir/'LOWER_MATRIX.json').write_text(json.dumps({'state_levels':lev,'sides':['X','Y'],'core_levels':list(range(-K,K+1)),'column_order':'X core ascending, then Y core ascending','denominator':str(S),'counts':[[[str(v) for v in row] for row in side] for side in M.tolist()]},indent=2)+'\n')
(outdir/'PERIOD_TRACE.json').write_text(json.dumps(trace,indent=2)+'\n')
result={'status':'certified' if beta else 'no_positive_common_mass','stage':stage,'prior_E_upper_counts':str(prior) if prior is not None else None,'prior_denominator':str(S),'exact_candidate_mean':0,'one_period_r_upper':str(data['r']),'beta_lower_rational':str(beta),'beta_lower_float_diagnostic':float(beta),'candidate_tail_exact':str(tv),'wave_tail_upper_rational':str(tau),'wave_tail_parts':parts,'residual_term_upper_rational':str(B*data['r']/beta) if beta else None,'posterior_W1_upper_rational':str(bound) if bound is not None else None,'posterior_W1_upper_float_diagnostic':float(bound) if bound is not None else None,'posterior_W1_upper_counts':str(bound_counts) if bound_counts is not None else None,'scientific_acceptance_W1':'1/1000','acceptance_met':bound is not None and bound<=F(1,1000),'optional_stop_met':bound is not None and bound<=F(1,10000),'elapsed_budget_seconds':time.monotonic()-started,'new_ODE_runs':0,'random_calls':0}
(outdir/'CERTIFICATE.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
