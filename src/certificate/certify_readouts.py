from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,Context,ROUND_FLOOR,ROUND_CEILING
from math import factorial,comb
import json,csv,hashlib,signal,resource
import numpy as np
from bounds import ROOT,OLD
signal.alarm(40);resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
for p,h in json.loads((ROOT/'protocol/READOUT_FREEZE.json').read_text())['files'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
plan=json.loads((ROOT/'protocol/READOUT_PLAN.json').read_text());forward=json.loads((ROOT/'results/SUMMARY.json').read_text());backward=json.loads((ROOT/'results/INDEPENDENT_BACKWARD_CHECK.json').read_text());assert forward['acceptance_met'] and backward['acceptance_met']
E=max(F(forward['best_W1_upper_rational']),F(backward['best_W1_upper_rational']));trace=json.loads((OLD/'results/posterior_attempt_v2/TRAJECTORY_LOWER_ENCLOSURE.json').read_text())['rows'];e=max(F(row['normalized_trajectory_W1_error_upper']) for row in trace);D=E+e;N=256
massbound=[1,4,16,128,1280];tailbound=[1,2,8,64,640]
H4=sum(comb(4,i)*massbound[i]*tailbound[4-i] for i in range(5));alpha4=4*sum(factorial(4)//(factorial(i)*factorial(j)*factorial(4-i-j))*massbound[i]*tailbound[j]*massbound[4-i-j] for i in range(5) for j in range(5-i));assert H4==4736 and alpha4==76288
feedback_error=2*D+F(2*H4,180*N**4);alpha_error=8*D+F(alpha4,180*N**4);spill_error=D
with np.load(OLD/'results/tail_h512_K24/raw_profiles.npz') as raw:
 levels=list(map(int,raw['levels']));X=[[F(float(v)) for v in row] for row in raw['x_128'][:N]];Y=[[F(float(v)) for v in row] for row in raw['y_128'][:N]]
AA=[[row[k]-row[k+1] for k in range(len(levels))] for row in X];BB=[[row[k]-row[k+1] for k in range(len(levels))] for row in Y]
assert all(sum(a)==1 and min(a)>=0 for a in AA+BB)
weights=[F(2 if j%2==0 else 4,3*N) for j in range(N)];assert sum(weights)==1
def state(index):
 shift,j=divmod(index,N);return AA[j],BB[j],X[j],Y[j],shift
def H(i,j):
 a,_,_,_,su=state(i);_,_,_,y,sv=state(j);total=F(0)
 for k,p in zip(levels,a):
  loc=k+su-sv-levels[0];tail=F(1) if loc<0 else F(0) if loc>=len(y) else y[loc];total+=p*tail
 return total
def sigma(index):
 a,b,_,_,shift=state(index);return sum((max(k+shift,0)*(p+q) for k,p,q in zip(levels,a,b)),F(0))/2
def decimal(z,rounding):
 c=Context(prec=30,rounding=rounding);return str(c.divide(Decimal(z.numerator),Decimal(z.denominator)))
def interval(value,error,lower=None,upper=None):
 l=value-error;u=value+error
 if lower is not None:l=max(l,F(lower))
 if upper is not None:u=min(u,F(upper))
 return {'estimate_rational':str(value),'error_upper_rational':str(error),'lower_rational':str(l),'upper_rational':str(u),'lower_decimal_outward':decimal(l,ROUND_FLOOR),'upper_decimal_outward':decimal(u,ROUND_CEILING),'estimate_float_diagnostic':float(value)}
alpha=4*sum((weights[j]*sum((a*y*b for a,y,b in zip(AA[j],Y[j],BB[j])),F(0)) for j in range(N)),F(0));feedback={};spill={}
for d0 in plan['feedback_d']:
 d=F(str(d0));off=d*N/2;assert off.denominator==1;off=int(off);val=sum((weights[j]*(H(j+off,j-off)-H(j-off,j+off)) for j in range(N)),F(0));feedback[str(d0)]=interval(val,feedback_error,-1,1)
for h0 in plan['spill_h']:
 for d0 in plan['spill_d']:
  h=F(str(h0));d=F(str(d0));i=(h+d/2)*N;j=(h-d/2)*N;assert i.denominator==j.denominator==1;val=sigma(int(i))+sigma(int(j));base=max(h+d/2,0)+max(h-d/2,0);spill[f'h={h0},d={d0}']=interval(val,spill_error,base)
result={'status':'certified_fluid_wave_readouts','profile_E_upper_rational':str(E),'cycle_e_upper_rational':str(e),'D_E_plus_e_rational':str(D),'profile_W1_acceptance_met':E<=F(1,1000),'quadrature':'exact rational periodic composite Simpson,N=256','fourth_derivative_bounds':{'feedback':2*H4,'origin':alpha4},'error_upper_rationals':{'feedback':str(feedback_error),'origin':str(alpha_error),'spill':str(spill_error)},'origin_alpha':interval(alpha,alpha_error,0,2),'feedback':feedback,'spill_over_q':spill,'finite_q_stochastic_error_certified':False,'new_ODE_runs':0,'random_calls':0}
(ROOT/'results/CERTIFIED_READOUTS.json').write_text(json.dumps(result,indent=2)+'\n')
with (ROOT/'results/CERTIFIED_READOUTS.csv').open('x') as f:
 w=csv.writer(f);w.writerow(['quantity','parameter','estimate_diagnostic','certified_lower','certified_upper'])
 for quantity,items in [('origin_alpha',{'0':result['origin_alpha']}),('feedback',feedback),('spill_over_q',spill)]:
  for key,val in items.items():w.writerow([quantity,key,val['estimate_float_diagnostic'],val['lower_decimal_outward'],val['upper_decimal_outward']])
print(json.dumps({'status':result['status'],'E':float(E),'D':float(D),'origin':result['origin_alpha'],'feedback_d1':feedback['1'],'spill_h0_d0':spill['h=0,d=0']}))
