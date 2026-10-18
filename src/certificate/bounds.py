"""Directed exact bounds. No ODE integration or random calls."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];OLD=REPO/'research/paper_profile_numerics_20261017';S=1<<80
def floor_scaled(x):return x.numerator*S//x.denominator
def ceil_scaled(x):return (x.numerator*S+x.denominator-1)//x.denominator
def exp_minus_lower(x):return floor_scaled(sum(((-x)**j/factorial(j) for j in range(8)),F(0)))
def exp_minus_upper(x):return ceil_scaled(sum(((-x)**j/factorial(j) for j in range(9)),F(0)))
def power_upper(base,n):
 a=S;b=base
 while n:
  if n&1:a=(a*b+S-1)//S
  n//=2
  if n:b=(b*b+S-1)//S
 return a
def check_anchors():
 s=F(89,100);p=s;total=F(0)
 for j in range(6):total+=p;p*=p
 upper=total+p/(1-p);assert upper<3
 exp_lower=sum((F(36,25)**j/factorial(j) for j in range(9)),F(0));assert exp_lower>F(21,5)
 return {'pooled_power_sum_at_0_89_upper_rational':str(upper),'pooled_power_sum_upper_float_diagnostic':float(upper),'exp_1_44_Taylor_lower_rational':str(exp_lower),'left_anchor_G_minus4_upper':'11/100','left_gamma':'18/25'}
UP={-1:S,0:S}
for j in range(1,81):
 k=j//2 if j%2==0 else (j-1)//2
 mark=F(5,4*k) if j%2==0 and k>=1 else F(3,4*k) if k>=1 else F(1)
 UP[j]=min(S,UP[j-1],(UP[j-1]*UP[j-2]+S-1)//S,ceil_scaled(mark))
GM={1:S}
for j in range(0,-81,-1):
 k=j//2 if j%2==0 else (j-1)//2
 mark=F(3,4*(1-k)) if j%2==0 else F(5,4*(1-k))
 GM[j]=min(S,GM[j+1],ceil_scaled(mark))
BNEG=ceil_scaled(1/sum((F(18,25)**j/factorial(j) for j in range(33)),F(0)))
def wave_tail(K):
 assert K>=3
 b=F(BNEG,S);neg=F(11,100)*(b**(2*K-4)+b**(2*K-5))/(1-b*b)
 pos=sum((F(UP[2*j]+UP[2*j+1],S) for j in range(K+1,17)),F(0)) if K<=16 else F(0)
 rest=F(UP[34]+UP[35],S-UP[33])
 if K>16:rest*=F(UP[33],S)**(K-16)
 return pos+rest+neg,{'positive':str(pos+rest),'negative':str(neg)}
class AnalyticFront:
 def __init__(self,J):self.J=J;self.base=exp_minus_upper(F(18,25*J));self.cache={}
 def lower(self,j,side,phase):
  J=self.J;r=(2*j+side)*J-2*phase;ceilr=(r+J-1)//J
  g=GM[min(1,max(-80,ceilr))]
  if r<=-4*J:g=min(g,ceil_scaled(F(11,100)))
  if r<=-6*J:
   n=-r-6*J
   if n not in self.cache:self.cache[n]=power_upper(self.base,n)
   g=min(g,(11*self.cache[n]+99)//100)
  return S-g
 def upper(self,j,side,phase):
  n=((2*j+side)*self.J-2*phase)//self.J
  return S if n<0 else UP[min(n,80)]
def inputs():
 status=json.loads((OLD/'results/posterior_attempt_v2/CERTIFICATE_STATUS.json').read_text());candidate=json.loads((OLD/'results/MEAN_ZERO_CANDIDATE.json').read_text())
 assert candidate['exact_candidate_mean']=='0' and status['status']=='validated_but_precision_target_failed'
 laws=[{int(k):F(v) for k,v in p.items()} for p in candidate['laws']]
 assert all(sum(p.values())==1 and min(p.values())>=0 for p in laws) and sum(sum(k*v for k,v in p.items()) for p in laws)==0
 trace=json.loads((OLD/'results/posterior_attempt_v2/TRAJECTORY_LOWER_ENCLOSURE.json').read_text())['rows'];correction=F(candidate['exact_W1_sum_change']);errors=[correction]+[F(row['normalized_trajectory_W1_error_upper']) for row in trace]
 with np.load(OLD/'results/tail_h512_K24/raw_profiles.npz') as raw:
  levels=list(map(int,raw['levels']));arrays=[raw['x_128'],raw['y_128']]
  low=[];high=[];mlo=[]
  for source in arrays:
   lows=[];highs=[];masses=[]
   for row,e in zip(source,errors):
    exact=[F(float(v)) for v in row];err=ceil_scaled(e);ps=[exact[k]-exact[k+1] for k in range(len(levels))];assert min(ps)>=0 and sum(ps)==1
    lows.append([max(0,floor_scaled(x)-err) for x in exact]);highs.append([min(S,ceil_scaled(x)+err) for x in exact]);masses.append([max(0,floor_scaled(x)-err) for x in ps])
   low.append(np.array(lows,dtype=object));high.append(np.array(highs,dtype=object));mlo.append(np.array(masses,dtype=object))
 return {'r':F(status['candidate_one_period_residual_upper']),'r_counts':ceil_scaled(F(status['candidate_one_period_residual_upper'])),'laws':laws,'levels':levels,'low':low,'high':high,'mass_low':mlo,'phase_nodes':257}
def candidate_tail(laws,K):return sum((sum((max(k-K,0)+max(-k-K-1,0))*v for k,v in p.items()) for p in laws),F(0))
