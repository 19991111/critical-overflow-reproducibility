"""Deterministic coefficient export from frozen certified wave; no ODE or sample."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Context,Decimal,ROUND_CEILING
import hashlib,json,math,numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
cert=json.loads((REPO/'research/paper_profile_certificate_20261017/results/CERTIFIED_READOUTS.json').read_text())
lines=(REPO/'research/paper_v4_release_20261017/results/WAVE_TABLE.hex').read_text().splitlines();assert lines[0]=='257 50 -24 24';table=np.array([[float.fromhex(v) for v in s.split()] for s in lines[1:]])
levels=np.arange(-24,25);phase=np.arange(257)/256
def laws(u):
 u=np.asarray(u);shift=np.floor(u).astype(int);j=(u-shift)*256;idx=np.floor(j).astype(int);t=j-idx;tails=(1-t[:,None])*table[idx]+t[:,None]*table[idx+1];return tails[:,:50],tails[:,50:],shift
def H(u,v):
 x,y,s=laws(u);xx,yy,t=laws(v);a=-np.diff(x,axis=1);ix=np.arange(49)[None,:]+s[:,None]-t[:,None];tails=yy[np.arange(len(u))[:,None],np.clip(ix,0,49)];tails=np.where(ix<0,1,np.where(ix>=50,0,tails));return np.sum(a*tails,axis=1)
weights=np.ones(257);weights[1:-1:2]=4;weights[2:-1:2]=2;weights/=768
G=[]
for j in range(4097):
 d=j/512;vals=H(phase+d/2,phase-d/2)-H(phase-d/2,phase+d/2);vals[-1]=vals[0];G.append(float(np.dot(weights,vals)))
def sigma(u):
 x,y,s=laws(u);a=-np.diff(x,axis=1)-np.diff(y,axis=1);return np.sum(a*np.maximum(levels[None,:]+s[:,None],0),axis=1)/2
dd=np.arange(16385)/512;S=sigma(dd/2)+sigma(-dd/2);S[-1]=16.
gf=list(map(F.from_float,G));sf=list(map(F.from_float,S));assert gf[0]==0 and all(-1<=x<=0 for x in gf)
gslopes=[512*(b-a) for a,b in zip(gf,gf[1:])];sslopes=[512*(b-a) for a,b in zip(sf,sf[1:])];assert all(-2<=x<=0 for x in gslopes) and all(0<=x<=F('0.5001') for x in sslopes) and sf[0]<F('.632')
nodeg=F('.00004861')+F(9472,180*256**4);globalg=nodeg+F(3,512**2);assert globalg<F('.000061')
tailg=max(1+gf[-1],abs(gf[-1]-F(cert['feedback']['8']['lower_rational'])),abs(gf[-1]-F(cert['feedback']['8']['upper_rational'])));assert tailg<F('.000061')
globals=F('.00002440')+F(1,4*512**2);assert globals<F('.000026')
# True readout tail at d32: exact right-front recursion plus directed left exponential sum.
up=Context(prec=70,rounding=ROUND_CEILING)
def dexpupper(r):return up.next_plus(up.exp(up.divide(Decimal(r.numerator),Decimal(r.denominator))))
b=F(dexpupper(F(-18,25)));left=F(11,100)*(b**26+b**25)/(1-b*b)
front={2:F(1),3:F(3,4),4:F(5,8),5:F(3,8)}
for j in range(6,80):
 z=min(front[j-1],front[j-1]*front[j-2]);front[j]=F((z.numerator*2**128+z.denominator-1)//z.denominator,2**128)
right=sum((front[j] for j in range(34,80)),F(0))+front[78]+front[79] # geometric remainder ratio<=1/2
assert front[78]<F(1,2) and (left+right)/2<F(1,10**8)
N=16384;delta_driver=F('.00013');records=[]
for lam,radius,meanerr in [(1,F('.010'),F('.0006')),(2,F('.020'),F('.0010'))]:
 L=2*lam;C0=2*lam*lam+12*lam;e=(F(10*L,17)+F(C0,2))/N+F(C0,2*N*128)
 assert (e+delta_driver)**2<=F(1,200)*(radius-lam*F('.000061')-F('.000001'))**2
 err=F(1,2)*(e+lam*F('.000061'))+F('.000026')+F('.00014')+F('.5001')*(delta_driver+F('.000001'))
 assert err<meanerr
 records.append({'lambda':lam,'Euler_L2_upper':str(e),'driver_L2_upper':str(delta_driver),'mean_total_bound':str(err),'mean_allowance':str(meanerr),'CDF_radius':str(radius),'CDF_coupling_probability_allowance':'1/200'})
tailref=F('.5001')**2/F('1.867')*F(dexpupper(-F('1.867')**2/(2*F('.5001')**2)));assert tailref<F('.00014')
def exp_lower(x):return sum((x**k/math.factorial(k) for k in range(100)),F(0))
assert exp_lower(F(7))>1000 and exp_lower(2*1024*F('.060')**2)>1000 and exp_lower(2*32768*F('.012')**2)>10000
assert F(8,1000)+F(4,10000)<F(1,40)
assert F(3800,2**52)<F(1,10**12) and 128*(F(1,10**6)+F(1,10**9))<delta_driver
# Fixed hex constants checked using rational series bounds.
def hx(s):
 mant,exp=s.lower().split('p');sg=-1 if mant.startswith('-') else 1;mant=mant.lstrip('+-').removeprefix('0x');a,b=mant.split('.');v=F(sg*int(a+b,16),16**len(b));e=int(exp);return v*(2**e if e>=0 else F(1,2**(-e)))
def atan(x,n):return sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(n)),F(0))
a=atan(F(1,5),80);aa=a+F(1,5)**161/161;c=atan(F(1,239),80);cc=c+F(1,239)**161/161;pilo=16*a-4*cc;pihi=16*aa-4*c;pi=hx('0xc.90fdaa22168c235p-2');assert pilo-F(1,2**62)<pi<pihi+F(1,2**62)
lnlo=2*sum((F(1,3)**(2*k+1)/(2*k+1) for k in range(100)),F(0));lnhi=lnlo+2*F(1,3)**201/(201*(1-F(1,9)));ln=hx('0xb.17217f7d1cf79acp-4');assert lnlo-F(1,2**64)<ln<lnhi+F(1,2**64)
with (ROOT/'results/TABLES.hex').open('x') as f:
 f.write('4097 16385 512\n')
 for x in G:f.write(float(x).hex()+'\n')
 for x in S:f.write(float(x).hex()+'\n')
out={'status':'passed','feedback_global_error':str(globalg),'feedback_tail_error':str(tailg),'feedback_allowance':'0.000061','readout_global_error':str(globals),'readout_tail_error_upper':str((left+right)/2),'readout_allowance':'0.000026','feedback_slopes_min_max':[float(min(gslopes)),float(max(gslopes))],'readout_slopes_min_max':[float(min(sslopes)),float(max(sslopes))],'reference_readout_tail_bias_upper':str(tailref),'Euler_and_mean_CDF_bounds':records,'DKW_union_upper':'21/2500','mean_family_error_upper':'2/125','new_ODE_runs':0,'reference_samples':0,'table_sha256':hashlib.sha256((ROOT/'results/TABLES.hex').read_bytes()).hexdigest()}
(ROOT/'results/TABLE_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['readout_tail_error_upper','feedback_global_error','feedback_tail_error','readout_global_error','reference_readout_tail_bias_upper']}))
