from pathlib import Path
from fractions import Fraction as F
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
cert=json.loads((REPO/'research/paper_profile_certificate_20261017/results/CERTIFIED_READOUTS.json').read_text());D=F(cert['D_E_plus_e_rational']);interp=D+F(1,256**2)
assert 2*interp+F(1,10**9)<F('0.00004861') and interp+F(1,10**7)<F('0.00002440')
with np.load(REPO/'research/paper_profile_numerics_20261017/results/tail_h512_K24/raw_profiles.npz') as raw:
 levels=raw['levels'];X=raw['x_128'];Y=raw['y_128'];assert X.shape==Y.shape==(257,50)
 with (ROOT/'results/WAVE_TABLE.hex').open('x') as f:
  f.write('257 50 -24 24\n')
  for j in range(257):f.write(' '.join(float(v).hex() for v in np.r_[X[j],Y[j]])+'\n')
out={'profile_plus_cycle_error':str(D),'linear_interpolation_W1_addition':'1/65536','uniform_feedback_allowance':'0.00004861','uniform_spill_allowance':'0.00002440','table_sha256':hashlib.sha256((ROOT/'results/WAVE_TABLE.hex').read_bytes()).hexdigest(),'new_ODE_runs':0}
(ROOT/'results/WAVE_LOOKUP_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
