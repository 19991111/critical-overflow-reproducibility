# Frozen end-to-end reference and statistical budget

This is a new follow-up comparison at q=512, kappa=1, lambda=1/2, not an amendment to the old primary checks. The targets are unconditional E[Z/q] versus E[S0(D1)] and the complete terminal D marginal CDF. The reference never substitutes the observed micro D in the readout. All constants below are fixed before new reference trajectories or new micro interval calculations. Internal mathematical derivation, not formalization or external peer review.

## 1. Explicit additive-noise Euler bound without an exponential constant

Write b=lambda*gbar. The proved front formula gives |F'|<=1, |F''|<=3, hence |gbar''|<=24. Also -2<=gbar'<=0, |gbar|<=1. Thus L=2lambda and C0=2lambda^2+12lambda bound |b'| and |b*b'+b''/2|. For h=1/N, hL<=1, the exact Euler map x+h*b(x) is nondecreasing and1-Lipschitz.

The exact local remainder has the Ito representation

R_k = integral_0^h(h-u)b'(X_ku)dB_u + integral_0^h(h-u)(b*b'+b''/2)(X_ku)du.

Its conditional mean has absolute value<=a=C0*h^2/2. Its conditional centered L2 norm is at most s=L*h^(3/2)/sqrt3+C0*h^2/2 (Ito isometry and the variance bound for a variable bounded in absolute value by C0*h^2/2). For the grid error, contraction implies v_(k+1)<=v_k+2a*sqrt(v_k)+a^2+s^2. Induction gives sqrt(v_N)<=sqrt(N)*s+N*a. Consequently

e_Euler <= (L/sqrt3+C0/2)/N + C0/(2*N*sqrtN).

N=16384, sqrtN=128. All scalar verification uses sqrt3>17/10, so it is rational and conservative. There is no fitted strong-rate constant.

Uniform drift table error delta_g contributes at most lambda*delta_g at T=1 by contraction of the TRUE Euler map. Noise-quantization/roundoff is compared through contraction of the APPROXIMATE map; its monotonicity and slope range are explicitly checked before sampling. Deterministic Euler arithmetic error is<=1e-6.

## 2. Global feedback and readout tables from the old certificate

Reuse the frozen257-phase wave table. Its uniform pointwise comparison allowances are4.861e-5 for g and2.440e-5 for S. For gbar nodes on[0,8] spaced1/512, use periodic composite Simpson with256 intervals. The true common-phase fourth derivative is<=9472, adding9472/(180*256^4). Linear interpolation in d adds24*(1/512)^2/8. The total is<0.000061. At d>=8 use the last node; monotonicity and the separately certified gbar(8) interval bound the whole tail, checked against the same0.000061 allowance. Odd reflection covers negative d. There is no extrapolation from a fitted curve.

For S0 use nodes[0,32] at1/512, even reflection and d/2 beyond32. Since sigma'=x0*y0 is between0and1, S0 is actually1/2-Lipschitz; its second derivative is<=2. Interpolation adds(1/512)^2/4, so total<0.000026. At32 set the node to16. The left front bound G(r)<=.11*exp(.72*(r+6)) for r<=-6 gives the two negative-part sums at mean16 starting at r=-32,-31, each with ratio exp(-1.44). The rational right-tail recursion bounds positive parts at mean-16 starting at F(34),F(35). Their half-sum is<1e-8. S0(d)-d/2 is nonnegative and nonincreasing for d>=0, so this tail bound holds for every d>=32. Every table slope is checked: g in[-2,0], S in[0,0.5001]. All guards are frozen; a failure stops the batch.

## 3. Gaussian driver and arithmetic

Each uniform is the midpoint of a2^-52 cell, from a disjoint fixed word. Couple it to a continuous uniform in that cell. Let Delta=2^-52 and R(u)=sqrt(-2log u). In the first cell, twice the two second moments bounds the integrated squared radius difference by4Delta*(1+log2+2log(1/Delta)). On interior cells, Cauchy--Schwarz gives Delta^2/2 times the integral of1/[u^2*(-2log u)]; split at1/2 to bound this by2Delta. The final cell contributes at most8Delta^2/(1-Delta). Angle error is<=pi*Delta and E R(midpoint)^2<=2 by convexity. Thus for either Box--Muller coordinate

E|Z-Z_mid|^2 <=100Delta*(1+log(1/Delta)) <1e-12.

The generator uses no opaque log/sin/cos approximation: log is a24-term atanh series after exact exponent reduction, sin/cos are16-term reduced-quadrant Taylor polynomials, and sqrt uses normalized Newton iteration. For input above1/2, log uses (u-1)/(u+1) directly, avoiding cancellation near1. Ranges: |normal|<9; mantissa series |t|<=1/3; trig argument<=pi/2. The log relative truncation is<=3^-48/[49*(1-1/9)], trig remainders are below1e-25. Binary64-mantissa long-double round-to-nearest and at most a few hundred bounded operations per elementary function give an intentionally loose1e-9 uniform coordinate arithmetic allowance. Pi is a fixed64-bit hex constant checked against rational Machin bounds. Subtractions near1 obey Sterbenz; normalized sqrt and the reduced log keep relative errors bounded. Deterministic edge/grid tests supplement these range/roundoff bounds, rather than certify an opaque libm.

Minkowski over N increments scaled1/sqrtN gives accumulated driver L2<=128*(1e-6+1e-9)<0.00013. This deliberately does not rely on cancellation of driver errors. Stored terminal hex values are exact long-double numbers. All state additions/interpolations, bounded by worst state magnitude9sqrtN+2<1155, fit within the separate1e-6 deterministic arithmetic allowance.

For a CDF coupling failure allowance tau=.005, radii r1=.010 and r2=.020 satisfy exactly

(e_Euler+.00013)^2 <=tau*(r_lambda-lambda*.000061-.000001)^2.

Therefore the exact SDE CDF is between F_computed(x-r)-tau and F_computed(x+r)+tau, globally in x. This horizontal bracket avoids an unproved density bound.

## 4. Strict unconditional means, including tails

Clip all mean observations at B=3, only for inference; retain every unclipped observation. For independent bounded observations V_i, products of1+lambda_bet*(V_i-m) or1+lambda_bet*(m-V_i) have expectation<=1 under the relevant one-sided null. Reject only when capital>=1000, giving one-sided error<=.001 by Markov. Invert on a1e-5 grid outwards. Micro stakes are1 for lower candidates m in[0,1] and2/5 for upper candidates m in[.5,3]; outside these domains the corresponding endpoint is automatically safe. Use exact integer products for micro. Reference stakes are1/8 on[0,3], with a conservative1e-6 log-capital allowance. These are elementary fixed-time betting intervals; the general method has precedent in [Waudby-Smith and Ramdas](https://arxiv.org/abs/2010.09686). No normal/Student approximation is used.

Changing one arrival mark can change the final labelled load L1 by at most2, hence spill by at most1 (equal total mass). Thus Z/q at n=q^2 has McDiarmid tail exp(-2t^2), and its1024-replicate mean upper bound is observed mean+.059 with failure<.001. For mu<=U<B, its clipping bias is bounded by

beta_micro(U) <= (1/4)/(B-U) * exp[-(B-U)^2/(2*(1/4))].

This uses the existing bounded-difference allocation coupling, not a new simulation. Directed65-digit Decimal exponential arithmetic encloses the tail bound; see the [official Decimal documentation](https://docs.python.org/3/library/decimal.html). Bounds are rounded outward.

For the ideal-Gaussian Euler process with approximate drift, contraction and oddness imply E Y1^2<=1. The terminal map of the N Gaussian increments is1-Lipschitz in Euclidean norm. The checked readout has Lipschitz constant0.5001 and value at0<=.632, so W=Shat(Y1) has mean<=1.133 and Gaussian upper tail with variance proxy0.5001^2. Its deterministic clipping bias at3 is<.00014. Coupling back to the computed clipped reference and to the true SDE mean gives the total symmetric reference allowances .0006(lambda1),.0010(lambda2), including profile, interpolation, Euler, driver, arithmetic and clipping. These inequalities are checked before trajectories. The reference mean is unconditional; no observed micro D enters it.

Eight micro betting endpoints, four micro coarse-mean upper bounds and four reference betting endpoints spend at most.016, within the frozen mean-family budget. Reference and micro samples can be dependent across lambda/kernel cells; union bounds do not require cross-cell independence. Individual replicate streams are assumed ideal independent streams, as in the old data; pseudorandom reproducibility does not prove this premise.

## 5. Distribution inference, decisions and cost gate

Four micro CDFs use DKW width.060 at1024 observations; two reference CDFs use.012 at32768. Exact rational Taylor lower bounds on exp verify their union error<.025. Combine these vertical bands with the fixed horizontal radii and tau. Compute the suprema of both directed shifted empirical-CDF differences exactly over every jump, including ties, to obtain an interval for the true Kolmogorov distance. Accept a practical comparison only when its full upper bound<=.20; detect mismatch only when its lower bound>.20. For means the full difference interval must lie in[-.05,.05]; an interval wholly outside detects mismatch. Otherwise the result is inconclusive.

Before any SDE trajectory, require micro mean intervals of width<=.04 and deterministic table/error guards. A fixed262144-pair arithmetic calibration on an integer grid is extrapolated to268435456 Box--Muller pairs, with factor1.5; require<=480 CPU seconds, leaving120 seconds within the600-second reference ceiling. Calibration creates no stochastic path and cannot change R,N,seeds or thresholds. Reference tables/calibration/sampling share the600-second CPU budget; all computation has1200 seconds wall and1GiB memory, one pinned core. If a gate/budget fails, stop and retain evidence, without a replacement batch or seed extension. The final reference mean interval must have width<=.010 and its reference CDF envelope width<=.10 before primary comparisons are promoted.
