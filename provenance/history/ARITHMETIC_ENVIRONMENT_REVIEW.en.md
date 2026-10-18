# Evidence for the arithmetic contract, and its limits

The received v2 analytic supplement supports the frozen coordinate1e-9 and joint state/readout1e-6 allowances under binary round-to-nearest, p>=64, finite normal operation ranges, exact radix-two frexp/ldexp, and the declared input-word model. It gives coordinate error<1e-12 and state rounding plus readout rounding/.5001<4e-12. We inspected its bounds against the frozen implementation:24 same-sign log-series terms; seven normalized Newton steps;16 Horner multiply-subtract stages with17 retained trig terms; exact dyadic h and sqrt(h); stored binary64 table values; and contractive approximate Euler maps. No trajectory, seed or statistical endpoint is recomputed.

The latest constant correction matters: the exporter checked membership around rational SERIES INTERVALS, not around exact pi or log2 directly. Including their widths gives |PI-pi|<4.01u and |LN2-ln2|<1.01u, u=2^-64. The reduced-angle claim is absolute error<5u. The later loose budgets remain unchanged. The received v1 is retained as superseded provenance; v2 is authoritative.

## Saved evidence from the actual original build and run

- The actual saved native command used g++14.2.0, -O3, -fno-fast-math and -ffp-contract=off. Its exit status was0. Those flags disable the transformations excluded by the supplement.
- The original deterministic selftest completed before sampling and recorded long_double_digits=64; its code checked radix2 and the type's round_style. The sampler uses its own log/sqrt/sincos implementations; comparisons with libm occur only in deterministic selftests.
- The original environment was Linux x86_64 on the existing AMD EPYC9V74 with glibc2.41. The frozen binary, source, table and inputs retain their old hashes. The source contains no fesetround call.

The type's round_style is not a dynamic fegetround reading. The old batch did not save x87/MXCSR control registers or per-operation underflow traces. Consequently the original runtime control state is not retroactively established by a new probe. The analytic conclusion remains explicitly conditional on the actual original arithmetic obeying the contract. This is a limitation of the recorded evidence, not an assertion that the old state was wrong.

## Current non-sampling environment probe

A separate tiny program, not the sampler, was compiled with the same optimization/arithmetic flags. It generated zero random words and invoked no microscopic or SDE trajectory. Its raw stdout/stderr and command are retained, with source hash, in logs/ and results/ENVIRONMENT_CONTRACT_CHECK.json.

It observed radix2, p=64, sizeof(long double)=16, IEC559 type support, exponent limits -16381 and16384, FE_TONEAREST, x87 control word895(0x037f) with64-bit precision and nearest rounding, and MXCSR8064(0x1f80) with nearest rounding and both FTZ and DAZ disabled. Representative frexp/ldexp roundtrips and a nearest-even half-ulp example passed. This supports the platform settings and wide normal exponent range used by the supplement. These representative checks are evidence, not a universal floating-point formalization or a historical register snapshot.

The supplement's range argument, not those examples alone, is what justifies its no-overflow/no-underflow premise for the relevant normal computations. Stored binary64 hexadecimal coefficients are exactly representable in the wider-significand long-double type. Exact power-of-two exponent operations are required to stay in the representable range. Zero produced by exact cancellation is treated with an absolute-error bound.

## Distinct statistical premise

The arithmetic proof does not establish random-word independence. The stochastic contract requires iid uniform words across replicate labels AND within every underlying stream, including each52-bit Box--Muller input field; the deterministic SplitMix implementation is a reproducible substitute. Neither disjoint seeds nor this environment probe proves physical or mathematical iid behavior. Even under that ideal word contract, post-hoc micro mean recipe selection still needs a separate coverage argument. No overall eight-comparison joint95% claim is restored by this supplement.
