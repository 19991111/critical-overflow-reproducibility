# SICOMP v2: theorem, figure and table map

The exact authoritative ZIP was read locally after transfer commit 2a19828175063a37fad1d68c6a8caa70f3af9448. Its 4,054,233 bytes have SHA-256 c863979b757ef5375b6a3537d8dfbe7a48608c57c51fdc4c8044f4f5e2350432. All 112 source-manifest entries and the complete 113-file set were verified. This replaces the provisional map. Final compiled labels confirm three main theorems, two corollaries, six figures, one article table and three supplement tables.

The [machine-readable map](../provenance/PAPER_MAP.json) names code, configuration, exact seed protocols, raw inputs and commands. Numerical inputs resolve through provenance/FILES.json to fixed original commits, archives, members and checksums.

## Mathematical results

| Result | Current statement and proof | Numerical relationship and limits |
| --- | --- | --- |
| Theorem 2.1, thm:universal | paper/manuscript/sections/model.tex:53; Section 3 in universality.tex; supplement process_estimates.tex | Fixed matching/product data illustrate two finite-scale kernels. They do not verify the theorem over all admissible kernels. The q^(-1/6) log²q bound is a kernel-comparison bound. |
| Theorem 2.2, thm:profile | model.tex:63; Section 4 in profile.tex; supplement Gini/process estimates | Saved profile, residual enclosure and forward/backward certificate code evaluate the profile. No new ODE solve; numerical tests do not verify existence, uniqueness or attraction proofs. |
| Theorem 2.3, thm:limit | model.tex:87; Section 5 in limit.tex and supplement auxiliary estimates | Prospective endpoints and the reused SDE reference assess finite-size predictions. The complete joint limit has no effective uniform finite-scale error bound here. |
| Corollary 2.4, cor:internalexcess | model.tex:124, proof at line 140 | Figure 1 illustrates the omitted internal contribution. It is not verification of the strict asymptotic conclusion. |
| Corollary 2.5, cor:capacity | model.tex:148, proof at line 164 | Capacity consequences follow from Theorem 2.3, convexity and first-Wasserstein convergence. No new experiment or seed; no algorithmic optimality claim. |

These are mathematical proofs with no simulation seed that can certify them. Model conditions include empty starts, a predetermined doubly stochastic cross kernel, left tie rule, epsilon*q tending to finite lambda, and integer capacity with c_q/m tending to positive kappa. Load-adaptive kernels are excluded.

## Six displayed figures

Run tools/build_paper.py with a new output directory. The exact upstream scripts are retained under tools/paper/ and invoked in a disposable copy after reconstructing their required inputs.

| Figure; label and page | Fixed inputs and script | Actual validation |
| --- | --- | --- |
| 1; fig:profileablation, p.2; readout_teaser | All 4096 prospective lambda1 matching endpoints, q512/kappa1; frozen S lookup; prepare_figures.py then prepare_teaser.py | MAEs 0.39633750915527344 and 0.017408271929607397; lookup allowance at most 0.000026. Both predictions use observed imbalance. No confidence or published-algorithm superiority claim. |
| 2; fig:model, p.4; model_scales | Model definitions; prepare_figures.py | Schematic, no experimental seed. Probabilities/scales and enlarged rendering checked. |
| 3; fig:mechanism, p.6; profile_mechanism | raw_profiles.npz, x_128/y_128, rows 0/128; frozen response nodes; prepare_mechanism.py | NPZ and lookup agree with fixed original inputs. All 196 atom rows and 4097 displayed response rows match. Connecting lines are not additional certificates. |
| 4; fig:kernels, p.7; kernel_universality | Exact I/4 and J/16; prepare_kernel_figure.py | Exact uniform row/column marginals 1/4. Schematic, not an experiment. Marginals alone do not imply independent loads. |
| 5; fig:finite, p.17; finite_scale | All original outcomes at q64/128/256/512, lambda1/2, both kernels; 128/128/256/1024 prefix labels; prepare_figures.py | All 16 plot rows match. Core audited 193536 checkpoints and reconstructed 9216 outcome metrics. RMS uses 16 phases; overflow maximum uses five prescribed offsets. Margins are not confidence bands on plotted means. |
| 6; fig:terminal, p.18; terminal_law | All prospective D/q endpoints and reused 32768-per-lambda SDE hex endpoints, 16384 Euler steps; prepare_figures.py | All 10004 CDF grid rows match. Difference panels use every empirical jump. This is the unconditional diffusion comparison. Core independently audited the exact rational CDF allowances. |

The old internal_contribution illustration remains in the original ZIP and the legacy batched generator may redraw it. It is not a seventh displayed figure and is not a default current-paper asset.

## Four tables

| Table; label and page | Inputs and command | Actual checks |
| --- | --- | --- |
| Article Table 1; tab:prospective, p.18 | Prospective MEAN_COMPARISONS.json, CDF_COMPARISONS.json, raw configurations and reused reference; core recompute plus paper build | All four outward mean intervals and CDF bounds enclose exact fractions. Four configurations share 4096 labels, giving 16384 records, not independent observations. The prospective family remains separate. |
| Supplement SM1; tab:certificate, SM6 | Four data/certificate/K*/CERTIFICATE.json and matrices; core certificates | Common mass displays are lower bounds and W1 bounds round upwards. Both initial failed stages remain. Full forward/backward propagation was rerun in the core, separately from the earlier attachment-only inspection. |
| Supplement SM2; tab:readouts, SM7 | Five selected CERTIFIED_READOUTS.json entries; core recompute plus paper build | All five outward intervals checked. Uniform D_* less than 9.038e-6 and fourth-stage W1 less than 9.034e-6 checked separately. These are not finite-q error bars. |
| Supplement SM3; tab:primary, SM8 | PRIMARY_FAILURE_RATES.json and complete original checkpoints/outcomes; core recompute | All ten criteria, certain/possible counts and outward intervals match. Nonzero failures and ambiguous classifications are retained. |

## Repaired dependencies and preserved provenance

The old workspace endpoint_comparison.csv was absent. The paper runner reconstructs 25600 per-endpoint rows (9216 original and 16384 prospective) from checksummed raw outcomes and the frozen S lookup, with the actual schema consumed by the authoritative script. It does not use copied plot data. The core's earlier same-named CSV is an aggregate report and is not used as this input.

Two small runtime tar containers contain only exact original members required by the legacy reader. They are not the complete historical return archives and do not claim their hashes. Member hashes and source provenance are recorded separately from new container hashes. The documented entrypoint requires no old workspace fallback.

The original revision notes referred to an absent verification/PACKAGE_CHECKS.json. Its absence is recorded in evidence/paper/ZIP_INTAKE.json; evidence/paper/PACKAGE_CHECKS.json is a new report. Old post-hoc corrections, failed certificates and missing historical dynamic SDE register snapshots remain traceable. No old confidence statement is restored.
