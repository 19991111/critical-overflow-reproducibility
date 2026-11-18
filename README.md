# Critical Overflow in Two-Choice Allocation: Reproducibility Repository

Public reproducibility code, certified numerical enclosures, fixed protocols, and simulation datasets for **Critical Overflow in Two-Choice Allocation: Coupling and Hidden Capacity Loss**.

## Versioned reproducibility bundle

Download release [`reproducibility-v2026.11`](https://github.com/19991111/critical-overflow-reproducibility/releases/tag/reproducibility-v2026.11) and its archive asset `critical_overflow_reproducibility_bundle.zip` (SHA-256 `549448d6410f399033da66011af7ff73d69151c0fd24943615a2c63edb44e78e`, 597 members). The bundle supplies:

- Certified rational profile and readout enclosures, arithmetic supplement, probability contract, and audit receipts (`evidence/` and `CERTIFICATE_REPRODUCE.txt`);
- Unconditional mean, coupling trade-off, linear-drift baseline, and CDF validation batches across fill horizons `kappa in {0.5, 1, 2}` and couplings `lambda in {1, 2}` at `q = 512`;
- Multi-scale boundary-layer experiments across `q in {512, 1024, 2048, 4096}`, support-aware lower bounds, frozen three-size probability forecasts, and the held-out `q = 4096` evaluation;
- The 24-configuration internal candidate-graph and tie-breaking study across `q in {128, 256, 512}` and its 60 simultaneous intervals;
- Read-only vector figure rendering scripts (`prepare_tcs_figures.py` and preserved baseline scripts) and file-level SHA-256 manifests (`MANIFEST.json`), with native Linux reproduction instructions in each experiment's `REPRODUCE.txt`.

The repository root provides the core C++17 / Python reproduction runner (`tools/reproduce.py`) for the original finite-scale, terminal-law, SDE reference, and rational matrix-certificate computations, with every selected frozen input verified by length and SHA-256 via [`provenance/FILES.json`](provenance/FILES.json).

## Problem and evidence

Two weakly coupled communities of `q = 2m` finite-capacity bins receive greedy two-choice allocations without departures under equal-side complete-bipartite internal sampling (`K_{m,m}`) and left-favored ties, while cross-community comparisons occur at rate `varepsilon_q` with `varepsilon_q q -> lambda`. At nominal fill `n_* = 2 q c_q ~ kappa q^2`, pooling each community into a single macro-resource of capacity `q c_q` (`|Delta|/2`) omits a positive order-`q` internal overflow contribution caused by intra-community coexistence of overloaded and underfilled bins. Evaluating the attracting periodic two-type lattice profile at the same normalized community demand imbalance predicts normalized overflow with vanishing mean absolute error, couples to a scalar nonlinear restoring diffusion for community imbalance (where lattice phase averages out of drift but is retained in the local overflow curve), extends unchanged across arbitrary predetermined doubly stochastic cross kernels, and contrasts with an order-`q^(3/2)` overflow regime under internal matchings. The retained experiments test these asymptotic predictions and finite-size boundary effects; numerical agreement illustrates finite-size behavior and is not a substitute for the mathematical proofs.

In the `kappa = 1` terminal-law batch (`2q` bins, capacity `c = q/2`, `n_* = q^2` arrivals at `q = 512`, `lambda in {1, 2}`, matching and product cross kernels), the four configurations share **4,096 replicate seed labels** (**16,384 terminal records** across 8,192 paired matching/product executions). The finite-scale checkpoint batch contains 9,216 trajectories and 193,536 saved checkpoints across `q in {64, 128, 256, 512}`. The reused SDE reference contains 32,768 paths per `lambda` (65,536 endpoint rows, 16,384 Euler steps per path) with shared label streams across `lambda`.

## Internal candidate-graph and tie-breaking study

The 24-cell topology and tie-rule study compares four uniform-marginal internal bipartite graphs (`K_{m,m}`, union of 16 random permutations, degree-two cycle, and disconnected matching) under left-favored and fair ties across `q in {128, 256, 512}` at `kappa = lambda = 1` (1,024 trajectories per cell; 24,576 trajectories in total). At `q = 512` with left-favored ties, sample mean normalized omitted internal overflow on the degree-two cycle and disconnected matching reaches 5.83 and 21.38 times the `K_{m,m}` baseline, respectively, aligning with the proved `Theta(q^(3/2))` matching law under weak cross coupling.

To inspect or replay this study standalone, download [`mechanism_robustness_v1.zip`](paper/fcs/mechanism_robustness_v1.zip) (also included inside `critical_overflow_reproducibility_bundle.zip`), extract it into a fresh directory, and run `python verify_return.py` to check saved file hashes, integer endpoint identities, cell counts, means, and Monte Carlo standard errors; see its `REPRODUCE.txt` for fixed-seed native Linux replay.

## Install and quick check

Requirements: native Linux x86-64, Python 3.11+ (tested on Python 3.12.14), `g++` with C++17, and binary `long double` with at least 64 significand bits. One pinned numerical dependency (`numpy`) is used; no GPU or CMake is required.

```sh
git clone https://github.com/19991111/critical-overflow-reproducibility.git
cd critical-overflow-reproducibility
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B tools/reproduce.py quick --out runs/quick
```

Every `--out` target must be a new directory; frozen inputs are verified by byte size and SHA-256 before execution, and only disposable copies inside the requested output directory are written.

## Recompute existing data

This mode generates no new trajectories or random words. It compiles the native verification tools, audits microscopic histograms and terminal configurations, recomputes the failure-rate table and terminal-law mean/CDF intervals, verifies exact rational readouts, and regenerates both numerical lookup tables against sealed records:

```sh
.venv/bin/python -B tools/reproduce.py recompute --out runs/recompute
```

Outputs include `PACKAGE_CHECKS.json`, `checks.json`, `stages.json`, logs, and `analysis/` (including `analysis/observed_imbalance_readout.csv` and `analysis/mechanism.json`).

## Replay the fixed stochastic batches

These commands replay all original seed labels on native Linux using unchanged C++ sources, fixed compiler flags (`-fno-fast-math`, plus `-ffp-contract=off` for the terminal-law and reference executables), frozen lookup bytes, and pinned single-thread execution:

```sh
.venv/bin/python -B tools/reproduce.py rerun --batch finite --out runs/finite
.venv/bin/python -B tools/reproduce.py rerun --batch prospective --out runs/prospective
.venv/bin/python -B tools/reproduce.py rerun --batch reference --out runs/reference
```

## Recompute the rational certificate matrices

```sh
.venv/bin/python -B tools/reproduce.py certificates --out runs/certificates
```

This replays the four forward directed-rounding matrix stages and the independent 64-bit backward propagation from the saved candidate trajectory and residual enclosure, checking exact rational bounds (`E < 9.034e-6`, `D_* = E + e < 9.038e-6`, and `s_* = S_0(0) in [0.63127654, 0.63129462]`). To run the complete core suite in one command:

```sh
.venv/bin/python -B tools/reproduce.py all --out runs/full
```

Execution records and resource measurements are documented in [`evidence/VALIDATION.md`](evidence/VALIDATION.md).

## Rebuild the preserved baseline document intake

To rebuild the preserved baseline manuscript intake under `paper/` on Debian/Ubuntu Linux:

```sh
sudo apt-get install --no-install-recommends texlive-xetex texlive-latex-extra texlive-science poppler-utils
.venv/bin/python -m pip install -r requirements-paper.txt
.venv/bin/python -B tools/build_paper.py --out runs/paper
```

## Methodological scope and assumptions

- Within-stream and across-label iid uniform random words are mathematical model assumptions; the counter-based SplitMix64 generator provides deterministic stream reproducibility, and shared seed labels across configurations are accounted for by union bounds.
- Uniform profile error `D_* < 9.038e-6` and fourth-stage `W_1` certificate bound `E < 9.034e-6` are distinct certified quantities verified against exact rational values in `data/certificate/CERTIFIED_READOUTS.json`.
- Unconditional diffusion validation uses independent microscopic and SDE streams, whereas the conditioned readout MAE diagnostic evaluates both `|Delta|/2` and `S_0(Delta/q)` at realized microscopic imbalances without a simultaneous confidence claim.
- Cross-kernel universality (`O(q^(-1/6) log^2 q)` coupling bound over `Theta(q^(5/3))`-step windows) applies to predetermined doubly stochastic cross kernels under `K_{m,m}` internal sampling; load-adaptive cross routing is excluded.
- Historical arithmetic and statistical correction records (including [`provenance/history/ARITHMETIC_SUPPLEMENT.v2.en.txt`](provenance/history/ARITHMETIC_SUPPLEMENT.v2.en.txt) and [`provenance/history/STATUS_CORRECTED.json`](provenance/history/STATUS_CORRECTED.json)) remain preserved verbatim for provenance auditability.

## Repository layout, historical provenance, and access

- `src/`, `configs/`, `data/`, `tools/`: frozen C++17 / Python sources, protocols, seed lists, compressed raw observations, lookup tables, and reproduction runners.
- `provenance/`: file-level SHA-256 map (`FILES.json`), original archive digests (`ARCHIVES.json`), baseline paper map (`PAPER_MAP.json`), and historical audit notes (`provenance/history/`).
- Historical snapshots (`source_inputs/`, `paper/`, `paper/fcs/`, and earlier release tag `peva-data-20261023`) are retained unchanged at their original paths and SHA-256 checksums so that all historical commits and provenance references remain verifiable; see [`docs/PAPER_MAP.md`](docs/PAPER_MAP.md) for the preserved baseline intake map.

Copyright and access: All rights reserved by the authors. Public read, download, and execution access is provided for peer review, scientific verification, and academic reproduction.
