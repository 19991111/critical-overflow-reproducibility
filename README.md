# Critical Overflow: reproducibility core

Portable CPU reproduction for **Critical Overflow in Weakly Coupled Two-Choice Allocation**. This private repository contains selected fixed experimental evidence, original numerical sources, exact seed lists, and a portable runner. It does not copy the Huawei business workspace or the full research history.

**Current scope:** fixed-evidence core. The authoritative uploaded SICOMP v2 ZIP was resolved in Library but could not be downloaded in this execution environment. Its manuscript source and `prepare_figures.py` have therefore **not been read or integrated here**. The paper-level mapping in [docs/PAPER_MAP.md](docs/PAPER_MAP.md) marks the supplied manuscript descriptions and the remaining integration work. The older v4.2 manuscript is not substituted for v2. No paper link is changed.

## Problem and evidence

Two weakly coupled communities receive two-choice allocations with a fixed left-candidate tie rule. Near total capacity, overflow depends on both the community imbalance and the internal load profiles. The mathematical work identifies the periodic load profile, the imbalance diffusion, and a phase-locked overflow readout; it also compares fixed balanced cross kernels. The retained experiments examine these mechanisms and finite-size behavior. Numerical agreement is not verification of a mathematical proof.

At the critical point used by the terminal experiments, there are `2q` bins, capacity `q/2` per bin, and `q²` arrivals. The prospective batch uses `q=512`, `κ=1`, `λ∈{1,2}`, matching/product kernels, and the same **4,096 labels across four configurations**: **16,384 terminal records**, not 16,384 independent samples. Its 8,192 paired executions include both kernels. The old microscopic batch has 9,216 trajectories and 193,536 saved checkpoints over four q scales. The reused SDE reference has 32,768 paths per λ, 65,536 endpoint rows, and 16,384 Euler steps per path; its λ cells also share their label streams.

## Install and quick check

Requirements: native Linux x86-64, Python 3.11+, `g++` with C++17, and a binary long-double format with at least 64 significand bits. The tested Python version is 3.12.14. One pinned numerical dependency is used. No GPU, CMake, or external research workspace is needed.

```sh
git clone https://github.com/19991111/critical-overflow-reproducibility.git
cd critical-overflow-reproducibility
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B tools/reproduce.py quick --out runs/quick
```

Run commands from any working directory using the script's full path if desired. Every output directory must be new. Frozen inputs are verified by size and SHA-256 before use; existing output directories are refused. Only disposable copies inside the requested output are changed.

## Recompute existing data

This mode generates no new model trajectories or random words. It rebuilds the three native tools, audits microscopic histograms and terminal configurations, recomputes the original failure-rate table and prospective mean/CDF comparisons, recomputes exact readouts, and regenerates both numerical lookup tables. It checks scientific results against the sealed records while retaining measured runtime fields separately.

```sh
.venv/bin/python -B tools/reproduce.py recompute --out runs/recompute
```

Outputs include `PACKAGE_CHECKS.json`, `checks.json`, `stages.json`, stdout/stderr logs, and `analysis/`. `analysis/observed_imbalance_readout.csv` and `analysis/mechanism.json` reconstruct the mechanism comparison from fixed prospective outcomes and the fixed lookup. `analysis/endpoint_comparison.csv` is rebuilt from the original microscopic outcomes using the explicit descriptive schema; compatibility with the unavailable manuscript plotting script is pending. This is more than copying retained plot data, but it is not yet a claim that all six v2 figures have been regenerated.

## Replay the fixed stochastic batches

These commands replay all original labels using unchanged sources, fixed compiler flags, lookup bytes and seeds. They do not add independent statistical evidence or choose replacements based on results. The reference replay retains its original microscopic precision gate and deterministic cost gate. All processes use one pinned CPU and one thread. The available CPU ID is recorded; a different machine cannot be described as the original hardware.

```sh
.venv/bin/python -B tools/reproduce.py rerun --batch finite --out runs/finite
.venv/bin/python -B tools/reproduce.py rerun --batch prospective --out runs/prospective
.venv/bin/python -B tools/reproduce.py rerun --batch reference --out runs/reference
```

The original finite batch uses `-fno-fast-math`; the prospective and reference sources additionally use `-ffp-contract=off`. Sources are preserved byte for byte under `src/`. The runner stores a pre-execution source/binary/seed freeze, compares regenerated raw scientific files by SHA-256, and checks the reference's complete progress labels separately from its machine-dependent CPU timings. Failed runs and partial output remain available. No seeds are generated, extended, screened, or replaced.

## Recompute the certificate matrices

```sh
.venv/bin/python -B tools/reproduce.py certificates --out runs/certificates
```

This replays the frozen forward stage sequence and independent 64-bit backward calculation in separate disposable directories. It compares all four retained forward matrices, stage rational bounds, threshold failures, and the backward rational summary. The archived input freezes are not edited to accept new runtime metadata. The earlier failed profile certificate and both initial failed acceptance stages remain in the fixed evidence. This command reuses the saved profile trajectory and residual enclosure; it does not rerun the ODE solver or establish the mathematical arguments that justify the certificate.

For the complete core reproduction, including data recomputation, all three stochastic replays, and matrix propagation:

```sh
.venv/bin/python -B tools/reproduce.py all --out runs/full
```

Measured results, command exits, runtime, memory high-water marks and remaining failures are recorded in [evidence/VALIDATION.md](evidence/VALIDATION.md). Paper compilation is pending receipt of the authoritative source. Its two TeX documents must eventually follow the supplied `BUILD.txt`, including the required cross-document multiple passes.

## Conditions and limits

- Ideal within-stream and across-label iid input words are mathematical assumptions. The fixed counter-based PRNG is a reproducible implementation, not a proof of iid randomness. Shared drivers across cells are intentional.
- The old post-hoc joint 95% mean/CDF statement was withdrawn. [provenance/history/STATUS_CORRECTED.json](provenance/history/STATUS_CORRECTED.json) and the correction commit remain traceable. Replaying numbers does not restore that statement. The prospective family and its unchanged conditions are kept separately.
- The supplied v2 manuscript description uses `D_* < 9.038e-6` for the uniform profile error and `9.034e-6` for the fourth-stage W1 certificate bound. They are distinct quantities; actual rational bounds are retained in the evidence.
- Only the terminal-law comparison tests an unconditional diffusion prediction. The mechanism MAE comparison conditions both predictions on observed imbalance, carries no confidence claim, and does not compare against published allocation algorithms.
- Fixed-kernel comparison `O(q^(-1/6) log²q)` is not a practical uniform error bound for the whole limit. Load-adaptive kernels are excluded. Qualitative attraction does not supply an effective finite-scale convergence rate.
- The historical SDE executable has no retained dynamic floating-point register snapshot. New runs cannot retroactively supply one. Platform arithmetic contracts and this limitation remain explicit.

## Layout, provenance, contribution

`src/` holds selected original numerical code; `configs/` holds original plans and labels; `data/` holds selected raw observations, configurations, lookup inputs and sealed comparisons, compressing large files losslessly; `tools/` contains the portable orchestration. `provenance/FILES.json` maps every selected frozen file to its fixed commit, archive/member, uncompressed checksum and stored checksum. `provenance/ARCHIVES.json` records the original archive checksums. `docs/` holds the mapping and selected proof material. Historical documents are retained as sources, not current paper entrypoints.

For a contribution, describe the change, use a new output directory, retain the fixed seed lists, and attach the relevant `PACKAGE_CHECKS.json` and logs. An additional experiment needs its own source/configuration freeze and clearly stated scope. Do not overwrite sealed data or infer proof correctness from code tests.

License: **pending author decision**. No third-party open-source license has been selected or granted by this repository setup.
