# Critical Overflow: reproducibility repository

Reproducibility materials for **Critical Overflow in Weakly Coupled Two-Choice Allocation**.

The current version is prepared for **Frontiers of Computer Science (FCS)**. Formal submission remains pending. Read the [main article](paper/fcs/paper_fcs.pdf), [supplement](paper/fcs/supplement_fcs.pdf), and [three-page highlights](paper/fcs/highlights_fcs.pptx). The [LaTeX source package](paper/fcs/critical_overflow_fcs_source.zip) includes build instructions. The [supporting records](paper/fcs/supporting_records_fcs.zip) contain the unchanged arithmetic contracts, protocols, and audit records referenced in the supplement. The [saved sensitivity study](paper/fcs/mechanism_robustness_v1.zip) supplies its fixed source, seeds, raw endpoints, and verification records. [File hashes](paper/fcs/FILES.json) identify this delivery.

The original numerical source code, fixed seeds, raw data, and Linux reproduction commands remain unchanged. Files under `paper/` outside `paper/fcs/` retain the validated SICOMP-format baseline used by the reproduction scripts. That baseline was not submitted to SICOMP. The [baseline theorem/figure/table map](docs/PAPER_MAP.md) refers to its original page numbers and source locations. The current FCS delivery includes revised proofs and presentation, a fixed sensitivity study completed on native Linux, and the finite-lambda matching comparison in Supplement SM9.2.

The original uploaded ZIP remains unchanged in `source_inputs/`, with SHA-256 `c863979b757ef5375b6a3537d8dfbe7a48608c57c51fdc4c8044f4f5e2350432`. All 112 source-manifest entries and the complete 113-file package were verified during baseline preparation. The following reproduction commands rebuild that fixed baseline; use the current FCS source package to compile the journal-format documents.

## Problem and evidence

Two weakly coupled communities receive two-choice allocations with a fixed left-candidate tie rule. Near total capacity, overflow depends on both the community imbalance and the internal load profiles. The mathematical work identifies the periodic load profile, the imbalance diffusion, and a phase-locked overflow readout; it also compares fixed balanced cross kernels. The retained experiments examine these mechanisms and finite-size behavior. Numerical agreement is not verification of a mathematical proof.

At the critical point used by the terminal experiments, there are `2q` bins, capacity `q/2` per bin, and `q²` arrivals. The prospective batch uses `q=512`, `κ=1`, `λ∈{1,2}`, matching/product kernels, and the same **4,096 labels across four configurations**: **16,384 terminal records**, not 16,384 independent samples. Its 8,192 paired executions include both kernels. The old microscopic batch has 9,216 trajectories and 193,536 saved checkpoints over four q scales. The reused SDE reference has 32,768 paths per λ, 65,536 endpoint rows, and 16,384 Euler steps per path; its λ cells also share their label streams.

## Current FCS sensitivity study

The saved study compares four internal candidate graphs under left-favored and fair ties at `q=128,256,512` and `kappa=lambda=1`.
Its 24 configurations contain 24,576 trajectories in total.
At `q=512` with left-favored ties, the cycle and disconnected-matching sample means of omitted overflow are 5.8 and 21.4 times the complete-bipartite mean.
These ratios compare unconditional sample means of the aggregation gap.
The main critical law uses complete-bipartite internal sampling; Supplement SM9.2 gives the separate order-`q^(3/2)` matching result under weak cross sampling.

Download [mechanism_robustness_v1.zip](paper/fcs/mechanism_robustness_v1.zip) and extract it into a new directory.
Run `python verify_return.py` there to verify the saved file hashes, integer endpoint identities, cell membership, means, and Monte Carlo standard errors.
The archive's `REPRODUCE.txt` describes fixed-seed replay on native Linux.
The LaTeX source package also contains this experiment under `reproducibility/mechanism_robustness_v1/`.
The simultaneous intervals for 24 means and 36 graph contrasts are supplied in `verification/sensitivity_intervals_v1/` inside that source package.
The original numerical records and historical failed runs remain at their existing repository paths and commits.

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

Outputs include `PACKAGE_CHECKS.json`, `checks.json`, `stages.json`, stdout/stderr logs, and `analysis/`. `analysis/observed_imbalance_readout.csv` and `analysis/mechanism.json` reconstruct the mechanism comparison from fixed prospective outcomes and the fixed lookup. The core's `analysis/endpoint_comparison.csv` contains aggregated descriptive fields. The paper runner below separately reconstructs the actual per-endpoint schema needed by the authoritative plotting script; it does not use those aggregates as a substitute.

## Rebuild the validated baseline figures and documents

In addition to the core requirements, install XeLaTeX, BibTeX, Poppler tools, and standard LaTeX packages used by the bundled SIAM class, including **algorithms** (`algorithm.sty`). Typical Debian packages are `texlive-xetex`, `texlive-latex-extra`, `texlive-science` and `poppler-utils`. No shell escape is enabled. The original pdfLaTeX compatibility claim fails on the tested engine, so this entry uses the same XeTeX engine family as the original Tectonic build.

On a Debian/Ubuntu system where these optional paper-build dependencies are missing:

```sh
sudo apt-get install --no-install-recommends texlive-xetex texlive-latex-extra texlive-science poppler-utils
```

```sh
.venv/bin/python -m pip install -r requirements-paper.txt
.venv/bin/python -B tools/build_paper.py --out runs/paper
```

This starts with an empty output tree, checks the authoritative ZIP and editable sources, reconstructs 25,600 original/prospective endpoint rows from fixed raw data, and regenerates every displayed figure. It checks all five plotted datasets against the authoritative retained datasets and all four tables against rational results, including outward rounding and nonzero failures. The two PDFs are built in alternating multiple passes until cross-document auxiliaries stabilize, then checked for unresolved references/citations, duplicate labels, missing characters, overfull boxes, embedded fonts and successful rendering of every page.

`runs/paper/paper/output/pdf/` holds the PDFs; `runs/paper/render/` holds every rendered page and contact sheets. `DOCUMENT_CHECKS.json`, `PLOT_DATA_CHECKS.json`, `checks.json`, `stages.json` and logs are the new acceptance evidence. The uploaded revision notes referred to a missing `verification/PACKAGE_CHECKS.json`; that absence is recorded, and the new top-level acceptance report is not represented as an original file.

The legacy plotting reader receives small runtime tar containers containing only byte-identical selected members. These are explicitly **not** the complete historical return archives and do not claim their archive hashes. Their members map back through `provenance/FILES.json`; original complete archive hashes remain in `provenance/ARCHIVES.json`. This keeps the repository small while rebuilding from raw outcomes, not copied plot data.

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

For the entire release, run the `all` command above followed by `build_paper.py` in two new output directories. Measured results, command exits, runtime, memory high-water marks and retained failed setup attempts are recorded in [evidence/VALIDATION.md](evidence/VALIDATION.md). The paper build follows the supplied cross-document pass dependency with native XeLaTeX and BibTeX; it retains the original class, source and bibliography.

## Conditions and limits

- Ideal within-stream and across-label iid input words are mathematical assumptions. The fixed counter-based PRNG is a reproducible implementation, not a proof of iid randomness. Shared drivers across cells are intentional.
- The old post-hoc joint 95% mean/CDF statement was withdrawn. [provenance/history/STATUS_CORRECTED.json](provenance/history/STATUS_CORRECTED.json) and the correction commit remain traceable. Replaying numbers does not restore that statement. The prospective family and its unchanged conditions are kept separately.
- The v2 manuscript uses `D_* < 9.038e-6` for the uniform profile error and `9.034e-6` for the fourth-stage W1 certificate bound. They are distinct quantities, checked separately against the retained rational values.
- Only the terminal-law comparison tests an unconditional diffusion prediction. The mechanism MAE comparison conditions both predictions on observed imbalance, carries no confidence claim, and does not compare against published allocation algorithms.
- Fixed-kernel comparison `O(q^(-1/6) log²q)` is not a practical uniform error bound for the whole limit. Load-adaptive kernels are excluded. Qualitative attraction does not supply an effective finite-scale convergence rate.
- The historical SDE executable has no retained dynamic floating-point register snapshot. New runs cannot retroactively supply one. Platform arithmetic contracts and this limitation remain explicit.

## Layout, provenance, contribution

`paper/fcs/` holds the current FCS delivery; the other files under `paper/` retain the validated baseline; `src/` holds selected original numerical code; `configs/` holds original plans and labels; `data/` holds selected raw observations, configurations, lookup inputs and sealed comparisons, compressing large files losslessly; `tools/` contains orchestration and exact upstream figure scripts. `provenance/FILES.json` maps every selected frozen numerical file to its fixed commit, archive/member, uncompressed checksum and stored checksum. `provenance/ARCHIVES.json` records original archive checksums; `provenance/manuscript/` records the authoritative paper intake. Historical documents are retained as sources, not current paper entrypoints.

For a contribution, describe the change, use a new output directory, retain the fixed seed lists, and attach the relevant `PACKAGE_CHECKS.json` and logs. An additional experiment needs its own source/configuration freeze and clearly stated scope. Do not overwrite sealed data or infer proof correctness from code tests.

License: **pending author decision**. No third-party open-source license has been selected or granted by this repository setup.
