# Measured clean-directory validation

On 2026-10-18 the fixed-evidence core passed **33 checks and 21 native execution stages** in a new source directory and a fresh Python virtual environment. The exact command was:

```sh
/workspace/critical-overflow-venv/bin/python -B /workspace/critical-overflow-clean-2/tools/reproduce.py all --out /workspace/critical-overflow-validation-2
```

The equivalent portable command is `.venv/bin/python -B tools/reproduce.py all --out runs/full`, with a new output directory. The runner locates all inputs relative to itself. It needs no original repository or original absolute path at runtime.

Hardware: AMD EPYC 9V74, Linux x86-64; host exposed five logical CPUs with a four-CPU cgroup quota and 16 GiB memory limit. **The test used CPU 0 only, one thread.** Python 3.12.14, newly installed NumPy 2.3.5, GCC 14.2.0. The tested kernels require Linux APIs and at least 64-bit binary long-double precision.

| Measured scope | Result |
| --- | --- |
| Full command wall time | 221.790054474 seconds |
| Sum of child CPU time | 218.310288 seconds |
| Child maximum resident-set high-water mark | 167,532 KiB, about 163.6 MiB; includes compilation |
| Original finite stochastic batch | 20.374402061 seconds wall; all raw files byte-identical |
| Prospective stochastic batch | 24.588994872 seconds wall; all raw files byte-identical |
| Reused-reference stochastic replay | 102.662560549 seconds wall; endpoint bytes identical; all 128 progress labels identical, measured CPU timings differ as expected |
| Four-stage forward certificate replay | 13.198776382 seconds wall; stage scientific fields and exact matrices match |
| Independent backward certificate replay | 16.863142170 seconds wall; rational summary matches |

All selected frozen files passed length and stored/uncompressed SHA-256 checks. All three exact original C++ sources compiled and their deterministic self-tests passed. The finite source emits inherited misleading-indentation warnings; compilation succeeded and its source bytes were not changed.

The raw-data audit independently reconstructed all **193,536** saved microscopic checkpoint identities and **9,216** outcome metrics. The prospective audit checked all **16,384** raw terminal configurations, the **8,192** progress records and the shared driver hashes for **4,096** labels across configurations. The old SDE mean grid, prospective rational mean/CDF intervals, certified readouts and both regenerated lookup tables matched the sealed values. The reference cost calibration and microscopic precision gates passed before its replay. No failed seeds were replaced.

From all 4,096 λ1 matching endpoints, the independently regenerated mechanism data gave mass-only MAE **0.39633750915527344** and profile-readout MAE **0.017408271929607397**, both using observed imbalance. These are mechanism diagnostics, not unconditional diffusion evidence, statistical confidence, or a published-algorithm comparison. The exact frozen failure-budget rational was separately checked to be strictly below **0.018138367441** under its stated assumptions.

Raw stage commands, exits, CPU/wall accounting and stdout/stderr are under [full-core/](full-core/). The resource field named `child_max_RSS_kib_highwater_so_far` is cumulative child high-water memory, not a separate per-stage peak. The overall wall time includes input checks, decompression, copying and comparisons; summed child CPU omits Python orchestration CPU.

An initial packaging check stopped **before compilation or experiments** because the large reference seed text was compressed under `configs/` but the CLI expected its plain path. The original exact text was exposed without changing a label or checksum, and the corrected core was then tested from a new source/output directory. The first failure report and explanation remain under `attempts/01-config-exposure/`.

## Scope limits

- At this earlier core-only acceptance, the v2 ZIP was unavailable. That block was subsequently resolved by the exact private GitHub transfer; the separate current-paper acceptance below supersedes this source-integration limitation.
- Matrix propagation **was rerun in this core validation**. It was not run during the main thread's earlier attachment-only check. The saved profile trajectory/residual enclosure was reused; no new ODE solve was performed.
- A replay reproduces the fixed data and does not create an independent statistical replication. Shared labels, old reference reuse and the withdrawn old post-hoc joint coverage statement retain their original interpretation.
- Historical dynamic SDE floating-point register snapshots remain absent. The current prospective executable reports its own current environment before word generation; this does not establish the missing historical state.
- These checks do not verify mathematical proofs or establish an effective uniform finite-scale error bound for the complete limiting theory.

## Authoritative SICOMP v2 paper acceptance

After the exact ZIP was transferred in private commit `2a19828175063a37fad1d68c6a8caa70f3af9448`, its 4,054,233 bytes, SHA-256, Git blob and complete 112-entry source manifest were checked locally. The original missing `verification/PACKAGE_CHECKS.json` remains recorded as absent; this section reports newly generated evidence.

The final runner was executed from `/workspace/critical-overflow-paper-clean-5` into a new `/workspace/critical-overflow-paper-validation-5` tree using the same Python virtual environment, extended with pinned plotting dependencies. It requires no original research-workspace path. The clean run passed **80 checks and 16 stages**, with **17.608926217 seconds wall**, **14.529903 seconds child CPU**, and **183,052 KiB** child memory high-water mark (about 178.8 MiB). It generated no random words, trajectories or new numerical integrations. The already validated stochastic/matrix core was not repeated for paper layout work. The current tracked TeX, plotting scripts and plotting inputs were each checked against the authoritative source before use, so local edits cannot be silently ignored.

All six current figures were regenerated from fixed raw outcomes and nodes. The missing endpoint CSV was reconstructed as 25,600 per-endpoint records with the authoritative reader's actual schema. Five plotted CSVs matched every retained numeric value exactly: 4,097 response rows, 16 finite-scale summaries, 196 profile atom rows, 4,096 readout endpoints and 10,004 CDF grid rows. None was substituted from copied plot-data snapshots. All four tables were independently checked against exact source values, including outward bounds, nonzero failure counts and certain/possible distinctions.

Native XeLaTeX plus BibTeX completed **three alternating cross-document rounds**, stabilizing both auxiliary files. Main paper **19 pages**, supplement **12 pages**, matching the input layout counts. All 42 main-document fonts and 26 supplement fonts are embedded. Final logs contain zero undefined references/citations, multiply defined labels, missing characters or overfull boxes. The sources, bibliography and official SIAM class remain byte-identical to the authoritative ZIP. Rendered PDFs differ from the original Tectonic PDF bytes because the build engine/dependency versions differ; no byte-identical PDF claim is made.

All 31 pages were rendered. The main and supplement overview sheets were reviewed; enlarged inspection covered main pages 2, 4, 6, 7, 17 and 18 (all six figures plus the main table) and supplement pages 6, 7 and 8 (all three tables). The final run's 31 rendered page images are pixel-identical to that reviewed build. No clipping or overlapping labels/legends was observed. This is layout inspection, not external review of the mathematics.

The tested toolchain was XeTeX 0.999996 / TeX Live 2025/dev, BibTeX 0.99d and Poppler 25.03.0. The local TeX installation lacked standard `algorithm.sty`; a byte-verified unmodified v0.1 dependency was supplied externally via normal `TEXINPUTS`. Its exact source/blob/hash and standard system-install alternative are documented in `TeX_RUNTIME_DEPENDENCY.json`. This is a compiler dependency, not an original research path or a repository license choice.

Setup failures are preserved: the original compatibility input fails under tested pdfLaTeX despite its no-effect comment; the initial XeLaTeX attempt lacked algorithms; BibTeX initially refused a parent-relative output under `openout_any=p`. The final runner uses XeLaTeX and runs BibTeX inside its output directory with derived bibliography paths. No TeX security policy or paper text was changed. Reports and logs are under `paper/`, and contact sheets under `paper/previews/`.

Remaining limits are mathematical/statistical rather than missing engineering steps: no formal proof verification, no new ODE solve, no effective whole-limit uniform finite-q bound, no proof that fixed PRNG words are iid, and no retroactive historical dynamic SDE register snapshot. The old withdrawn post-hoc joint coverage statement remains withdrawn.
