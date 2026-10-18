# Reproduce the frozen follow-up

The original microscopic dataset and v4.1 are read-only inputs, locked in `protocol/INPUT_LOCK.json`. This directory contains no microscopic model executable or new microscopic seeds. All new randomness is confined to the independently namespaced, prefrozen SDE reference.

The sole reference command, run from the repository root, was:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B research/paper_sde_end_to_end_20261017/src/run_pipeline.py
```

Do not rerun it in this release. Outputs/logs already exist and the `/tmp` binary is also fixed. A fresh independent reproduction uses a clean isolated checkout and destination, preserving the exact seed list, sources, coefficients, decision rules and resource budgets; it must not replace this completed batch. Stage commands, CPU/wall accounting, exit statuses and stdout/stderr are in `logs/`.

The stages were deterministic certified-table export, native compilation, deterministic elementary-function/lookup checks, exact old-micro mean precision precheck, fixed arithmetic cost calibration, the single SDE reference batch, reference confidence intervals, and end-to-end comparisons. The compiler flags, binary hash, lookup hash, seed hash and protocol hash are fixed before the first reference trajectory in `protocol/RUN_INPUT_FREEZE.json`. The earlier code/plan/seeds were frozen in `protocol/SOURCE_FREEZE.json` before any new mean calculation or reference path.

Raw reference terminal values and readouts are exact long-double hex records in `results/reference_endpoints.csv`; every replicate has both lambda rows and a driver hash. The progression file records every256 completed labels. No trajectory or seed is selected from partial data. Raw Gaussian increments are reconstructible from the frozen counter generator and seed, rather than stored as billions of words.

All reference outputs were frozen before the independent read-only audit in `protocol/REFERENCE_DATA_FREEZE.json`. Its verifier source and input-freeze hash are locked in `protocol/INDEPENDENT_CHECK_FREEZE.json`; its full report is `results/INDEPENDENT_INTERVAL_CHECK.json`. This check generates no reference or micro paths. The tool completion is retained separately; there was no retry to obtain its log.

Human-readable findings and every error component are in `EXPERIMENT_REPORT.zh.md` and `ANALYTIC_ERROR_BUDGET.en.md`. Machine-readable exact interval endpoints are in `END_TO_END_MEAN_COMPARISONS.json`, `END_TO_END_D_CDF_COMPARISONS.json` and `REFERENCE_QUALITY.json` under `results/`. The new statistical family is a follow-up family, not retrospective validation of old primary criteria. Archive and readback scripts perform storage verification only and never execute any numerical kernel.
