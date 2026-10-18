# v4 paper release

The new full manuscript is `paper/paper_core_v4.pdf` (editable root: `paper/paper_core_v4.tex`). The v3 files are old locked inputs and were not modified.

From this release directory, rebuild the PDF without generating data:

```bash
mkdir -p logs
python3 -B src/build_document.py independent_build
```

The attempt name must be unused. This preserves three-pass pdfLaTeX logs. Requirements: Python 3, pdfLaTeX with Latin Modern, AMS, booktabs, graphicx, microtype, enumitem, hyperref and geometry. A direct standalone build also works:

```bash
cd paper
pdflatex -interaction=nonstopmode -halt-on-error paper_core_v4.tex
pdflatex -interaction=nonstopmode -halt-on-error paper_core_v4.tex
pdflatex -interaction=nonstopmode -halt-on-error paper_core_v4.tex
```

The editable source bundle includes the figures, so Matplotlib is not needed to compile. To regenerate the figures from archived observations, run from the release directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B src/build_paper_assets.py
```

This requires NumPy/Matplotlib and the existing result files. It calls no model, ODE or SDE. Initial asset generation printed a Fontconfig cache warning; both vector PDFs were produced, embedded fonts passed document QA, and no figure computation failed.

The actual one-time stochastic command was run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B research/paper_v4_release_20261017/src/run_validation.py
```

It completed once. **Do not run it again in this release.** The runner refuses existing logs and outcomes, and the recorded `/tmp` binary also prevents reuse. All sources, seeds, numerical lookup and binary hashes were frozen before trajectories. Exact stage commands and statuses are in `logs/*.STATUS.json`. A fresh reproduction belongs in a separate clean checkout/output context, preserving the original seed list, decision rules and budgets; it is not an invitation to replace this run or extend its seeds.

To check archived observations without generating trajectories:

```bash
python3 -B src/audit_observations.py
```

The release already contains its output. Archive creation and remote verification do not invoke this audit or any experiment automatically. `EXPERIMENT_REPORT.zh.md`, `CURRENT_EVIDENCE.zh.md`, the frozen plan and machine-readable results state the scope and remaining limits.
