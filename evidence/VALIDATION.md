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

- The authoritative v2 ZIP remains unavailable in this execution environment. Its source was not locally read, its six figures/four tables were not rendered, and neither TeX document was compiled here. The supplied mapping remains explicitly provisional.
- Matrix propagation **was rerun in this core validation**. It was not run during the main thread's earlier attachment-only check. The saved profile trajectory/residual enclosure was reused; no new ODE solve was performed.
- A replay reproduces the fixed data and does not create an independent statistical replication. Shared labels, old reference reuse and the withdrawn old post-hoc joint coverage statement retain their original interpretation.
- Historical dynamic SDE floating-point register snapshots remain absent. The current prospective executable reports its own current environment before word generation; this does not establish the missing historical state.
- These checks do not verify mathematical proofs or establish an effective uniform finite-scale error bound for the complete limiting theory.
