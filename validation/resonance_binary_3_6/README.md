# Controlled resonance check: binary impurity spacings (3, 6)

**Date:** 2026-10-16  
**Scope:** A focused finite-size comparison at nominal energies E = -1 and E = +1. This is a validation note, not a claim of complete physical verification or a new discovery.

## Why this check was needed

The earlier archived energy/localization-length comparison showed approximately **96% relative differences** at E = +/-1 when the new simulation used N = 10000. We recovered the **original mean-log-transmission traces** used for the archived fitted values at those two energies. Each original trace contains **5000 values**, indexed n = 0,...,4999. Applying the historical `scipy.stats.linregress` rule to these raw values reproduces the archived xi values to floating-point precision.

The thesis may describe N = 10000 elsewhere; the raw traces establish the fit input length **for these two specific archived points only**. They do not establish the chain length of all historical simulations or why these traces have 5000 values.

## Parameters and methods

- Model: binary impurity-to-impurity spacings `s1=3`, `s2=6`, probability `p=0.5`.
- Impurity potential `U=0.1`; nominal energies `E=-1`, `E=+1`.
- New C++ runs: `N=5000` and `N=10000`, `realizations=1000`, `seed=12345`, `threads=1`.
- Historical fit: all **5000 recorded points** with site indices `0,...,4999`; `xi = -2/slope`.
- New C++ default fit: all **N+1 points** with site indices `0,...,N` (includes the endpoint).
- The new fit on its first 5000 points instead of all 5001 points gives `xi=1745.077524` at either energy, versus `1745.410603` when including the endpoint. The endpoint convention explains a small additional difference.
- The original RNG seed, exact number of realizations and run invocation are **not recovered** from the historical traces. They should not be assumed identical to the new run.

## Measured results

| Nominal E | Historical trace rows | Historical fitted xi | New xi, N=5000 | Difference at N=5000 | New xi, N=10000 |
|---:|---:|---:|---:|---:|---:|
| -1 | 5000 | 1743.289257 | 1745.410603 | 0.1217% | 3411.048012 |
| +1 | 5000 | 1742.922324 | 1745.410603 | 0.1428% | 3411.048012 |

| Nominal E | Historical linear-fit R² | New R², N=5000 | New R², N=10000 |
|---:|---:|---:|---:|
| -1 | 0.821949 | 0.821858 | 0.795563 |
| +1 | 0.821911 | 0.821858 | 0.795563 |

Relative differences use the historical xi as denominator. Full-precision values are in [`resonance_summary.csv`](resonance_summary.csv).

## Interpretation and limits

The apparent ~96% discrepancy at these two energies largely disappears **when the fitting length is matched**. This is strong evidence that the earlier discrepancy was primarily a finite-size/fitting-window comparison issue rather than an obvious implementation error in the new code. The fitted log-transmission curves near these special energies are **not well described by a single straight line** (R² is approximately 0.82 for N=5000 and 0.80 for N=10000). The resulting finite-length `xi=-2/slope` is sensitive to the fitting interval and should not automatically be interpreted as an asymptotic localization length.

This controlled check does **not** prove the general asymptotic transport law, verify all energies/models, or establish that all historical thesis calculations used N=5000. In particular, the full thesis methodology and historical random seeds are not reconstructed here. The numerical archive remains unchanged.

## Files and provenance

- `historical_traces/`: unmodified original raw mean-ln(T) traces supplied by the researcher for nominal E = +/-1.
- `new_runs/`: run configurations, summaries and complete new C++ mean-ln(T) traces for both N values and both energies.
- `comparison_trace_E*.csv`: aligned old and new N=5000 curves.
- `overlay_E*.png`: plots comparing original traces and new N=5000 curves.
- `resonance_summary.csv`: measured xi and fit R² for all four controlled new runs and both original traces.
- `check_fits.py`: independently recomputes slopes, xi and R² from the included raw traces using only Python's standard library.

## Recompute fits from checked-in raw data

From the repository root:

```bash
python validation/resonance_binary_3_6/check_fits.py
```

This checks the historical and new fitted values and confirms the number of data points. It does **not** rerun the C++ Monte Carlo simulations.

## Rerun the C++ simulation

Build the C++ program according to the root README. On Linux/macOS, for example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

Then run the four configurations (output directories are intentionally new):

```bash
./build/quantum_transport --model binary --N 5000 --U 0.1 --s1 3 --s2 6 --p 0.5 --E -1 --realizations 1000 --threads 1 --seed 12345 --out validation/work/resonance_N5000_Em1
./build/quantum_transport --model binary --N 5000 --U 0.1 --s1 3 --s2 6 --p 0.5 --E 1 --realizations 1000 --threads 1 --seed 12345 --out validation/work/resonance_N5000_Ep1
./build/quantum_transport --model binary --N 10000 --U 0.1 --s1 3 --s2 6 --p 0.5 --E -1 --realizations 1000 --threads 1 --seed 12345 --out validation/work/resonance_N10000_Em1
./build/quantum_transport --model binary --N 10000 --U 0.1 --s1 3 --s2 6 --p 0.5 --E 1 --realizations 1000 --threads 1 --seed 12345 --out validation/work/resonance_N10000_Ep1
```

For Windows multi-configuration builds, use the appropriate executable path, commonly `build/Release/quantum_transport.exe`, and `ctest --test-dir build -C Release --output-on-failure`.
