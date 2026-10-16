# Controlled resonance check: binary impurity spacings (2, 5)

Date: 2026-10-16

**Scope:** A targeted numerical and transfer-matrix analysis of the exact resonances at `E = -1` and `E = +1`, for a one-dimensional nearest-neighbor tight-binding chain with `U = 0.1` and binary impurity-to-impurity spacings `s = 2` or `5` (`p = 0.5`). This is not a validation of every model or energy.

## Why this check was needed

The archived MSc localization-length values at nominal `E = +/-1` are on the order of `10^9`, and the new `N=10000` simulation gives different values of the same order. The archived raw mean-log-transmission files were recovered: **both have 5000 rows**, and the original post-processing convention `xi = -2 / slope` reproduces the archived fitted values. The original C++ code's energy grid is offset from exact integer energies by `eps=1e-8`, so the archived file names indicate **nominal** rather than necessarily mathematically exact energies.

The key result is that **neither large fitted xi is a meaningful finite estimate of an asymptotic localization length at the exact resonance**. The log-transmission traces fluctuate around a small, bounded negative value rather than decreasing linearly with the number of sites.

## Numerical comparison

All new runs use the C++ code in this repository, with `realizations=1000`, `seed=12345`, `threads=1`, `U=0.1`, `s1=2`, `s2=5`, `p=0.5`, and exact `E=-1` or `E=+1`. The C++ trace includes site indices `0,...,N` (N+1 rows), whereas the archived files have 5000 rows indexed `0,...,4999`.

| Nominal E | Archived fit xi (5000 rows) | Archived R² | New fit xi (N=5000) | New R² |
|---|---:|---:|---:|---:|
| -1 | 1.532809239e9 | 0.00005469 | 3.463409703e9 | 0.00001137 |
| +1 | 2.711107146e9 | 0.00001790 | 1.785420784e9 | 0.00004029 |

**Do not interpret the large fitted xi values as measured physical localization lengths.** Here `R²` is effectively zero, and the fitted slopes are dominated by bounded oscillations rather than a linear trend. Unlike the Binary (3,6) case, matching N does *not* make the huge fitted xi values match: that is expected because dividing by a nearly zero, unstable slope magnifies tiny changes.

The raw mean-log-transmission curves agree far better than their fitted xi values suggest:

| Nominal E | Archived mean log T over sites 500..4999 | New N=5000 mean over same sites | Pointwise RMSE over 5000 sites |
|---|---:|---:|---:|
| -1 | -0.002085506 | -0.002085894 | 0.000080712 |
| +1 | -0.002384800 | -0.002384830 | 0.000091143 |

Increasing the new chain length by a factor of twenty does not cause mean log transmission to decay:

| E | New mean log T at N=5000 | at N=10000 | at N=100000 |
|---|---:|---:|---:|
| -1 | -0.002179633 | -0.002135130 | -0.002076564 |
| +1 | -0.002339023 | -0.002383970 | -0.002319550 |

The N=100000 linear fit even has a **positive** tiny slope at E=-1; this further illustrates why `-2/slope` is not an appropriate estimator in this special case. The numerical values in the tables are directly computed from the checked-in raw traces and summary files.

## Exact transfer-matrix explanation (ideal model)

Use the free-site and impurity-site transfer matrices

```math
F(E)=\begin{pmatrix}
E & -1\\
1 & 0
\end{pmatrix},
\qquad
A(E,U)=\begin{pmatrix}
E-U & -1\\
1 & 0
\end{pmatrix}.
```

At `E=+1`, `F^3=-I`; at `E=-1`, `F^3=+I`. Since the two allowed spacings differ by three sites, the matrices advancing from one impurity to the next are

```math
B_2=AF,\qquad B_5=AF^4=\pm AF=\pm B_2.
```

More specifically, `B_5=-B_2` at `E=+1`, while `B_5=B_2` at `E=-1`.

Thus, apart from bounded terminal factors, the random product of impurity-to-impurity matrices equals `+/- (B_2)^m`. In this convention,

```math
\det B_2=1,\qquad \mathrm{tr}(B_2)=E(E-U)-2.
```

For `U=0.1`, the trace of `B_2` is `-1.1` at `E=+1` and `-0.9` at `E=-1`. Both lie strictly between `-2` and `2`, so the eigenvalues of `B_2` are distinct, lie on the unit circle, and its powers are bounded.

Therefore, the log transmission does not grow linearly in magnitude with N, and the **asymptotic Lyapunov exponent is zero** at the two *exact* resonance energies in this idealized model. The asymptotic localization length is divergent, rather than a finite number such as `1e9`.

This transfer-matrix argument relies on the precise spacings, energy, potential, and nearest-neighbor model. It does not prove delocalization over a finite energy interval or under arbitrary model modifications.
## Files and provenance

- `historical_traces/`: original uploaded raw mean-log-transmission files, copied without changing their contents.
- `new_runs/N5000_E*/`, `N10000_E*/`, `N100000_E*/`: raw C++ traces, `summary.csv`, and `run_config.txt` for each run.
- `comparison_trace_E*.csv`: aligned original and new N=5000 traces and their pointwise difference.
- `overlay_E*.png`, `difference_E*.png`: graphical comparisons.
- `resonance_summary.csv`: compact numerical summary.
- `check_fits.py`: independent re-fit of all raw traces and exact-matrix identity check using only Python's standard library.

The original raw traces do not record the historical RNG seed, exact simulation invocation, or the number of realizations. The nominal energy naming also cannot establish that the historical energy was exactly +/-1. These uncertainties are preserved rather than guessed.

## Recompute the checked-in results

From the repository root:

```bash
python validation/resonance_binary_2_5/check_fits.py
```

This recomputes numerical fits and checks the transfer-matrix identities. It does not rerun the Monte Carlo simulation.

To regenerate the six C++ runs on Linux/macOS, first build and test the project as in the root README, then use the commands below (choose fresh output directories):

```bash
./build/quantum_transport --model binary --N 5000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E -1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N5000_Eminus1

./build/quantum_transport --model binary --N 5000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E 1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N5000_Eplus1

./build/quantum_transport --model binary --N 10000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E -1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N10000_Eminus1

./build/quantum_transport --model binary --N 10000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E 1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N10000_Eplus1

./build/quantum_transport --model binary --N 100000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E -1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N100000_Eminus1

./build/quantum_transport --model binary --N 100000 --U 0.1 --s1 2 --s2 5 --p 0.5 --E 1 --realizations 1000 --threads 1 --seed 12345 --out check_binary25_N100000_Eplus1
```

On Windows, use the appropriate executable path (often `build/Release/quantum_transport.exe`) and the Windows CTest option described in the root README. Output directories should be new, or `--overwrite` must be supplied intentionally.

## Interpretation and limits

This check independently explains the apparently enormous and inconsistent fitted xi values at the two exact resonances of the ideal Binary (2,5) model. It supports the implementation in this specific controlled case, but **does not verify all energies, all disorder models, the entire thesis methodology, or finite-size scaling away from the exact resonances**. The archived numerical data have not been modified.
