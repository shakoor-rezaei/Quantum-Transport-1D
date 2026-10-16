# Quantum transport in disordered one-dimensional chains

A reproducible numerical study of transmission and localization in the one-dimensional nearest-neighbor tight-binding model. Impurity positions are generated from a probability distribution for the integer spacing between successive impurities.

This repository is based on the MSc research project **Quantum Transport in Lévy Glass** (Shakoor Rezaei, University of Tehran, 2025; supervisor: Dr Reza Sepehrinia). The three standard distributions reproduce the model families studied in [Sepehrinia, arXiv:2210.01377](https://arxiv.org/abs/2210.01377). The two generalized zeta options extend the choice of allowed spacings and their weights. They are not presented as results of that article.

## Model

With hopping set to one, the discrete Schrödinger equation is

```math
\psi_{n+1}+\psi_{n-1}+U_n\psi_n=E\psi_n.
```

A sampled spacing $s$ produces one block $[0,\ldots,0,U]$, containing $s-1$ zero-potential sites and one impurity. Blocks are appended until the chain contains exactly $N$ sites; the last block is truncated if needed. Leads are clean. The simulation uses energies $E\in(-2,2)$.

| `--model` | Allowed spacings | Probability rule |
|---|---|---|
| `binary` | $s_1,s_2$ | $P(s_1)=p,\ P(s_2)=1-p$ |
| `exponential` | $1,2,3,\ldots$ | $P(s)\propto\exp[-\mu(s-1)]$ |
| `zeta` | $1,2,3,\ldots$ | $P(s)\propto s^{-(1+\alpha)}$ |
| `zeta-rank` | $s_r=s_0+(r-1)d$ | $P(s_r)\propto r^{-(1+\alpha)}$ |
| `zeta-value` | $s_r=s_0+(r-1)d$ | $P(s_r)\propto s_r^{-(1+\alpha)}$ |

**Important:** `d` is the *arithmetic spacing increment*, not a geometric ratio. For `s0=1, d=2` the allowed spacings are `1, 3, 5, 7, ...`. The spacing `3` has weight $2^{-(1+\alpha)}$ in `zeta-rank`, but $3^{-(1+\alpha)}$ in `zeta-value`. Each law is normalized separately.

All infinite-support laws are truncated and renormalized at `--max-rank` (default `100000`). For small `alpha`, especially `alpha <= 1`, finite cutoff and finite-size effects can be substantial. Run cutoff and system-size convergence studies before interpreting asymptotic behavior. The reported mean spacing is for the truncated distribution.

## Numerical observables

For each energy, a seeded set of disorder realizations yields the mean log transmission $\langle\ln T_n\rangle$ after `n = 0,...,N` sites. With the thesis convention,

```math
\langle\ln T_n\rangle\simeq C-2\gamma n,\qquad \xi=\gamma^{-1}=-2/\text{slope}.
```

Fitting bounds `--fit-start` and `--fit-end` are fractions of `N` (default: the whole trace). Near-zero slopes are reported as `xi=inf`; a significantly positive slope is reported as `xi=nan`, since a negative inverse localization length is not a physical conclusion. Very large finite-chain `xi` alone does **not** establish a transparent state in the thermodynamic limit.

Units: energy and impurity potential are in units of hopping; $\xi$ is in lattice spacings. No eV or Å conversion is implied.

## Build

Requires a C++17 compiler, CMake >= 3.16, and Python 3 with NumPy and Matplotlib.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
ctest --test-dir build -C Release --output-on-failure
python -m pip install -r requirements.txt
```

On Windows, the executable may be `build/Release/quantum_transport.exe`; on Linux/macOS or a single-configuration generator it is generally `build/quantum_transport`.

Direct compile on Linux/macOS if CMake is unavailable:

```bash
g++ -O3 -std=c++17 -pthread -Iinclude src/main.cpp src/distributions.cpp src/transport.cpp -o quantum_transport
```

## Quick start

```bash
./build/quantum_transport --model binary --s1 2 --s2 5 --p 0.5 \
  --N 600 --U 0.05 --realizations 200 --emin -1.5 --emax 1.5 \
  --dE 0.5 --threads 4 --out results_binary

python analysis/plot_results.py results_binary --energy 0.5
python analysis/compare_second_order.py results_binary
```

`results_binary` must be new, or add `--overwrite`. Use separate output directories for each parameter set.

Standard exponential model:

```bash
./build/quantum_transport --model exponential --mu 0.5 --N 1000 \
  --U 0.2 --realizations 1000 --out results_exponential
```

Standard zeta model:

```bash
./build/quantum_transport --model zeta --alpha 2 --N 1000 \
  --U 0.2 --realizations 1000 --out results_zeta
```

Generalized zeta with rank weights:

```bash
./build/quantum_transport --model zeta-rank --alpha 2 --s0 1 --d 2 \
  --N 1000 --U 0.2 --realizations 1000 --out results_rank
```

Generalized zeta with spacing-value weights:

```bash
./build/quantum_transport --model zeta-value --alpha 2 --s0 1 --d 2 \
  --N 1000 --U 0.2 --realizations 1000 --out results_value
```

For a single energy, replace the energy-grid options with `--E 1`. Run `--help` for all options.

## Files produced

- `run_config.txt`: model, simulation parameters, seed and truncated mean spacing.
- `distribution_preview.csv`: first 40 supported spacings and their normalized probabilities.
- `sample_chain.csv`: first 500 sites of one realization.
- `energy_00000.dat`, etc.: $N+1$ values of $\langle\ln T_n\rangle$, starting at $n=0$.
- `summary.csv`: energies, final $\langle\ln T_N\rangle$, slope, Lyapunov exponent, localization length and fit quality.
- `plots/`: graphs from `plot_results.py`.
- `comparison/`: second-order comparison from `compare_second_order.py`.


**GitHub upload note:** To keep the browser-upload package small, precomputed per-energy
trace files in `examples/demo_*` are omitted. The example `summary.csv` files,
parameter configurations and plots remain available. Detailed traces can be
regenerated by running the C++ simulator with the documented parameters.
All thesis reference data and validation comparisons are retained.

The independent random seed of a realization does not depend on energy or thread count; identical parameters and seed reproduce the same output in different thread configurations. The C++ code computes transmission; Python reads the existing output without resampling any chains.

## Reproducing the MSc numerical results

The archived MSc numerical output and historical plots are preserved under
[`validation/`](validation/README.md), with SHA-256 checksums in the report.
A separate execution of this C++ version was performed for the four standard
models at `N=10000`, `realizations=1000`, `seed=12345`, `dE=0.05`,
`E=-1.9..1.9`. The 77 generated energies were matched to the 400 archived
energy points; these are **not claimed to be bitwise-identical stochastic runs**.

| Model | Paired energies | Median absolute relative difference, selected region |
|---|---:|---:|
| Binary 2,5 | 77 | 2.90% |
| Binary 3,6 | 77 | 3.00% |
| Exponential | 77 | 1.21% |
| Standard zeta | 77 | 1.30% |

Comparison-region rules, 90th percentiles, individual CSVs and the resonance
exceptions are given in [`validation/figures/comparison_report.txt`](validation/figures/comparison_report.txt).
In particular, the old and new binary `(3,6)` finite-chain fitted values at
`E=±1` are not equal, and the binary `(2,5)` huge values at those energies have
very poor fit quality. The agreement figures above intentionally **do not**
include those singular neighborhoods. No unverified assertion of agreement at
all energies is made.

Reproduce the four runs and comparison locally (after building):

```bash
python validation/reproduce.py --exe build/quantum_transport --threads 4
```

On Windows, the executable may instead be
`build/Release/quantum_transport.exe`.

## Validation and limitations

`ctest` checks clean-chain transmission, the exact one-impurity result, distribution-weight ratios, the standard-zeta limit of both generalized schemes, and finite numerical output. The additional CLI test checks deterministic results across thread counts and consistency between the traces and summary:

```bash
python tests/test_cli.py build/quantum_transport
```

`analysis/compare_second_order.py` implements the **second-order**, not fourth-order, weak-disorder prediction of the cited paper:

```math
\gamma_2(E)=\frac{U^2}{8\langle s\rangle\sin^2 k}\,
\frac{1-|A_2|^2}{|1-A_2|^2},\qquad
A_2=\sum_s P(s)e^{-2iks},\quad E=2\cos k.
```

Its probability sum uses the same finite cutoff as the simulator. This approximation may fail near band edges, resonant/anomalous energies and outside weak disorder. It must not be used as a universal exact benchmark.

For comparison with the article, relevant published parameter sets include binary `(s1,s2,p,U)=(2,5,0.5,0.1)` and `(3,6,0.5,0.1)`, exponential `(mu,U)=(0.5,0.2)`, and zeta `(alpha,U)=(2,0.2)`. **The small demo in `examples/demo_binary` is an illustrative numerical check, not a reproduction of those full published figures.** Before presenting a figure as reproduced, verify its system size, averaging count, fit interval, cutoff and uncertainties.

The MSc implementation and subsequent revisions were the starting points for this refactor. The refactor improves reproducibility and numerical range; it does not, by itself, establish a new physical result or prove correctness for every parameter regime.

## Reference

R. Sepehrinia, *Localization in a one-dimensional alloy with an arbitrary distribution of spacing between impurities: Application to Lévy glass*, arXiv:2210.01377 (2022).
