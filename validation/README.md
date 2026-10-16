# MSc numerical-results archive and reproducibility check

The four files in `thesis_data/` are unmodified copies of the submitted two-column
energy/localization-length files (the figures in `thesis_figures/` are the supplied
historical plots). The file `unclassified/results_unidentified.dat` is preserved,
but its model and parameters have not been established; it is not used for validation.

The independent reruns in `reproduced/` use the current C++ source. They include
configuration files, complete energy summaries, and a representative mean-log-
transmission trace at E=0. The large full set of individual traces is excluded from
the GitHub package, but can be regenerated using `reproduce.py`.

Run from the repository root, after building:

```bash
python validation/reproduce.py --exe build/quantum_transport --threads 4
```

On Windows, use `--exe build/Release/quantum_transport.exe` when applicable.
The script writes fresh outputs into `validation/work/` and then compares them
with the unchanged thesis files. The checked-in comparison from the verified run
is in `figures/comparison_report.txt`, together with the CSV files and charts.

## Interpretation

The four cases match the article's model parameters (Sepehrinia, 2022,
Figs. 2, 3, 5, 7). These calculations are *not* the author's fourth-order
perturbation curves; they are simulations compared with the student's earlier
numerical results. The separate second-order check can be generated with
`analysis/compare_second_order.py` on a complete C++ run. The 2022 paper also
shows numerical and fourth-order comparisons, which remain a distinct benchmark.

The old files have 400 energy rows, including E=-2 where the clean-lead formula
is singular; the new run evaluates 77 points on -1.9..1.9 at dE=0.05.
The original RNG seeds, exact fit interval and averaging counts are not encoded
in those two-column files. The thesis describes N=10000 for the four relevant
figures and 5000/10000 realizations elsewhere, but does not associate each
archived output with a uniquely recovered run configuration.

The comparison statistics are deliberately restricted to finite, positive
lengths and |E| <= 1.75; near E=+-1 the binary cases are treated separately.
The binary (2,5) peaks are extremely sensitive to near-zero fitted slopes, and
the two binary (3,6) points at E=+-1 differ markedly from the old archive.
Those differences are visible in the full comparison CSVs. Agreement away from
singular regions is evidence for the implementation in those tested conditions,
not universal proof of accuracy, exact divergence or new discoveries.
