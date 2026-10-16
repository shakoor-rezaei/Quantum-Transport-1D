#!/usr/bin/env python3
"""Recompute finite-length linear-fit xi and R^2 from archived resonance traces.

Uses only the Python standard library. No simulation is performed here.
"""
from pathlib import Path
import csv
import math

ROOT = Path(__file__).resolve().parent


def fit(values):
    n = len(values)
    xbar = (n - 1) / 2.0
    ybar = sum(values) / n
    xx = sum((i - xbar) ** 2 for i in range(n))
    xy = sum((i - xbar) * (v - ybar) for i, v in enumerate(values))
    slope = xy / xx
    intercept = ybar - slope * xbar
    total = sum((v - ybar) ** 2 for v in values)
    residual = sum((v - slope * i - intercept) ** 2 for i, v in enumerate(values))
    r2 = 1.0 - residual / total if total else float("nan")
    return -2.0 / slope, r2


def load_trace(path):
    return [float(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def check(actual, expected, label, tolerance=5e-6):
    if not math.isclose(actual, expected, abs_tol=tolerance, rel_tol=1e-9):
        raise AssertionError(f"{label}: got {actual:.12g}, expected {expected:.12g}")


def main():
    with (ROOT / "resonance_summary.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        e = int(row["energy"])
        print(f"\nNominal E={e:+d}")
        old = load_trace(ROOT / "historical_traces" / f"energy_{e}_pot_0.1.dat")
        assert len(old) == 5000
        xi, r2 = fit(old)
        check(xi, float(row["thesis_xi"]), "historical xi")
        check(r2, float(row["thesis_fit_r2"]), "historical R2")
        print(f"  historical: {len(old)} points, xi={xi:.6f}, R2={r2:.6f}")
        for N in (5000, 10000):
            trace = load_trace(ROOT / "new_runs" / f"N{N}_E{e:+d}" / "energy_00000.dat")
            assert len(trace) == N + 1
            xi, r2 = fit(trace)
            check(xi, float(row[f"new_xi_N{N}"]), f"new xi N={N}")
            check(r2, float(row[f"new_fit_r2_N{N}"]), f"new R2 N={N}")
            print(f"  new N={N}: {len(trace)} points, xi={xi:.6f}, R2={r2:.6f}")
    print("\nAll archived fit checks passed.")


if __name__ == "__main__":
    main()
