#!/usr/bin/env python3
"""Compare numerical Lyapunov exponents with the second-order formula."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_config(folder):
    entries = {}
    for raw in (folder / "run_config.txt").read_text().splitlines():
        if "=" in raw:
            key, value = raw.split("=", 1)
            entries[key] = value
    return entries


def spacing_distribution(c):
    model = c["model"]
    if model == "binary":
        return np.array([int(c["s1"]), int(c["s2"])]), np.array([
            float(c["p"]), 1.0 - float(c["p"])
        ])
    rank = np.arange(1, int(c["max_rank"]) + 1, dtype=float)
    if model == "exponential":
        support = rank
        log_weights = -float(c["mu"]) * (rank - 1)
    else:
        support = int(c["s0"]) + (rank - 1) * int(c["d"])
        base = support if model == "zeta-value" else rank
        log_weights = -(1.0 + float(c["alpha"])) * np.log(base)
    log_weights -= log_weights.max()
    weights = np.exp(log_weights)
    return support, weights / weights.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    args = ap.parse_args()
    c = read_config(args.run)
    support, weights = spacing_distribution(c)
    avg_s = float(np.dot(support, weights))
    with (args.run / "summary.csv").open(newline="") as f:
        data = list(csv.DictReader(f))
    out = args.run / "comparison"
    out.mkdir(exist_ok=True)
    numeric = np.array([float(row["gamma"]) for row in data])
    energies = np.array([float(row["energy"]) for row in data])
    gamma2 = np.full_like(energies, np.nan)
    U = float(c["U"])
    for i, E in enumerate(energies):
        k = float(np.arccos(E / 2))
        a2 = np.dot(weights, np.exp(-2j * k * support))
        denominator = abs(1.0 - a2) ** 2
        if denominator > 1e-12 and abs(np.sin(k)) > 1e-8:
            gamma2[i] = (U * U / (8.0 * avg_s * np.sin(k) ** 2)
                         * (1.0 - abs(a2) ** 2) / denominator)
    with (out / "gamma_comparison.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["energy", "gamma_numerical", "gamma_2nd_order"])
        writer.writerows(zip(energies, numeric, gamma2))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(energies, numeric, ".", label="Numerical")
    ax.plot(energies, gamma2, "-", label="Second-order perturbation")
    ax.set(xlabel="Energy / hopping", ylabel=r"Lyapunov exponent $\gamma$",
           title="Numerical and second-order results")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "gamma_comparison.png", dpi=170)
    plt.close(fig)
    print(f"Comparison saved to: {out}")
    print("The expansion may fail near band edges, resonances, or at large |U|.")


if __name__ == "__main__":
    main()
