#!/usr/bin/env python3
"""Plot localization length and mean log transmission from a simulation run."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path, help="Simulation output directory")
    ap.add_argument("--energy", type=float, help="Draw <ln T_n> for the nearest energy")
    ap.add_argument("--ymax", type=float, help="Optional visible upper limit for xi(E)")
    args = ap.parse_args()
    with (args.run / "summary.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("Empty summary.csv")
    energies = np.array([float(r["energy"]) for r in rows])
    xi = np.array([float(r["xi"]) for r in rows])
    means = np.array([float(r["mean_lnT_N"]) for r in rows])
    out = args.run / "plots"
    out.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(energies, means, ".-")
    ax.set(xlabel="Energy / hopping", ylabel=r"$\langle\ln T_N\rangle$",
           title="Mean log transmission after N sites")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "mean_lnT_N_vs_energy.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    finite = np.isfinite(xi) & (xi > 0)
    ax.plot(energies[finite], xi[finite], ".-")
    if args.ymax is not None:
        if args.ymax <= 0:
            raise SystemExit("--ymax must be positive")
        ax.set_ylim(0, args.ymax)
        ax.set_title("Localization length (display range clipped)")
    else:
        ax.set_title("Localization length")
    ax.set(xlabel="Energy / hopping", ylabel=r"$\xi$ (lattice sites)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "localization_length_vs_energy.png", dpi=170)
    plt.close(fig)

    if args.energy is not None:
        idx = int(np.argmin(np.abs(energies - args.energy)))
        row = rows[idx]
        y = np.loadtxt(args.run / row["file"], ndmin=1)
        n = int(row["N"])
        if len(y) != n + 1:
            raise ValueError("Trace length does not match N+1")
        if not np.isclose(y[-1], float(row["mean_lnT_N"]), rtol=1e-10, atol=1e-12):
            raise ValueError("Trace and summary disagree")
        x = np.arange(n + 1)
        start, end = int(row["fit_start_site"]), int(row["fit_end_site"])
        slope, gamma = float(row["slope"]), float(row["gamma"])
        intercept = float(np.mean(y[start:end + 1] - slope * x[start:end + 1]))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(x, y, label=r"$\langle\ln T_n\rangle$")
        ax.plot(x[start:end + 1], slope * x[start:end + 1] + intercept,
                label=f"Fit; gamma={gamma:.5g}")
        ax.set(xlabel="Number of sites n", ylabel=r"$\langle\ln T_n\rangle$",
               title=f"Energy={energies[idx]:.8g}")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out / "mean_lnT_vs_sites.png", dpi=170)
        plt.close(fig)
    print(f"Plots saved to: {out}")


if __name__ == "__main__":
    main()
