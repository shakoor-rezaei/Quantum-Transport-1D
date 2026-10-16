#!/usr/bin/env python3
"""Compare archived MSc results with a separate reproducible C++ run."""

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASES = {
    "binary_2_5": "Binary spacings 2 and 5",
    "binary_3_6": "Binary spacings 3 and 6",
    "exponential": "Exponential spacing",
    "zeta": "Standard zeta spacing",
}
ROOT = Path(__file__).resolve().parent


def read_run(path):
    with (path / "summary.csv").open(newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    return {round(float(row["energy"]), 7): row for row in records}


def comparison(name, label, runs_dir, figures_dir):
    source = ROOT / "thesis_data" / (name + ".dat")
    old = np.loadtxt(source)
    run = read_run(runs_dir / name)
    records = []
    for e, thesis_xi in old:
        row = run.get(round(float(e), 7))
        if row is None:
            continue
        xi_new = float(row["xi"])
        fit_r2 = float(row["r_squared"])
        eligible = (np.isfinite(thesis_xi) and np.isfinite(xi_new)
                    and thesis_xi > 0 and xi_new > 0)
        # The summary metric omits band edges and the narrow binary resonances.
        stable = (eligible and abs(e) <= 1.75 and
                  ("binary" not in name or abs(abs(e) - 1) >= 0.1))
        delta = abs(xi_new - thesis_xi) / thesis_xi if eligible else float("nan")
        records.append((float(e), thesis_xi, xi_new, fit_r2, delta, stable))
    if not records:
        raise RuntimeError(f"No common energies for {name}")

    all_data = np.array(records, dtype=float)
    selected = all_data[all_data[:, 5] > 0.5]
    if len(selected) < 10:
        raise RuntimeError(f"Too few comparison points for {name}")
    errs = selected[:, 4]
    with (figures_dir / (name + "_comparison.csv")).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["energy", "thesis_xi", "new_xi", "new_fit_r2",
                         "absolute_relative_difference", "included_in_summary"])
        for e, ref, got, r2, diff, stable in records:
            writer.writerow([f"{e:.8g}", ref, got, r2, diff, int(stable)])

    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    mask = np.isfinite(old[:, 1]) & (old[:, 1] > 0)
    # Hide extreme resonance amplitudes on the display only, never in the CSV.
    values = old[mask & (abs(old[:, 0]) <= 1.75), 1]
    if "binary" in name:
        values = old[mask & (abs(old[:, 0]) <= 1.75) &
                     (abs(abs(old[:, 0]) - 1) >= 0.1), 1]
    ymax = np.percentile(values, 98) * 1.3
    ax.plot(old[:, 0], old[:, 1], label="Archived MSc data", linewidth=1.2)
    ax.plot(all_data[:, 0], all_data[:, 2], ".", markersize=4,
            label="Reproduced C++ run")
    ax.set(xlabel="Energy / hopping", ylabel="Localization length / sites",
           title=label, xlim=(-2, 2), ylim=(0, ymax))
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / (name + "_overlay.png"), dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    ax.plot(selected[:, 0], 100 * errs, ".-")
    ax.set(xlabel="Energy / hopping", ylabel="Absolute difference (%)",
           title=label + " — selected finite-chain comparison")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures_dir / (name + "_relative_difference.png"), dpi=170)
    plt.close(fig)
    return {
        "name": name, "label": label,
        "archived_points": len(old), "new_points": len(run),
        "paired_points": len(all_data), "selected_points": len(selected),
        "median_percent": float(np.median(errs) * 100),
        "p90_percent": float(np.percentile(errs, 90) * 100),
        "maximum_percent": float(np.max(errs) * 100),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "resonance": [(r[0], r[1], r[2], r[3]) for r in records
                      if "binary" in name and abs(abs(r[0]) - 1) < 1e-7],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-dir", type=Path, default=ROOT / "reproduced")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = [comparison(name, label, args.runs_dir, args.output_dir)
               for name, label in CASES.items()]
    lines = [
        "Archived thesis results versus an independent C++ rerun",
        "",
        "Rerun: N=10000, 1000 realizations/energy, seed=12345, dE=0.05,",
        "energies -1.9..1.9. Archived samples span -2..1.99 with dE=0.01.",
        "Other archived run metadata (seeds, fit intervals, sample counts)",
        "are not recoverable from the two-column files alone.",
        "",
        "Comparison domain: |E|<=1.75, finite positive xi; binary cases also",
        "exclude ||E|-1|<0.1 from summary statistics. All paired values",
        "including singular energies are retained in comparison CSV files.",
        "These differences measure agreement between finite Monte Carlo runs,",
        "not a rigorous confidence interval or proof of thermodynamic behavior.",
        "",
    ]
    for r in results:
        lines.extend([
            f"{r['label']} ({r['name']})",
            f"  Archived/new/paired/selected points: {r['archived_points']}/"
            f"{r['new_points']}/{r['paired_points']}/{r['selected_points']}",
            f"  Median |delta xi|/xi_thesis: {r['median_percent']:.2f}%",
            f"  90th percentile: {r['p90_percent']:.2f}%",
            f"  Maximum in selected domain: {r['maximum_percent']:.2f}%",
            f"  Archive SHA-256: {r['sha256']}",
        ])
        if r["resonance"]:
            for E, old, new, r2 in r["resonance"]:
                lines.append(f"  Resonance E={E:g}: old xi={old:.8g}, "
                             f"new xi={new:.8g}, new fit R^2={r2:.5g}")
        lines.append("")
    dest = args.output_dir / "comparison_report.txt"
    dest.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("Charts and tables:", args.output_dir)


if __name__ == "__main__":
    main()
