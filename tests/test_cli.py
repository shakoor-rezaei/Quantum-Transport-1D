#!/usr/bin/env python3
"""Check CLI output, reproducibility, and basic data consistency."""

import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path


exe = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path("build/quantum_transport").resolve()


def run(model, folder, threads=1, extra=()):
    subprocess.run([str(exe), "--model", model, "--N", "60", "--U", "0.2",
                    "--realizations", "15", "--emin", "-0.5", "--emax", "0.5",
                    "--dE", "0.5", "--threads", str(threads), "--out", str(folder),
                    *extra], check=True, capture_output=True, text=True)
    with (folder / "summary.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    assert all(math.isfinite(float(row["mean_lnT_N"])) for row in rows)
    for row in rows:
        trace = [float(x) for x in (folder / row["file"]).read_text().split()]
        assert len(trace) == 61
        assert math.isclose(trace[-1], float(row["mean_lnT_N"]), abs_tol=1e-12)
    return [((folder / row["file"]).read_bytes(), float(row["gamma"])) for row in rows]


with tempfile.TemporaryDirectory() as folder:
    root = Path(folder)
    for model in ("binary", "exponential", "zeta", "zeta-rank", "zeta-value"):
        extra = ("--s0", "1", "--d", "2", "--max-rank", "1000") if "zeta-" in model else (
            "--max-rank", "1000")
        a = run(model, root / (model + "_1"), 1, extra)
        b = run(model, root / (model + "_2"), 2, extra)
        assert a == b, f"Thread-dependent output: {model}"
    clean = run("binary", root / "clean", 1, ("--U", "0",))
    assert all(abs(gamma) < 1e-10 for _, gamma in clean)
print("CLI checks passed.")
