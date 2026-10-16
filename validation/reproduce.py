#!/usr/bin/env python3
"""Rerun the four standard cases and compare with archived MSc data."""

import argparse
import subprocess
import sys
from pathlib import Path


BASE = Path(__file__).resolve().parent
CASES = {
    "binary_2_5": ["--model", "binary", "--s1", "2", "--s2", "5", "--p", "0.5", "--U", "0.1"],
    "binary_3_6": ["--model", "binary", "--s1", "3", "--s2", "6", "--p", "0.5", "--U", "0.1"],
    "exponential": ["--model", "exponential", "--mu", "0.5", "--U", "0.2"],
    "zeta": ["--model", "zeta", "--alpha", "2", "--U", "0.2"],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--runs-dir", type=Path, default=BASE / "work")
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.is_file():
        parser.error(f"Cannot find executable: {exe}")
    if args.threads < 1:
        parser.error("threads must be positive")
    for name, model_args in CASES.items():
        dest = args.runs_dir.resolve() / name
        argv = [str(exe), *model_args, "--N", "10000", "--realizations", "1000",
                "--emin", "-1.9", "--emax", "1.9", "--dE", "0.05",
                "--seed", "12345", "--threads", str(args.threads),
                "--max-rank", "100000", "--out", str(dest), "--overwrite"]
        print("Running", name, flush=True)
        subprocess.run(argv, check=True)
    subprocess.run([sys.executable, str(BASE / "compare_thesis.py"),
                    "--runs-dir", str(args.runs_dir.resolve()),
                    "--output-dir", str(args.runs_dir.resolve() / "comparison")], check=True)


if __name__ == "__main__":
    main()
