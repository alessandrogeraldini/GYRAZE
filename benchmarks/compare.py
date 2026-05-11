#!/usr/bin/env python3
"""
Compare phi_n_DS.txt and phi_n_MP.txt from a GYRAZE run against reference files.

Usage:
    python3 compare.py <run_dir> <ref_dir>

    run_dir  — directory containing phi_n_DS.txt and phi_n_MP.txt from the run
    ref_dir  — directory containing the reference phi_n_DS.txt and phi_n_MP.txt
               (e.g. tests/pwall5.0 or tests/jwall0.0)

Exits 0 and prints "correct" if all values agree within 1% relative tolerance.
Exits 1 and lists every failing value otherwise.
"""

import sys
from pathlib import Path

TOL = 0.01  # 1% relative tolerance
ABS_FLOOR = 1e-10  # avoid divide-by-zero for near-zero reference values


def load(path: Path) -> list[list[float]]:
    rows = []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError as e:
            sys.exit(f"ERROR: parse error on line {i} of {path}: {e}")
    return rows


def compare_file(run_path: Path, ref_path: Path, label: str) -> list[str]:
    run = load(run_path)
    ref = load(ref_path)

    failures = []

    if len(run) != len(ref):
        failures.append(
            f"{label}: row count mismatch — got {len(run)}, ref {len(ref)}"
        )
        return failures

    for row_i, (run_row, ref_row) in enumerate(zip(run, ref)):
        if len(run_row) != len(ref_row):
            failures.append(
                f"{label} row {row_i}: column count mismatch "
                f"— got {len(run_row)}, ref {len(ref_row)}"
            )
            continue
        for col_j, (got, expected) in enumerate(zip(run_row, ref_row)):
            scale = max(abs(expected), ABS_FLOOR)
            if abs(got - expected) / scale > TOL:
                failures.append(
                    f"{label} row {row_i} col {col_j}: "
                    f"got {got:.8g}, ref {expected:.8g}, "
                    f"rel err {abs(got - expected) / scale:.3%}"
                )

    return failures


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("Usage: compare.py <run_dir> <ref_dir>")

    run_dir = Path(sys.argv[1])
    ref_dir = Path(sys.argv[2])

    for d, label in [(run_dir, "run_dir"), (ref_dir, "ref_dir")]:
        if not d.is_dir():
            sys.exit(f"ERROR: {label} '{d}' is not a directory")

    failures = []
    for fname in ("phi_n_DS.txt", "phi_n_MP.txt"):
        run_file = run_dir / fname
        ref_file = ref_dir / fname
        for f, name in [(run_file, "run"), (ref_file, "ref")]:
            if not f.exists():
                sys.exit(f"ERROR: {name} file not found: {f}")
        failures.extend(compare_file(run_file, ref_file, fname))

    if not failures:
        print("correct")
    else:
        print(f"FAILED ({len(failures)} value(s) outside {TOL:.0%} tolerance):")
        for msg in failures:
            print(f"  {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
