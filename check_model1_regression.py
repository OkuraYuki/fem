#!/usr/bin/env python3
"""Lightweight regression check for the known small model_1 case."""

import math
import re
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(1)


def parse_summary(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    out = {}
    for key in ("final_step", "final_max_diff", "dt", "nstep(max)"):
        m = re.search(rf"^{re.escape(key)}:\s+([^\n]+)$", text, flags=re.MULTILINE)
        if m:
            out[key] = m.group(1).strip()
    return out


def main() -> int:
    repo = Path(__file__).resolve().parent
    summary = repo / "output" / "model_1" / "analysis_summary.txt"
    tempa = repo / "output" / "model_1" / "tempa.dat010"

    if not summary.exists():
        fail(f"missing summary: {summary}")
    if not tempa.exists():
        fail(f"missing final step output: {tempa}")

    stats = parse_summary(summary)
    if int(stats.get("final_step", "0")) != 10:
        fail(f"unexpected final_step: {stats.get('final_step')}")

    final_max_diff = float(stats.get("final_max_diff", "nan"))
    if not (8.0e-4 <= final_max_diff <= 2.0e-3):
        fail(f"final_max_diff out of expected band: {final_max_diff:.6e}")

    tokens = tempa.read_text(encoding="utf-8").split()
    if not tokens:
        fail("tempa.dat010 is empty")

    vals = []
    for t in tokens:
        v = float(t.replace("D", "E"))
        if not math.isfinite(v):
            fail("non-finite value found in tempa.dat010")
        vals.append(v)

    vmin = min(vals)
    vmax = max(vals)
    if vmin < -1e-6 or vmax > 1.0 + 1e-6:
        fail(f"value range out of expected bounds: vmin={vmin:.6e}, vmax={vmax:.6e}")

    has_fixed_1v = any(abs(v - 1.0) < 1e-12 for v in vals)
    if not has_fixed_1v:
        fail("expected fixed 1.0V node was not found")

    print("[OK] model_1 regression check passed")
    print(f"      final_step=10, final_max_diff={final_max_diff:.6e}, vmin={vmin:.6e}, vmax={vmax:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
