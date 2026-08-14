#!/usr/bin/env python3
"""
Non-rigorous generator for candidate intervals.

The computation here uses ordinary floating point arithmetic. Its output is a
starting point for certify.py, not a proof.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List

QUANTITIES = ("gamma", "beta_plus", "beta_minus")
MILESTONES = {10, 50, 100, 500, 1000, 2000, 5000}


def bisect_root(f) -> float:
    """Floating point bisection for a root of f, with f(0)<0."""
    lo, hi = 0.0, 1.0
    while f(hi) <= 0.0:
        hi *= 2.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if f(mid) <= 0.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def G(p: float, q: float) -> float:
    """Floating point value of the scalar function G(p,q)."""
    if p > q:
        p, q = q, p

    if p == 0.0 and q == 0.0:
        return 1.0

    if p == 0.0:
        sq = math.sqrt(q)
        u = bisect_root(lambda t: sq * t - q * math.log1p(t / sq) - 1.0)
        return q + sq * u

    sp = math.sqrt(p)
    sq = math.sqrt(q)
    spq = math.sqrt(p * q)

    def residual(u: float) -> float:
        eu = math.exp(u)
        ratio = (sp * eu + sq) / (sp + sq * eu)
        return 2.0 * spq * math.sinh(u) - (p + q) * u + (p - q) * math.log(ratio) - 1.0

    u = bisect_root(residual)
    return p + q + 2.0 * spq * math.cosh(u)


def symmetric_indices(k: int):
    for i in range(1, k // 2 + 1):
        yield i, 2
    if k % 2 == 1:
        yield k // 2 + 1, 1


def compute_gamma(K: int) -> List[float]:
    x = [0.0] * (K + 1)
    for k in range(1, K + 1):
        total = 0.0
        for i, weight in symmetric_indices(k):
            total += weight * G(x[i - 1], x[k - i])
        x[k] = total / k
        if k in MILESTONES:
            print(f"k={k:5d}, gamma[k]/k^2={x[k] / (k*k):.12f}", flush=True)
    return x


def compute_beta_plus(K: int) -> List[float]:
    x = [0.0] * (K + 1)
    for k in range(1, K + 1):
        best = 0.0
        for i, _weight in symmetric_indices(k):
            best = max(best, G(x[i - 1], x[k - i]))
        x[k] = best
        if k in MILESTONES:
            print(f"k={k:5d}, beta_plus[k]/k^2={x[k] / (k*k):.12f}", flush=True)
    return x


def compute_beta_minus(K: int, full_scan_prefix: int = 300, full_scan_all: bool = False) -> List[float]:
    """
    Compute beta_minus candidate values.

    For the generated K=5000 file, the minimizing split is i=2 in the checked
    prefix and is then extrapolated. This is only a way to propose candidates;
    certify.py verifies the recurrence rigorously from the intervals. Use
    --beta-minus-full-scan to compute the floating point minimum over all splits.
    """
    x = [0.0] * (K + 1)
    for k in range(1, K + 1):
        if k == 1:
            x[k] = 1.0
        elif full_scan_all or k <= full_scan_prefix:
            best = float("inf")
            best_i = None
            for i, _weight in symmetric_indices(k):
                value = G(x[i - 1], x[k - i])
                if value < best:
                    best = value
                    best_i = i
            x[k] = best
            if k >= 3 and best_i != 2:
                print(f"warning: beta_minus floating minimizer at k={k} was i={best_i}, not i=2", flush=True)
        else:
            # Symmetry gives the same value for i=k-1.
            x[k] = G(x[1], x[k - 2])

        if k in MILESTONES:
            print(f"k={k:5d}, beta_minus[k]/k^2={x[k] / (k*k):.12f}", flush=True)
    return x


def values_to_intervals(values: List[float], rel_pad: float, abs_pad: float) -> List[dict]:
    intervals = []
    for k, value in enumerate(values):
        radius = abs_pad + rel_pad * k * k
        intervals.append({
            "lo": f"{max(0.0, value - radius):.17g}",
            "hi": f"{value + radius:.17g}",
        })
    intervals[0] = {"lo": "0", "hi": "0"}
    if len(intervals) > 1:
        intervals[1] = {"lo": "1", "hi": "1"}
    return intervals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--K", type=int, required=True)
    parser.add_argument("--out", default="intervals.json")
    parser.add_argument("--quantity", choices=("all",) + QUANTITIES, default="all")
    parser.add_argument("--rel-pad", default="1e-9", help="interval radius includes rel_pad*k^2")
    parser.add_argument("--beta-minus-rel-pad", default="1e-8", help="beta_minus interval radius includes beta_minus_rel_pad*k^2")
    parser.add_argument("--abs-pad", default="1e-12", help="interval radius includes this constant")
    parser.add_argument("--beta-minus-full-scan", action="store_true", help="use the O(K^2) floating scan for beta_minus")
    parser.add_argument("--beta-minus-full-scan-prefix", type=int, default=300, help="prefix checked before the fast beta_minus extrapolation")
    args = parser.parse_args()

    rel_pad = float(args.rel_pad)
    beta_minus_rel_pad = float(args.beta_minus_rel_pad)
    abs_pad = float(args.abs_pad)

    requested = QUANTITIES if args.quantity == "all" else (args.quantity,)
    intervals: Dict[str, List[dict]] = {}

    for quantity in requested:
        print(f"computing {quantity}", flush=True)
        if quantity == "gamma":
            values = compute_gamma(args.K)
        elif quantity == "beta_plus":
            values = compute_beta_plus(args.K)
        else:
            values = compute_beta_minus(
                args.K,
                full_scan_prefix=args.beta_minus_full_scan_prefix,
                full_scan_all=args.beta_minus_full_scan,
            )
        quantity_rel_pad = beta_minus_rel_pad if quantity == "beta_minus" else rel_pad
        intervals[quantity] = values_to_intervals(values, quantity_rel_pad, abs_pad)

    output = {
        "K": args.K,
        "intervals": intervals,
        "metadata": {
            "note": "generated by floating point code; certify.py performs the rigorous check",
            "radius": "abs_pad + rel_pad*k^2",
            "abs_pad": args.abs_pad,
            "rel_pad": args.rel_pad,
            "beta_minus_rel_pad": args.beta_minus_rel_pad,
            "beta_minus_generation": "fast i=2 extrapolation after a floating full-scan prefix unless --beta-minus-full-scan is used",
        },
    }

    Path(args.out).write_text(json.dumps(output, indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
