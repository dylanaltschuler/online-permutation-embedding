#!/usr/bin/env python3
"""
Rigorous interval checker for the dynamic programs in the manuscript.

Input:  intervals.json, containing rational interval enclosures for
        gamma_k, beta_k^+, and beta_k^- for 0 <= k <= K.
Output: a pass/fail check of the recurrence inclusions, plus a JSON file
        slacks.json recording the certified recurrence slacks and the finite
        bounds on the limiting constants, including the lower bound on c_typ
        obtained from the beta_k^- intervals.

The only non-standard package is python-flint. It gives access to Arb, which
performs ball arithmetic with rigorous outward rounding.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from flint import arb, ctx, fmpq

QUANTITIES = ("gamma", "beta_plus", "beta_minus")
MANUSCRIPT_GOALS = {
    "gamma": {"lower": None, "upper": "0.49967"},
    "beta_plus": {"lower": "0.50547", "upper": "0.50568"},
    "beta_minus": {"lower": None, "upper": "0.48867"},
}
CTYP_LOWER_CUTOFF = 312
CTYP_LOWER_GOAL = "0.48934"


# ---------------------------------------------------------------------------
# Reading exact decimal/rational input
# ---------------------------------------------------------------------------


def exact_number(x) -> fmpq:
    """Convert a JSON number/string to an exact rational number for FLINT."""
    r = Fraction(str(x).strip())
    return fmpq(r.numerator, r.denominator)


def ball(x) -> arb:
    """The Arb ball containing the exact number x."""
    return arb(exact_number(x))


def proved_nonnegative(x: arb) -> bool:
    """True only when Arb proves x >= 0."""
    return x.lower() >= 0


def proved_positive(x: arb) -> bool:
    """True only when Arb proves x > 0."""
    return x.lower() > 0


def stop(message: str):
    raise SystemExit("CERTIFICATION FAILED: " + message)


def arb_to_string(x):
    return x.str(40) if hasattr(x, "str") else x


# ---------------------------------------------------------------------------
# Certified bounds for the scalar function G(p,q)
# ---------------------------------------------------------------------------
#
# The recurrence uses
#
#   G(p,q) = min_{0 <= a < b <= 1}
#            (1 + p log(b/a) + q log((1-a)/(1-b))) / (b-a).
#
# The checker first finds an approximate scalar root by floating point
# bisection. This approximation is not trusted. It is only used to suggest a
# small interval [u_lo,u_hi]. Arb then verifies the required signs at the two
# endpoints, and only after those sign checks does the program use the bracket.


def boundary_residual(r: arb, u: arb) -> arb:
    sr = r.sqrt()
    return sr * u - r * (1 + u / sr).log() - 1


def interior_residual(p: arb, q: arb, u: arb) -> arb:
    sp = p.sqrt()
    sq = q.sqrt()
    eu = u.exp()
    ratio = (sp * eu + sq) / (sp + sq * eu)
    return 2 * (p * q).sqrt() * u.sinh() - (p + q) * u + (p - q) * ratio.log() - 1


def bisect_float_root(f) -> float:
    """Ordinary floating point bisection for a root of f, with f(0)<0."""
    lo, hi = 0.0, 1.0
    while f(hi) <= 0.0:
        hi *= 2.0
    for _ in range(45):
        mid = (lo + hi) / 2.0
        if f(mid) <= 0.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def approximate_u(p: arb, q: arb) -> float:
    """Floating point guess for the scalar root u."""
    pf = float(p)
    qf = float(q)

    if pf == 0.0 or qf == 0.0:
        r = max(pf, qf)
        sr = math.sqrt(r)
        return bisect_float_root(lambda u: sr * u - r * math.log1p(u / sr) - 1.0)

    sp = math.sqrt(pf)
    sq = math.sqrt(qf)
    spq = math.sqrt(pf * qf)

    def residual(u: float) -> float:
        eu = math.exp(u)
        ratio = (sp * eu + sq) / (sp + sq * eu)
        return 2.0 * spq * math.sinh(u) - (pf + qf) * u + (pf - qf) * math.log(ratio) - 1.0

    return bisect_float_root(residual)


def certified_u_bracket(p: arb, q: arb) -> Tuple[arb, arb]:
    """Find [u_lo,u_hi] and prove residual(u_lo)<=0<=residual(u_hi)."""
    u0 = approximate_u(p, q)
    width = max(1e-11, 1e-11 * abs(u0))

    if p == 0 or q == 0:
        r = q if p == 0 else p
        residual = lambda u: boundary_residual(r, u)
        lower_floor = 0.0
    else:
        residual = lambda u: interior_residual(p, q, u)
        lower_floor = 1e-300

    for _ in range(70):
        u_lo = arb(str(max(lower_floor, u0 - width)))
        u_hi = arb(str(u0 + width))
        if residual(u_lo).upper() <= 0 and residual(u_hi).lower() >= 0:
            return u_lo, u_hi
        width *= 2.0

    stop("could not certify a scalar root bracket")


def G_bounds(p: arb, q: arb) -> Tuple[arb, arb]:
    """Return certified lower and upper Arb bounds for G(p,q)."""
    if p == 0 and q == 0:
        return arb(1), arb(1)

    u_lo, u_hi = certified_u_bracket(p, q)

    if p == 0 or q == 0:
        r = q if p == 0 else p
        sr = r.sqrt()
        return r + sr * u_lo, r + sr * u_hi

    spq = (p * q).sqrt()
    return p + q + 2 * spq * u_lo.cosh(), p + q + 2 * spq * u_hi.cosh()


def cheap_G_lower_bound(p: arb, q: arb) -> arb:
    """
    Fast rigorous lower bound for G used to screen beta_minus lower checks.

    The bound G(p,q) >= (sqrt(p)+sqrt(q))^2 is Lemma G-boundary-lower in the
    manuscript; when one argument is zero we use the stronger boundary version
    G(0,q) >= (sqrt(q)+1/2)^2, and symmetrically.
    """
    if p == 0 and q == 0:
        return arb(1)
    if p == 0:
        return (q.sqrt() + arb(1) / 2) ** 2
    if q == 0:
        return (p.sqrt() + arb(1) / 2) ** 2
    return (p.sqrt() + q.sqrt()) ** 2


# ---------------------------------------------------------------------------
# Loading intervals
# ---------------------------------------------------------------------------


@dataclass
class IntervalInput:
    K: int
    intervals: Dict[str, List[Tuple[arb, arb]]]
    metadata: dict


def parse_interval_list(raw_intervals: Sequence[dict], quantity: str, K: int) -> List[Tuple[arb, arb]]:
    if K != len(raw_intervals) - 1:
        stop(f"K does not match the number of intervals for {quantity}")

    intervals = [(ball(x["lo"]), ball(x["hi"])) for x in raw_intervals]
    if intervals[0] != (arb(0), arb(0)):
        stop(f"the interval for {quantity}[0] must be [0,0]")
    if K >= 1 and (intervals[1][0] > 1 or intervals[1][1] < 1):
        stop(f"the interval for {quantity}[1] must contain 1")

    for k, (L, U) in enumerate(intervals):
        if not proved_nonnegative(U - L):
            stop(f"empty interval for {quantity} at k={k}")
    return intervals


def load_interval_input(filename: str, max_k: Optional[int] = None) -> IntervalInput:
    data = json.loads(Path(filename).read_text())

    # New manuscript-aligned schema:
    #   {"K": K, "intervals": {"gamma": [...], "beta_plus": [...], ...}}
    if isinstance(data.get("intervals"), dict):
        K_file = int(data["K"])
        K = min(K_file, max_k) if max_k is not None else K_file
        intervals = {}
        for quantity, raw in data["intervals"].items():
            if quantity not in QUANTITIES:
                stop(f"unknown quantity {quantity!r}")
            intervals[quantity] = parse_interval_list(raw[: K + 1], quantity, K)
        return IntervalInput(K=K, intervals=intervals, metadata=data.get("metadata", {}))

    # Backward-compatible old schema: one quantity per file.
    quantity = data.get("quantity")
    if quantity not in QUANTITIES:
        stop("quantity must be one of gamma, beta_plus, beta_minus")
    K_file = int(data["K"])
    K = min(K_file, max_k) if max_k is not None else K_file
    return IntervalInput(
        K=K,
        intervals={quantity: parse_interval_list(data["intervals"][: K + 1], quantity, K)},
        metadata=data.get("metadata", {}),
    )


# ---------------------------------------------------------------------------
# Recurrence checks
# ---------------------------------------------------------------------------


def symmetric_indices(k: int) -> Iterable[Tuple[int, int]]:
    """
    Yield (i, weight) for the symmetric sum over i=1,...,k.

    Since G(p,q)=G(q,p), terms i and k+1-i are equal.  The weight is 2 for a
    paired term and 1 for the middle term when k is odd.
    """
    for i in range(1, k // 2 + 1):
        yield i, 2
    if k % 2 == 1:
        yield k // 2 + 1, 1


def update_minimum(old: Optional[arb], new: arb) -> arb:
    """Keep the Arb ball whose lower endpoint is smallest."""
    if old is None or new.lower() < old.lower():
        return new
    return old


def check_gamma_step(k: int, intervals: List[Tuple[arb, arb]]) -> Tuple[arb, arb, dict]:
    lower_sum = arb(0)
    upper_sum = arb(0)

    for i, weight in symmetric_indices(k):
        L_left, L_right = intervals[i - 1][0], intervals[k - i][0]
        U_left, U_right = intervals[i - 1][1], intervals[k - i][1]
        g_lower, _ = G_bounds(L_left, L_right)
        _, g_upper = G_bounds(U_left, U_right)
        lower_sum += weight * g_lower
        upper_sum += weight * g_upper

    lower_rhs = lower_sum / k
    upper_rhs = upper_sum / k
    return lower_rhs - intervals[k][0], intervals[k][1] - upper_rhs, {}


def check_beta_plus_step(k: int, intervals: List[Tuple[arb, arb]]) -> Tuple[arb, arb, dict]:
    lower_best = None
    upper_best = None

    for i, _weight in symmetric_indices(k):
        L_left, L_right = intervals[i - 1][0], intervals[k - i][0]
        U_left, U_right = intervals[i - 1][1], intervals[k - i][1]
        g_lower, _ = G_bounds(L_left, L_right)
        _, g_upper = G_bounds(U_left, U_right)

        # For the lower bound, it is enough to find one large candidate.
        if lower_best is None or g_lower.lower() > lower_best.lower():
            lower_best = g_lower

        # For the upper bound, we must dominate every candidate.
        if upper_best is None or g_upper.upper() > upper_best.upper():
            upper_best = g_upper

    return lower_best - intervals[k][0], intervals[k][1] - upper_best, {}


def beta_minus_preferred_indices(k: int) -> List[int]:
    """Likely minimising splits, used only to prove an upper inclusion quickly."""
    if k == 2:
        return [1]
    candidates = [2, 1, k // 2 + 1]
    return [i for i in candidates if 1 <= i <= k // 2 + 1]


def check_beta_minus_step(k: int, intervals: List[Tuple[arb, arb]]) -> Tuple[arb, arb, dict]:
    L_k, U_k = intervals[k]
    min_lower_slack = None
    lower_exact_terms = 0
    lower_screened_terms = 0

    # Lower recurrence inclusion: L_k <= min_i G(L_{i-1}, L_{k-i}).
    # Hence every split must be lower-bounded by L_k.  Most splits are certified
    # by a cheap analytic lower bound; the few remaining near-minimizers are
    # certified by the scalar G interval certificate.
    for i, _weight in symmetric_indices(k):
        L_left, L_right = intervals[i - 1][0], intervals[k - i][0]
        cheap_slack = cheap_G_lower_bound(L_left, L_right) - L_k
        if proved_nonnegative(cheap_slack):
            lower_screened_terms += 1
            min_lower_slack = update_minimum(min_lower_slack, cheap_slack)
            continue

        g_lower, _ = G_bounds(L_left, L_right)
        lower_exact_terms += 1
        exact_slack = g_lower - L_k
        min_lower_slack = update_minimum(min_lower_slack, exact_slack)

    # Upper recurrence inclusion: min_i G(U_{i-1}, U_{k-i}) <= U_k.
    # A single split with G upper endpoint <= U_k is sufficient.  The split i=2
    # is the intended candidate for k>=3; if it ever fails, scan all splits.
    upper_slack = None
    upper_exact_terms = 0
    tried = set()
    for i in beta_minus_preferred_indices(k):
        tried.add(i)
        U_left, U_right = intervals[i - 1][1], intervals[k - i][1]
        _, g_upper = G_bounds(U_left, U_right)
        upper_exact_terms += 1
        upper_slack = U_k - g_upper
        if proved_nonnegative(upper_slack):
            return min_lower_slack, upper_slack, {
                "lower_screened_terms": lower_screened_terms,
                "lower_exact_terms": lower_exact_terms,
                "upper_exact_terms": upper_exact_terms,
            }

    # Robust fallback: search the full symmetric list for any upper certificate.
    for i, _weight in symmetric_indices(k):
        if i in tried:
            continue
        U_left, U_right = intervals[i - 1][1], intervals[k - i][1]
        _, g_upper = G_bounds(U_left, U_right)
        upper_exact_terms += 1
        candidate_slack = U_k - g_upper
        if upper_slack is None or candidate_slack.lower() > upper_slack.lower():
            upper_slack = candidate_slack
        if proved_nonnegative(candidate_slack):
            return min_lower_slack, candidate_slack, {
                "lower_screened_terms": lower_screened_terms,
                "lower_exact_terms": lower_exact_terms,
                "upper_exact_terms": upper_exact_terms,
            }

    return min_lower_slack, upper_slack, {
        "lower_screened_terms": lower_screened_terms,
        "lower_exact_terms": lower_exact_terms,
        "upper_exact_terms": upper_exact_terms,
    }


def check_one_step(quantity: str, k: int, intervals: List[Tuple[arb, arb]]) -> Tuple[arb, arb, dict]:
    if quantity == "gamma":
        return check_gamma_step(k, intervals)
    if quantity == "beta_plus":
        return check_beta_plus_step(k, intervals)
    if quantity == "beta_minus":
        return check_beta_minus_step(k, intervals)
    stop(f"unknown quantity {quantity}")


@dataclass
class CheckResult:
    min_lower_slack: object
    min_upper_slack: object
    counters: dict


def merge_counters(total: dict, update: dict) -> None:
    for key, value in update.items():
        total[key] = total.get(key, 0) + value


def check_recurrence_serial(quantity: str, K: int, intervals: List[Tuple[arb, arb]]) -> CheckResult:
    min_lower_slack = None
    min_upper_slack = None
    counters = {}

    for k in range(2, K + 1):
        lower_slack, upper_slack, step_counters = check_one_step(quantity, k, intervals)

        if not proved_nonnegative(lower_slack):
            stop(f"{quantity}: lower recurrence inequality failed at k={k}; slack={lower_slack}")
        if not proved_nonnegative(upper_slack):
            stop(f"{quantity}: upper recurrence inequality failed at k={k}; slack={upper_slack}")

        min_lower_slack = update_minimum(min_lower_slack, lower_slack)
        min_upper_slack = update_minimum(min_upper_slack, upper_slack)
        merge_counters(counters, step_counters)

        if k in {10, 20, 50, 100, 500, 1000, 2000, 5000}:
            print(f"{quantity}: checked recurrence inequalities for every 2 <= k <= {k}", flush=True)

    return CheckResult(min_lower_slack, min_upper_slack, counters)


# The next functions are for optional parallel checking.

_WORK_QUANTITY = None
_WORK_INTERVALS = None


def init_worker(filename: str, quantity: str, precision: int, max_k: Optional[int]):
    global _WORK_QUANTITY, _WORK_INTERVALS
    ctx.prec = precision
    loaded = load_interval_input(filename, max_k=max_k)
    _WORK_QUANTITY = quantity
    _WORK_INTERVALS = loaded.intervals[quantity]


def check_block(block: Tuple[int, int]) -> dict:
    start, end = block
    min_lower_slack = None
    min_upper_slack = None
    counters = {}

    for k in range(start, end + 1):
        lower_slack, upper_slack, step_counters = check_one_step(_WORK_QUANTITY, k, _WORK_INTERVALS)

        if not proved_nonnegative(lower_slack):
            return {"ok": False, "message": f"{_WORK_QUANTITY}: lower recurrence inequality failed at k={k}"}
        if not proved_nonnegative(upper_slack):
            return {"ok": False, "message": f"{_WORK_QUANTITY}: upper recurrence inequality failed at k={k}"}

        min_lower_slack = update_minimum(min_lower_slack, lower_slack)
        min_upper_slack = update_minimum(min_upper_slack, upper_slack)
        merge_counters(counters, step_counters)

    return {
        "ok": True,
        "start": start,
        "end": end,
        "lower_slack": min_lower_slack.str(40),
        "upper_slack": min_upper_slack.str(40),
        "lower_endpoint": float(min_lower_slack.lower()),
        "upper_endpoint": float(min_upper_slack.lower()),
        "counters": counters,
    }


def make_blocks(K: int, block_size: int) -> List[Tuple[int, int]]:
    blocks = []
    start = 2
    while start <= K:
        end = min(K, start + block_size - 1)
        blocks.append((start, end))
        start = end + 1
    return blocks


def check_recurrence_parallel(filename: str, quantity: str, K: int, precision: int, jobs: int, block_size: int, max_k: Optional[int]) -> CheckResult:
    blocks = make_blocks(K, block_size)
    min_lower = None
    min_upper = None
    min_lower_endpoint = None
    min_upper_endpoint = None
    counters = {}

    with mp.Pool(processes=jobs, initializer=init_worker, initargs=(filename, quantity, precision, max_k)) as pool:
        for result in pool.imap_unordered(check_block, blocks):
            if not result["ok"]:
                stop(result["message"])

            if min_lower_endpoint is None or result["lower_endpoint"] < min_lower_endpoint:
                min_lower_endpoint = result["lower_endpoint"]
                min_lower = result["lower_slack"]
            if min_upper_endpoint is None or result["upper_endpoint"] < min_upper_endpoint:
                min_upper_endpoint = result["upper_endpoint"]
                min_upper = result["upper_slack"]
            merge_counters(counters, result.get("counters", {}))

            print(f"{quantity}: checked recurrence inequalities for every {result['start']} <= k <= {result['end']}", flush=True)

    return CheckResult(min_lower, min_upper, counters)


# ---------------------------------------------------------------------------
# Bounds for the limiting constants
# ---------------------------------------------------------------------------


def constant_bounds(quantity: str, K: int, intervals: List[Tuple[arb, arb]]) -> Tuple[str, arb, arb, dict]:
    if quantity == "gamma":
        total = arb(0)
        for j in range(K + 1):
            total += intervals[j][0].sqrt()
        lower = (2 * total / ((K + 1) * (K + 2))) ** 2
        upper = intervals[K][1] / (K * K)
        return "c_gamma", lower, upper, {"manuscript_use": "upper bounds c_typ"}

    if quantity == "beta_plus":
        lower = intervals[K][0] / ((K + 1) * (K + 1))
        upper = intervals[K][1] / (K * K)
        return "c_plus", lower, upper, {}

    if quantity == "beta_minus":
        lower = arb(1) / 4
        upper = intervals[K][1] / (K * K)
        return "c_minus", lower, upper, {"lower_bound_source": "analytic bound c_minus >= 1/4 from the manuscript"}

    stop(f"unknown quantity {quantity}")


def goal_pair(quantity: str, args) -> Tuple[Optional[str], Optional[str]]:
    if args.manuscript_goals:
        goals = MANUSCRIPT_GOALS[quantity]
        return goals["lower"], goals["upper"]
    return args.lower_goal, args.upper_goal


def check_goals(quantity: str, lower_bound: arb, upper_bound: arb, args) -> dict:
    lower_goal, upper_goal = goal_pair(quantity, args)
    out = {
        "lower_goal": lower_goal or "not checked",
        "lower_goal_slack": "not checked",
        "upper_goal": upper_goal or "not checked",
        "upper_goal_slack": "not checked",
    }

    if lower_goal is not None:
        lower_goal_slack = lower_bound - ball(lower_goal)
        if not proved_positive(lower_goal_slack):
            stop(f"{quantity}: could not prove lower bound > {lower_goal}")
        out["lower_goal_slack"] = lower_goal_slack
        print(f"{quantity}: certified lower bound is above {lower_goal}; slack={lower_goal_slack}")

    if upper_goal is not None:
        upper_goal_slack = ball(upper_goal) - upper_bound
        if not proved_positive(upper_goal_slack):
            stop(f"{quantity}: could not prove upper bound < {upper_goal}")
        out["upper_goal_slack"] = upper_goal_slack
        print(f"{quantity}: certified upper bound is below {upper_goal}; slack={upper_goal_slack}")

    return out


def ctyp_lower_bound_from_beta_minus(
    intervals: List[Tuple[arb, arb]], cutoff: int = CTYP_LOWER_CUTOFF
) -> arb:
    """Certified finite lower bound on c_typ from beta_j^- for 1 <= j <= cutoff."""
    if cutoff >= len(intervals):
        stop(f"the c_typ lower bound requires beta_minus intervals through k={cutoff}")

    total = arb(0)
    for j in range(1, cutoff + 1):
        total += intervals[j][0].sqrt()
    return (2 * total / ((cutoff + 1) * (cutoff + 2))) ** 2


def check_ctyp_lower_goal(intervals: List[Tuple[arb, arb]], args) -> dict:
    """Compute the beta_minus-derived lower bound and optionally check its manuscript goal."""
    lower_bound = ctyp_lower_bound_from_beta_minus(intervals)
    out = {
        "passed": True,
        "cutoff": CTYP_LOWER_CUTOFF,
        "lower_bound": lower_bound,
        "lower_goal": CTYP_LOWER_GOAL if args.manuscript_goals else "not checked",
        "lower_goal_slack": "not checked",
        "source": "certified lower endpoints for beta_minus and the finite lower bound proved in the manuscript",
    }

    if args.manuscript_goals:
        slack = lower_bound - ball(CTYP_LOWER_GOAL)
        if not proved_positive(slack):
            stop(f"c_typ: could not prove lower bound > {CTYP_LOWER_GOAL}")
        out["lower_goal_slack"] = slack
        print(f"c_typ: certified lower bound is above {CTYP_LOWER_GOAL}; slack={slack}")

    return out


def write_slacks(filename: str, summary: dict) -> None:
    def convert(obj):
        if hasattr(obj, "str"):
            return obj.str(40)
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    Path(filename).write_text(json.dumps(convert(summary), indent=2) + "\n")


# ---------------------------------------------------------------------------
# Command line interface
# ---------------------------------------------------------------------------


def parse_only(value: Optional[str], available: Iterable[str]) -> List[str]:
    available = list(available)
    if value is None or value == "all":
        return available
    requested = [x.strip() for x in value.split(",") if x.strip()]
    for quantity in requested:
        if quantity not in available:
            stop(f"requested quantity {quantity!r} is not present in the interval file")
    return requested


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("intervals_json", help="input interval file; normally intervals.json")
    parser.add_argument("--prec", type=int, default=128, help="Arb precision in bits")
    parser.add_argument("--jobs", type=int, default=1, help="number of parallel worker processes")
    parser.add_argument("--block-size", type=int, default=25, help="number of k-values per parallel block")
    parser.add_argument("--max-k", type=int, default=None, help="check only intervals through this k")
    parser.add_argument("--only", default="all", help="comma-separated subset: gamma,beta_plus,beta_minus; default all")
    parser.add_argument("--lower-goal", default=None, help="optional lower goal; intended when checking one quantity")
    parser.add_argument("--upper-goal", default=None, help="optional upper goal; intended when checking one quantity")
    parser.add_argument("--manuscript-goals", action="store_true", help="check the numerical goals stated in the manuscript")
    parser.add_argument("--slacks", default="slacks.json", help="output JSON file; default slacks.json")
    args = parser.parse_args()

    ctx.prec = args.prec
    loaded = load_interval_input(args.intervals_json, max_k=args.max_k)
    quantities = parse_only(args.only, loaded.intervals.keys())

    if (args.lower_goal is not None or args.upper_goal is not None) and len(quantities) != 1 and not args.manuscript_goals:
        stop("--lower-goal/--upper-goal without --manuscript-goals is only allowed when checking one quantity")

    print(f"loaded {args.intervals_json}; K={loaded.K}; Arb precision={ctx.prec} bits")
    print("checking quantities:", ", ".join(quantities))

    summary = {
        "input_file": str(args.intervals_json),
        "K": loaded.K,
        "arb_precision_bits": ctx.prec,
        "passed": True,
        "quantities": {},
        "derived_bounds": {},
    }

    for quantity in quantities:
        intervals = loaded.intervals[quantity]
        print(f"\n{quantity}: starting recurrence certification")
        if args.jobs <= 1:
            result = check_recurrence_serial(quantity, loaded.K, intervals)
        else:
            print(f"{quantity}: using {args.jobs} worker processes", flush=True)
            result = check_recurrence_parallel(
                args.intervals_json, quantity, loaded.K, args.prec, args.jobs, args.block_size, args.max_k
            )

        constant_name, lower_bound, upper_bound, notes = constant_bounds(quantity, loaded.K, intervals)
        goal_data = check_goals(quantity, lower_bound, upper_bound, args)

        print(f"{quantity}: PASS")
        print(f"{quantity}: minimum lower-recurrence slack: {result.min_lower_slack}")
        print(f"{quantity}: minimum upper-recurrence slack: {result.min_upper_slack}")
        print(f"{quantity}: certified {constant_name} lower bound: {lower_bound}")
        print(f"{quantity}: certified {constant_name} upper bound: {upper_bound}")

        summary["quantities"][quantity] = {
            "passed": True,
            "min_lower_recurrence_slack": result.min_lower_slack,
            "min_upper_recurrence_slack": result.min_upper_slack,
            "counters": result.counters,
            f"{constant_name}_lower": lower_bound,
            f"{constant_name}_upper": upper_bound,
            **goal_data,
            **notes,
        }

        if quantity == "beta_minus" and loaded.K >= CTYP_LOWER_CUTOFF:
            ctyp_data = check_ctyp_lower_goal(intervals, args)
            summary["derived_bounds"]["c_typ_lower_from_beta_minus"] = ctyp_data
            print(f"beta_minus: certified c_typ lower bound: {ctyp_data['lower_bound']}")

    if args.manuscript_goals and "beta_minus" in quantities and loaded.K < CTYP_LOWER_CUTOFF:
        stop(f"the manuscript c_typ lower goal requires --max-k >= {CTYP_LOWER_CUTOFF}")

    write_slacks(args.slacks, summary)
    print(f"\nwrote {args.slacks}")


if __name__ == "__main__":
    main()
