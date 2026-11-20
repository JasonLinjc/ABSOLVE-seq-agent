"""Conditional binomial power for the ABSOLVE-seq paired and manuscript designs."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

import numpy as np

from analysis_report import check_coverage, read_table

EFFECTS = (0.01, 0.1)
DRAWS = 1000
SEED = 20261121
ALPHA = 0.05
CUTOFF = NormalDist().inv_cdf(1 - ALPHA)
ARMS = ("NoEP", "HiFi", "WT")


def pair_wald(control, treated, depths):
    """Closed-form Firth fit for a saturated two-group binomial model."""
    p0 = (control + 0.5) / (depths[0] + 1)
    p1 = (treated + 0.5) / (depths[1] + 1)
    log_odds_ratio = np.log(p1 / (1 - p1)) - np.log(p0 / (1 - p0))
    variance = 1 / (depths[0] * p0 * (1 - p0)) + 1 / (depths[1] * p1 * (1 - p1))
    return log_odds_ratio / np.sqrt(variance)


def panel_wald(outcomes, depths, design):
    """Mean bias-reduced logistic score and model-based Wald z, one fit per draw."""
    proportions = (outcomes + 0.5) / (depths + 1)
    beta = np.log(proportions / (1 - proportions)) @ np.linalg.pinv(design).T
    for _ in range(150):
        eta = beta @ design.T
        prob = np.exp(-np.logaddexp(0, -eta))
        weights = depths * prob * (1 - prob)
        info = design.T @ (weights[:, :, None] * design)
        covariance = np.linalg.inv(info)
        leverage = weights * np.sum((design @ covariance) * design, axis=2)
        score = (outcomes - depths * prob + leverage * (0.5 - prob)) @ design
        step = (covariance @ score[:, :, None])[:, :, 0]
        largest = np.max(np.abs(step), axis=1, keepdims=True)
        if np.all(largest < 1e-7):
            se = np.sqrt(np.diagonal(covariance, axis1=1, axis2=2))
            return beta[:, 1:3] / se[:, 1:3]
        beta += step / np.maximum(largest, 1)
    raise ValueError("A simulated donor-adjusted Firth fit did not converge")


def estimate(target, strategy, depths, baseline, editors, design, draws, seed):
    rows = []
    for effect in EFFECTS:
        key = f"{seed}/{target}/{strategy}/{effect}"
        rng = np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(key.encode()).digest()[:16], "big")))
        probabilities = baseline + (np.arange(len(depths)) % (1 + len(editors)) != 0) * effect / 100
        status = "estimated"
        hits = [None] * len(editors)
        if np.any(depths == 0):
            status = "zero denominator"
        elif np.any(probabilities > 1):
            status = "effect exceeds probability range"
        else:
            outcomes = rng.binomial(depths, probabilities, size=(draws, len(depths)))
            if design is None:
                frequency = np.ones(draws, dtype=int)
                scores = pair_wald(outcomes[:, 0], outcomes[:, 1], depths)[:, None]
            else:
                unique, frequency = np.unique(outcomes, axis=0, return_counts=True)
                scores = panel_wald(unique, depths, design)
            hits = (scores > CUTOFF).T @ frequency
        for editor, hit in zip(editors, hits):
            power = None if hit is None else int(hit) / draws
            if power is None:
                low = high = None
            else:
                z = NormalDist().inv_cdf(0.975)
                center = power + z**2 / (2 * draws)
                half = z * np.sqrt(power * (1 - power) / draws + z**2 / (4 * draws**2))
                scale = 1 + z**2 / draws
                # At 0/n and n/n the Wilson endpoint is exact. Avoid cancellation at the boundary.
                low = 0.0 if hit == 0 else (center - half) / scale
                high = 1.0 if hit == draws else (center + half) / scale
            rows.append(dict(target=target, strategy=strategy, editor=editor, effect_pp=effect,
                             power=power, mc_low95=low, mc_high95=high, detections=None if hit is None else int(hit),
                             simulations=0 if hit is None else draws, status=status))
    return rows


def count(row, numerator, denominator):
    k, n = int(row[numerator]), int(row[denominator])
    if not 0 <= k <= n:
        raise ValueError("Invalid binomial counts")
    return k, n


def estimate_panel(folder, draws, seed):
    comparisons = read_table(folder / "panel_comparison.tsv")
    check_coverage(comparisons, panel=True)
    cells = defaultdict(dict)
    counts = read_table(folder / "library_counts.tsv")
    for row in counts:
        arm = row["arm"]
        key = (row["donor"], row["replicate"], arm)
        group = cells[row["target"], row["strategy"]]
        if key in group or arm not in ARMS:
            raise ValueError("Duplicate or unknown library in panel counts")
        group[key] = count(row, "indel", "denominator")
    expected = {(r["target"], r["strategy"]) for r in comparisons}
    if set(cells) != expected:
        raise ValueError("Target or strategy mismatch between panel and counts")
    fields = read_table(folder / "model_field_comparison.tsv")
    fitted = {(r["target"], r["strategy"], r["field"]): float(r["agent"])
              for r in fields if r["field"].startswith("predicted-") and r["field"].endswith("-NoEP")}
    results = []
    libraries = set(next(iter(cells.values())))
    donors = sorted({d for d, _, _ in libraries})
    if len(libraries) != 26 or len(donors) != 3:
        raise ValueError("Expected the 26-library three-donor panel")
    design = np.array([[1, a == "HiFi", a == "WT", d == donors[1], d == donors[2]]
                       for d in donors for a in ARMS], dtype=float)
    for (target, strategy), group in sorted(cells.items()):
        if set(group) != libraries:
            raise ValueError("Inconsistent library design across targets")
        depths = np.array([sum(n for (dn, _, arm), (_, n) in group.items() if dn == d and arm == a)
                           for d in donors for a in ARMS])
        baseline = np.repeat([fitted[target, strategy, f"predicted-{d}-NoEP"] for d in donors], len(ARMS))
        if not np.all((baseline >= 0) & (baseline <= 1)):
            raise ValueError("Invalid donor-specific NoEP estimate")
        results.extend(estimate(target, strategy, depths, baseline, ARMS[1:], design, draws, seed))
    return results


def estimate_pair(folder, draws, seed):
    pairs = read_table(folder / "paired_target_rates.tsv")
    check_coverage(pairs)
    results = []
    for row in sorted(pairs, key=lambda r: (r["target"], r["strategy"])):
        k, nc = count(row, "NoEP_indel", "NoEP_denominator")
        _, nt = count(row, "HiFi_indel", "HiFi_denominator")
        baseline = (k + 0.5) / (nc + 1)
        results.extend(estimate(row["target"], row["strategy"], np.array([nc, nt]),
                                baseline, ("HiFi",), None, draws, seed))
    return results


def calculate(folder, panel, draws=DRAWS, seed=SEED):
    if draws < 1 or seed < 0:
        raise ValueError("Use positive simulations and a nonnegative seed")
    rows = estimate_panel(folder, draws, seed) if panel else estimate_pair(folder, draws, seed)
    with (folder / "power_estimates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    inputs = ["library_counts.tsv", "model_field_comparison.tsv", "panel_comparison.tsv"] if panel else ["paired_target_rates.tsv"]
    receipt = dict(simulations=draws, seed=seed, effects_pp=EFFECTS, alpha=ALPHA,
                   test="one-sided Firth Wald, nominal alpha", numpy=np.__version__,
                   generator="PCG64, SHA-256 seed per target/strategy/effect",
                   interval="95% Wilson interval for Monte Carlo detection probability",
                   method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   design="treatment + donor" if panel else "two-group conditional binomial",
                   baseline="fitted donor-specific NoEP" if panel else "Firth-smoothed NoEP fraction",
                   input_sha256={name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in inputs})
    (folder / "power_provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--panel", action="store_true")
    args = parser.parse_args()
    from result_figures import render
    from verify_results import verify_power
    if any((args.folder / name).exists() for name in ("power_estimates.tsv", "power_provenance.json", "figures")):
        parser.error("Power or figure outputs already exist. Use a new results directory.")
    render(args.folder, calculate(args.folder, args.panel), args.panel)
    verify_power(args.folder)
