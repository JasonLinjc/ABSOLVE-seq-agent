"""Write the English run report beside its supporting result tables."""
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

from verify_results import BASELINES

STRATEGIES = ("Unfiltered", "Dedup", "Dedud", "Dedud+Dedup")


def read_table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_coverage(rows, panel=False):
    keys = {(r["target"], r["strategy"], r["editor"] if panel else "") for r in rows}
    targets = {r["target"] for r in rows}
    editors = ("WT", "HiFi") if panel else ("",)
    expected = {(t, s, e) for t in targets for s in STRATEGIES for e in editors}
    if len(targets) != 79 or keys != expected or len(keys) != len(rows):
        raise ValueError("Incomplete or duplicated target comparisons in report input")


def number(value):
    return "undefined" if value is None or value == "" else f"{float(value):.6g}"


def rate(k, n):
    return number(100 * k / n if n else None)


def paired_report(folder, mode):
    rows = read_table(folder / "paired_target_rates.tsv")
    check_coverage(rows)
    settings = {(r["donor"], r["replicate"], r.get("maxi_method", "unrecorded"), r.get("target_filter", "unrecorded")) for r in rows}
    if len(settings) != 1:
        raise ValueError("Inconsistent pair or filter settings")
    donor, rep, maxi, target_filter = settings.pop()
    lines = [f"# ABSOLVE-seq paired analysis: {donor} / {rep}", "",
             f"Read mode: {mode}. Plasmid barcode rule: {maxi}. Target filter: {target_filter}.",
             "The table contains 316 unique comparisons across 79 targets and four strategies.",
             "Differences are HiFi minus NoEP in percentage points (pp)."]
    if (folder / "qc.json").is_file():
        qc = json.loads((folder / "qc.json").read_text())
        lines += ["", "## Read retention", "", "| Arm | Library | Raw pairs | Prepared pairs | Retained (%) |",
                  "| --- | --- | ---: | ---: | ---: |"]
        for arm in ("NoEP", "HiFi"):
            q = qc[arm]
            lines.append(f"| {arm} | {q['sample_id']} | {q['pairs']:,} | {q['prepared_pairs']:,} | "
                         f"{rate(q['prepared_pairs'], q['pairs'])} |")
    lines += ["", "## Results by strategy", "",
              "Arm rates pool the counts across targets. The median uses targets with a defined paired difference.",
              "A zero denominator gives an undefined rate.", "",
              "| Strategy | NoEP pooled (%) | HiFi pooled (%) | Median target Δ (pp) | Undefined Δ |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for strategy in STRATEGIES:
        selected = [r for r in rows if r["strategy"] == strategy]
        pooled = [rate(sum(int(r[arm+"_indel"]) for r in selected),
                       sum(int(r[arm+"_denominator"]) for r in selected)) for arm in ("NoEP", "HiFi")]
        deltas = [float(r["delta_HiFi_minus_NoEP_pp"]) for r in selected if r["delta_HiFi_minus_NoEP_pp"]]
        lines.append(f"| {strategy} | {pooled[0]} | {pooled[1]} | {number(median(deltas) if deltas else None)} "
                     f"| {len(selected) - len(deltas)} |")
    lines += ["", "## Manuscript focal target: 1617_OT_0040_ALT", "",
              "| Strategy | NoEP edited / total | NoEP (%) | HiFi edited / total | HiFi (%) | Δ (pp) |",
              "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in rows:
        if r["target"] == "1617_OT_0040_ALT":
            lines.append(f"| {r['strategy']} | {r['NoEP_indel']} / {r['NoEP_denominator']} | {number(r['NoEP_rate_pct'])} "
                         f"| {r['HiFi_indel']} / {r['HiFi_denominator']} | {number(r['HiFi_rate_pct'])} "
                         f"| {number(r['delta_HiFi_minus_NoEP_pp'])} |")
    return lines, ["paired_target_rates.tsv", "strategy_summary.tsv", "qc.json", "library_counts.tsv",
                   "sample_inventory.tsv", "provenance.json"]


def panel_report(folder, mode):
    rows = read_table(folder / "panel_comparison.tsv")
    check_coverage(rows, panel=True)
    differences = [r for r in rows if r["agent_hit"] != r["submitted_hit"]]
    max_delta = max(abs(float(r["delta_pp"]) - float(r["manuscript_delta_pp"])) for r in rows)
    max_q = max(abs(float(r["q_bh"]) - float(r["manuscript_q_bh"])) for r in rows)
    lines = ["# ABSOLVE-seq panel analysis", "", f"Analysis: {mode}.",
             f"Primary calls agree with the submitted supplements for {len(rows)-len(differences)} of 632 comparisons.",
             "The primary call requires edited minus NoEP Δ > 0.1 pp and within-editor, within-strategy BH q < 0.05.",
             f"The largest absolute differences from the submitted estimates are {number(max_delta)} pp for Δ and {number(max_q)} for q.",
             "", "| Strategy | Editor | Agent positives | Submitted positives | Discordant calls |",
             "| --- | --- | ---: | ---: | ---: |"]
    for strategy in STRATEGIES:
        for editor in ("WT", "HiFi"):
            subset = [r for r in rows if r["strategy"] == strategy and r["editor"] == editor]
            if any(r[k] not in {"TRUE", "FALSE"} for r in subset for k in ("agent_hit", "submitted_hit")):
                raise ValueError("Invalid panel call")
            lines.append(f"| {strategy} | {editor} | {sum(r['agent_hit']=='TRUE' for r in subset)} "
                         f"| {sum(r['submitted_hit']=='TRUE' for r in subset)} "
                         f"| {sum(r['agent_hit']!=r['submitted_hit'] for r in subset)} |")
    if differences:
        lines += ["", "## Discordant calls", "",
                  "| Target | Strategy | Editor | Agent Δ (pp) | Submitted Δ (pp) | Agent q | Submitted q |",
                  "| --- | --- | --- | ---: | ---: | ---: | ---: |"]
        for r in differences:
            lines.append(f"| {r['target']} | {r['strategy']} | {r['editor']} | {number(r['delta_pp'])} "
                         f"| {number(r['manuscript_delta_pp'])} | {number(r['q_bh'])} | {number(r['manuscript_q_bh'])} |")
    return lines, ["panel_comparison.tsv", "comparison_summary.tsv", "library_counts.tsv", "model_field_comparison.tsv",
                   "donor_weights.tsv", "sample_inventory.tsv", "provenance.json"]


def write_report(folder: Path, mode, *, panel=False, run_config=None):
    lines, files = panel_report(folder, mode) if panel else paired_report(folder, mode)
    profiles = ("panel-refit",) if panel else ("pair-fast", "pair-historical")
    for profile in profiles:
        reference = BASELINES[profile]
        if all((folder / name).is_file() and hashlib.sha256((folder / name).read_bytes()).hexdigest() == digest
               for name, digest in reference.items()):
            lines += ["", f"Benchmark: all reference tables in `{profile}` match the accepted SHA-256 digests."]
            break
    if (folder / "power_provenance.json").is_file():
        power = read_table(folder / "power_estimates.tsv")
        receipt = json.loads((folder / "power_provenance.json").read_text())
        stem = "Figure3b" if panel else "paired_effect_power"
        book = "Figure3b_source_data" if panel else stem
        lines += ["", "## Estimated detection power", "",
                  f"Each simulated comparison uses the observed denominators and {receipt['baseline']} background. "
                  f"Power is the fraction of {receipt['simulations']:,} binomial simulations passing a one-sided "
                  f"Wald test at nominal α = {receipt['alpha']}. The assumed increases are 0.01 and 0.1 percentage points. "
                  "This measures conditional detection power before multiple-testing correction. "
                  "The TSV gives a 95% Monte Carlo interval and a reason for any unavailable estimate.",
                  "", "| Strategy | Editor | Assumed increase (pp) | Median power | Targets with estimates |",
                  "| --- | --- | ---: | ---: | ---: |"]
        for strategy in STRATEGIES:
            for editor in (("HiFi", "WT") if panel else ("HiFi",)):
                for effect in receipt["effects_pp"]:
                    values = [float(r["power"]) for r in power if r["strategy"] == strategy and r["editor"] == editor
                              and float(r["effect_pp"]) == effect and r["status"] == "estimated"]
                    med = f"{median(values):.1%}" if values else "undefined"
                    lines.append(f"| {strategy} | {editor} | {effect} | {med} | {len(values)} / 79 |")
        lines += ["", f"![Dedud effects and count-based power](figures/{stem}.png)", "",
                  f"[Figure PDF](figures/{stem}.pdf) · [Excel source data](figures/{book}.xlsx)"]
        files.extend(["power_estimates.tsv", "power_provenance.json", "figures/plotted_data.tsv"])
    lines += ["", "## Supporting files", ""]
    lines.extend(f"- [{name}]({name})" for name in files if (folder / name).is_file())
    if run_config:
        lines.append(f"- [Run configuration](../{run_config})")
    (folder / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def complete_report(folder, mode, *, panel=False, run_config=None, image=None):
    """Generate power, figure and source data inside the caller's staging directory."""
    import subprocess
    import sys
    script = Path(__file__).with_name("power_analysis.py").resolve()
    command = [sys.executable, str(script), str(folder.resolve())]
    if panel:
        command.append("--panel")
    if image:
        binds = ",".join(str(p) for p in (script.parent, folder.resolve()))
        command = ["apptainer", "exec", "--cleanenv", "--bind", binds, image, "python", *command[1:]]
    subprocess.run(command, check=True)
    write_report(folder, mode, panel=panel, run_config=run_config)
