"""Render a pair summary or manuscript Figure 3b with newly simulated power."""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from openpyxl import Workbook

from analysis_report import read_table

FIGURE_SKILL = Path(__file__).resolve().parent / ".agents/skills/absolve-seq-figures/scripts"
sys.path.insert(0, str(FIGURE_SKILL))
import manuscript_figures as manuscript


def panel_figure(folder, power):
    panel = manuscript.unique(read_table(folder / "panel_comparison.tsv"), ["target", "strategy", "editor"])
    fields = manuscript.unique(read_table(folder / "model_field_comparison.tsv"), ["target", "strategy", "field"])
    order = read_table(manuscript.ORDER)
    rows = []
    for editor, suffix in manuscript.EDITORS.items():
        for i, label in enumerate(order, 1):
            t = label["target"]
            comp = panel[t, "Dedud", editor]
            control = float(comp["control_rate_pct"])
            rows.append(dict(target=t, label=label["label"], order=i, guide=t[:4], strategy="Dedud",
                             editor=editor, delta_pp=float(comp["delta_pp"]),
                             lower_pp=100 * float(fields[t, "Dedud", "confLow-" + suffix]["agent"]) - control,
                             upper_pp=100 * float(fields[t, "Dedud", "confHigh-" + suffix]["agent"]) - control,
                             q_bh=float(comp["q_bh"]), power001=power[t, editor, 0.01],
                             power01=power[t, editor, 0.1], positive_control=t in manuscript.CONTROLS))
    sources = {name: manuscript.digest(folder / name) for name in
               ("panel_comparison.tsv", "model_field_comparison.tsv", "power_estimates.tsv", "power_provenance.json")}
    sources[manuscript.ORDER.name] = manuscript.digest(manuscript.ORDER)
    manuscript.save_figure(folder / "figures", rows, sources, "agent refit", "simulated from this panel's counts and fitted NoEP rates")


def pair_figure(folder, power):
    rows = sorted((r for r in read_table(folder / "paired_target_rates.tsv") if r["strategy"] == "Dedud"),
                  key=lambda r: r["target"])
    source = []
    for i, r in enumerate(rows, 1):
        delta = r["delta_HiFi_minus_NoEP_pp"]
        source.append(dict(order=i, target=r["target"], strategy="Dedud", editor="HiFi",
                           delta_pp=float(delta) if delta else None,
                           power001=power[r["target"], "HiFi", 0.01], power01=power[r["target"], "HiFi", 0.1]))
    dest = folder / "figures"
    dest.mkdir()
    with (dest / "plotted_data.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, source[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(source)
    book = Workbook()
    sheet = book.active
    sheet.title = "Paired Dedud"
    sheet.append(list(source[0]))
    for row in source:
        sheet.append(list(row.values()))
    sheet.freeze_panes = "C2"
    sheet.column_dimensions["B"].width = 25
    book.save(dest / "paired_effect_power.xlsx")
    x = np.arange(1, len(source) + 1)
    fig, (effect_ax, power_ax) = plt.subplots(2, figsize=(9, 5.6), sharex=True, layout="constrained")
    effect_ax.axhline(0, color="#aab3b8", linewidth=0.8)
    effect_ax.scatter(x, [r["delta_pp"] for r in source], s=15, color="#175c75")
    effect_ax.set_yscale("symlog", linthresh=0.1)
    effect_ax.set_ylabel("HiFi − NoEP (pp)\nSymmetric log scale")
    for column, effect, color, marker in (("power001", "0.01 pp", "#137f78", "o"), ("power01", "0.1 pp", "#b65d23", "x")):
        power_ax.scatter(x, [r[column] for r in source], color=color, s=18, marker=marker, label=effect)
    power_ax.set(ylim=(-0.05, 1.12), ylabel="Detection probability", xlabel="Target index in the source table")
    power_ax.legend(title="Assumed effect", loc="upper right", ncol=2)
    fig.suptitle(f"{rows[0]['donor']} / {rows[0]['replicate']}   |   Dedud paired effects and conditional power")
    fig.savefig(dest / "paired_effect_power.png", dpi=180)
    fig.savefig(dest / "paired_effect_power.pdf")
    plt.close(fig)


def render(folder, rows, panel):
    power = {(r["target"], r["editor"], r["effect_pp"]): r["power"] for r in rows if r["strategy"] == "Dedud"}
    if panel:
        panel_figure(folder, power)
    else:
        pair_figure(folder, power)
