"""Regenerate ABSOLVE-seq Figure 3b and its Excel source data."""
import argparse
import csv
from tempfile import TemporaryDirectory
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter
import numpy as np

ORDER = Path(__file__).resolve().parents[1] / "references/target_order.tsv"
EDITORS = {"WT": "WTCas9", "HiFi": "HiFi"}
COLORS = {"WT": "#00852a", "HiFi": "#e72a98"}
CONTROLS = {"1450_OT_0000_REF", "1617_OT_0000_REF", "1617_OT_0040_ALT"}
POWER = LinearSegmentedColormap.from_list("power", ["#2626ab", "#b3b3b3", "#cc1616"])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def unique(rows, keys):
    indexed = {tuple(row[k] for k in keys): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError("Duplicate keys for " + str(keys))
    return indexed


def adjust(pvalues):
    """Benjamini-Hochberg within a single editor and strategy."""
    ordered = sorted(pvalues, key=pvalues.get, reverse=True)
    q, bound = {}, 1.0
    for index, key in enumerate(ordered):
        p = pvalues[key]
        if not 0 <= p <= 1:
            raise ValueError("Invalid submitted p-value")
        bound = min(bound, p * len(ordered) / (len(ordered)-index))
        q[key] = bound
    return q


def load_panel(args):
    labels = tsv(ORDER)
    targets = [r["target"] for r in labels]
    if len(targets) != 79 or len(set(targets)) != 79:
        raise ValueError("Expected 79 distinct targets in manuscript order")
    paths = [ORDER]
    panel = fields = {}
    if args.results:
        comparison = args.results/"panel_comparison.tsv"
        model = args.results/"model_field_comparison.tsv"
        paths += [comparison, model]
        panel = unique(tsv(comparison), ["target","strategy","editor"])
        fields = unique(tsv(model), ["target","strategy","field"])
        if len(panel) != 632:
            raise ValueError("Expected a complete 632-row panel")
    records = []
    strategy = "Dedud"
    path = args.tables / "SupplementaryTable2.xlsx"
    paths.append(path)
    book = load_workbook(path, read_only=True, data_only=True)
    try:
        it = iter(book.active.values)
        next(it)
        header = next(it)
        rows = unique([dict(zip(header,r)) for r in it if r and r[0]], ["barcode"])
    finally:
        book.close()
    if set(rows) != {(t,) for t in targets}:
        raise ValueError("Workbook targets differ from the submitted design")
    for editor, suffix in EDITORS.items():
        submitted_q = adjust({t:float(rows[t,]["asymPval-treat"+suffix]) for t in targets})
        for position, label in enumerate(labels):
            t = label["target"]
            src = rows[t,]
            keys = ["predicted-NoEP", "predicted-"+suffix, "confLow-"+suffix, "confHigh-"+suffix]
            marginal = [float(src[k]) for k in keys]
            q = submitted_q[t]
            if args.results:
                comp = panel[t,strategy,editor]
                if "manuscript_delta_pp" in comp:
                    if abs(float(comp["manuscript_delta_pp"])-100*(marginal[1]-marginal[0]))>1e-8 or abs(float(comp["manuscript_q_bh"])-q)>1e-8:
                        raise ValueError("Model and workbook describe different submitted estimates")
                    for k, original in zip(keys, marginal):
                        if abs(float(fields[t,strategy,k]["submitted"])-original)>1e-8:
                            raise ValueError("Model-field reference mismatch for " + t)
                marginal = [float(fields[t,strategy,k]["agent"]) for k in keys]
                q = float(comp["q_bh"])
                if abs(float(comp["delta_pp"])-100*(marginal[1]-marginal[0]))>1e-8:
                    raise ValueError("Inconsistent agent estimates for " + t)
            background, estimate, low, high = marginal
            powers = [float(src[f"power{threshold}Treat{suffix}"]) for threshold in ("0.01","0.1")]
            if not (0 <= low <= estimate <= high <= 1 and 0 <= background <= 1 and 0 <= q <= 1
                    and all(0 <= v <= 1 for v in powers)):
                raise ValueError("Invalid model values for " + t)
            records.append(dict(target=t, label=label["label"], order=position+1, guide=t[:4],
                                strategy=strategy, editor=editor, delta_pp=100*(estimate-background),
                                lower_pp=100*(low-background), upper_pp=100*(high-background),
                                q_bh=q, power001=powers[0], power01=powers[1], positive_control=t in CONTROLS))
    return records, {p.name:digest(p) for p in paths}


def editing_plot(rows):
    fig = plt.figure(figsize=(17, 7.7))
    guides = ["1450", "1617"]
    groups = [[r for r in rows if r["guide"]==guide and r["editor"]=="WT"] for guide in guides]
    grid = fig.add_gridspec(5, 2, height_ratios=[1.4,1.45,1.5,2.55,1.1],
                            width_ratios=[len(g) for g in groups], hspace=.09, wspace=.025,
                            left=.084, right=.99, bottom=.115, top=.925)
    lower = min(-.04, math.floor(100*min(r["lower_pp"] for r in rows))/100)
    limits = [(0.4,100), (.04,.4), (lower,.04)]
    lookup = unique(rows, ["target","editor"])
    for col, (guide, group) in enumerate(zip(guides, groups)):
        positions = np.arange(len(group))
        width = len(group)
        for band, (lo,hi) in enumerate(limits):
            ax = fig.add_subplot(grid[band,col])
            for x in range(0,width,2):
                ax.axvspan(x-.5,x+.5,color="#f1f1f1",lw=0)
            for editor, offset in [("WT",-.16),("HiFi",.16)]:
                points = [lookup[r["target"],editor] for r in group]
                for x,r in zip(positions,points):
                    if r["upper_pp"] >= lo and r["lower_pp"] <= hi:
                        ax.vlines(x+offset, max(lo,r["lower_pp"]), min(hi,r["upper_pp"]),color=COLORS[editor],lw=.9)
                visible = [(x+offset,r["delta_pp"]) for x,r in zip(positions,points) if lo <= r["delta_pp"] <= hi]
                if visible:
                    x, y = zip(*visible)
                    ax.scatter(x, y, s=9, color=COLORS[editor], zorder=3)
            for threshold in (.01,.1):
                if lo<=threshold<=hi:
                    ax.axhline(threshold,color="#cf2036",linestyle="--",lw=.85)
            if lo<0:
                ax.axhline(0,color=".48",lw=.6)
            ax.set(xlim=(-.5,width-.5), ylim=(lo,hi),xticks=[])
            if band==0:
                ax.set_title("sgRNA "+guide,fontsize=12,fontweight="bold")
                ax.set_yticks([.4,20,40,60,80,100])
                highlights = [(0, "On-target BCL11A +" + ("55" if col==0 else "58"), 92)]
                if col==1:
                    highlights.append((next(i for i,r in enumerate(group) if r["target"]=="1617_OT_0040_ALT"), "OT-rs114518452 (CPS1)", 49))
                for x, label, y in highlights:
                    ax.text(x+.5,y,label,fontsize=7,va="center")
            elif band==1:
                ax.set_yticks([.04,.1,.2,.3,.4])
                if not col:
                    ax.set_ylabel("Edited − unedited (percentage points)")
            if col:
                ax.tick_params(labelleft=False, left=False)
            ax.tick_params(labelsize=8,length=2)
        names = fig.add_subplot(grid[3,col])
        names.set_xlim(-.5,width-.5)
        names.set_ylim(0,1)
        names.axis("off")
        for i,r in enumerate(group):
            names.text(i,.98,r["label"],ha="center",va="top",rotation=90,fontsize=7)
        ax = fig.add_subplot(grid[4,col])
        for i,r in enumerate(group):
            for editor in EDITORS:
                record = lookup[r["target"],editor]
                for y,value in enumerate([record["q_bh"]<.05, record["power001"], record["delta_pp"]>.01,
                                          record["power01"], record["delta_pp"]>.1]):
                    if y in (1, 3):
                        color = "#ededed" if value is None else POWER(value)
                    else:
                        color = COLORS[editor] if value else "white"
                    vertices = ([(i-.5,y),(i-.5,y+1),(i+.5,y+1)] if editor=="HiFi"
                                else [(i-.5,y),(i+.5,y),(i+.5,y+1)])
                    ax.add_patch(Polygon(vertices, facecolor=color,edgecolor=".15",linewidth=.32))
        ax.set(xlim=(-.5,width-.5),ylim=(0,5),xticks=[])
        ax.set_yticks(np.arange(5)+.5, ["BH q < 0.05","Power at 0.01%","Δ > 0.01 pp","Power at 0.1%","Δ > 0.1 pp"] if not col else [])
        ax.tick_params(length=0,labelsize=8)
    handles = [Line2D([],[],ls="",marker="o",color=v,label=("SpCas9" if k=="WT" else k)) for k,v in COLORS.items()]
    fig.legend(handles=handles,ncol=2,frameon=False,loc="lower left",bbox_to_anchor=(.08,.024))
    bar=fig.colorbar(plt.cm.ScalarMappable(cmap=POWER,norm=plt.Normalize(0,1)),cax=fig.add_axes([.36,.048,.13,.012]),orientation="horizontal")
    bar.set_label("Statistical power",fontsize=8)
    bar.ax.tick_params(labelsize=7)
    fig.suptitle("Figure 3b   |   Dedud", fontsize=15)
    return fig


def write_workbook(path, rows, receipt):
    book = Workbook()
    data = book.active
    data.title = "Figure 3b"
    columns = list(rows[0])
    data.append(columns)
    for row in rows:
        data.append([row[key] for key in columns])
    data.freeze_panes = "C2"
    data.auto_filter.ref = data.dimensions
    for index, width in enumerate([23, 29, 9, 10, 12, 12, 22, 22, 22, 18, 14, 14, 20], 1):
        data.column_dimensions[get_column_letter(index)].width = width
    for cells in data.iter_rows(min_row=2, min_col=7, max_col=12):
        for cell in cells:
            cell.number_format = "0.00000E+00" if cell.column == 10 else "0.00000000"

    info = book.create_sheet("Notes")
    info.append(["Field", "Definition or value"])
    for item in [
        ("estimates", receipt["estimates"]),
        ("power", receipt["power"]),
        ("order", "Submitted target order within the two sgRNA blocks."),
        ("editor", "WT denotes SpCas9. Each editor is compared with NoEP."),
        ("delta_pp", "Edited minus NoEP marginal rate in percentage points."),
        ("lower_pp / upper_pp", "Treatment confidence bounds shifted by the NoEP point estimate, in percentage points."),
        ("q_bh", "Benjamini-Hochberg q-value within the editor and Dedud strategy."),
        ("power001 / power01", "Simulated power at effects of 0.01 and 0.1 percentage points, respectively."),
        ("heatmap", "Effect flags use strict > 0.01 and > 0.1 pp. The significance flag uses q < 0.05."),
        ("plotted_data.tsv SHA-256", receipt["plotted_data_sha256"]),
        *((name + " SHA-256", value) for name, value in receipt["input_sha256"].items()),
    ]:
        info.append(item)
    info.column_dimensions["A"].width = 43
    info.column_dimensions["B"].width = 102
    for sheet in (data, info):
        sheet.sheet_view.showGridLines = False
        for cell in sheet[1]:
            cell.fill = PatternFill("solid", fgColor="204458")
            cell.font = Font(color="FFFFFF", bold=True)
        sheet.row_dimensions[1].height = 23
    book.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure", choices=["3b"])
    parser.add_argument("--tables", type=Path, required=True, help="Directory containing SupplementaryTable2.xlsx")
    parser.add_argument("--results", type=Path, help="Use the completed agent panel refit for estimates and BH q-values")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    rows, sources = load_panel(args)
    save_figure(args.out, rows, sources, "agent refit" if args.results else "submitted workbook",
                "submitted candidate-specific simulations")


def save_figure(output, rows, sources, estimates, power):
    fig = editing_plot(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as tmp:
        dest = Path(tmp) / "result"
        dest.mkdir()
        fig.savefig(dest / "Figure3b.pdf", metadata={"CreationDate": None})
        fig.savefig(dest / "Figure3b.png", dpi=170)
        plt.close(fig)
        data = dest / "plotted_data.tsv"
        with data.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, rows[0], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        receipt = dict(figure="3b", targets=79, rows=len(rows), input_sha256=sources,
                       plotted_data_sha256=digest(data),
                       estimates=estimates,
                       power=power, coordinate_order="submitted Figure 3b")
        write_workbook(dest / "Figure3b_source_data.xlsx", rows, receipt)
        (dest / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
        dest.rename(output)


if __name__ == "__main__":
    main()
