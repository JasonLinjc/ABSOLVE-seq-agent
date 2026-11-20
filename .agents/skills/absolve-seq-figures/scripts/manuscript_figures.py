"""Regenerate the submitted ABSOLVE-seq editing and power figures."""
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
from matplotlib.patches import Polygon, Patch
from openpyxl import load_workbook
import numpy as np

ORDER = Path(__file__).resolve().parents[1] / "references/target_order.tsv"
FIGURES = {"3b": ("Dedud", 2), "S10": ("Unfiltered", 3), "S11": ("Dedup", 4),
           "S12": ("Dedud+Dedup", 5)}
STRATEGIES = ["Unfiltered", "Dedup", "Dedud", "Dedud+Dedup"]
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


def load_panel(args, chosen):
    labels = tsv(ORDER)
    targets = [r["target"] for r in labels]
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
    for strategy, number in chosen:
        path = args.tables/f"SupplementaryTable{number}.xlsx"
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


def load_summary_archive(folder):
    """Follow the original visualize.ipynb inputs for both halves of Supplementary Figure 13."""
    stems = ["raw","raw_dedup","dedudS_filtered_v3","dedudS_filtered_v3_dedup"]
    labels = tsv(ORDER)
    targets = {(row["target"],) for row in labels}
    rows = []
    paths = [ORDER]
    for strategy, stem in zip(STRATEGIES, stems):
        path = folder / ("summaryResult_" + stem + "_fdr.tsv")
        archive = unique(tsv(path), ["offtarget"])
        paths.append(path)
        if archive.keys() != targets:
            raise ValueError("Unexpected original figure target set")
        powers = {}
        for threshold in ["0.01", "0.1"]:
            power_file = folder / f"powerAnalysis_{stem}_{threshold}.csv"
            paths.append(power_file)
            with power_file.open() as stream:
                data = unique(list(csv.DictReader(stream)), ["sequence"])
            if data.keys() != targets:
                raise ValueError("Unexpected original power target set")
            powers[threshold] = data
        for editor, suffix in EDITORS.items():
            for position, label in enumerate(labels):
                t = label["target"]
                q = float(archive[t,]["adj_pval_"+suffix])
                delta = float(archive[t,]["Delta"+suffix])
                p1, p2 = [float(powers[key][t,]["Treat"+suffix]) for key in ["0.01","0.1"]]
                if not (0 <= q <= 1 and 0 <= p1 <= 1 and 0 <= p2 <= 1 and math.isfinite(delta)):
                    raise ValueError("Invalid archived figure estimates")
                rows.append(dict(target=t, label=label["label"], order=position+1, guide=t[:4],
                                 strategy=strategy,editor=editor,delta_pp=delta,q_bh=q,power001=p1,
                                 power01=p2,positive_control=t in CONTROLS))
    return rows, {p.name:digest(p) for p in paths}


def editing_plot(rows, name):
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
                    color = POWER(value) if y in (1,3) else COLORS[editor] if value else "white"
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
    fig.suptitle(f"{'Figure 3b' if name=='3b' else 'Supplementary Figure '+name[1:]}   |   {rows[0]['strategy']}",fontsize=15)
    return fig


def summary_plot(rows):
    fig, axes = plt.subplots(2,2,figsize=(10,7.4),layout="constrained")
    for j,(delta,key) in enumerate([(.1,"power01"),(.01,"power001")]):
        bars, dots = axes[j]
        for si, strategy in enumerate(STRATEGIES):
            for i, editor in enumerate(EDITORS):
                selected = [r for r in rows if r["strategy"]==strategy and r["editor"]==editor]
                x = si+(i-.5)*.32
                hits = [r for r in selected if r["q_bh"]<.05 and r["delta_pp"]>delta]
                bars.bar(x,len(hits),width=.31,color=COLORS[editor],alpha=.72)
                bars.bar(x,sum(r["positive_control"] for r in hits),width=.31,color=COLORS[editor],ec=".25",lw=.7)
                for k,r in enumerate(selected):
                    control=r["positive_control"]
                    dots.scatter(x+(((k*31)%79)/78-.5)*.23,r[key],s=15 if control else 7,
                                 marker="^" if control else "o",color="red" if control else COLORS[editor])
        bars.set_ylim(0,35)
        bars.set_ylabel(f"Targets with Δ > {delta}% and BH q < 0.05")
        dots.set_ylabel(f"Power to detect Δ = {delta}%")
        dots.set_ylim(-.025,1.035)
        for value,style in [(.8,"-"),(.9,"--")]:
            dots.axhline(value,color="#cc3646",ls=style,lw=.8)
        for ax in (bars,dots):
            ax.set_xticks(range(4),STRATEGIES,rotation=18)
            ax.spines[["top","right"]].set_visible(False)
        bars.set_title("a" if j==0 else "b",loc="left",fontweight="bold")
    fig.suptitle("Supplementary Figure 13   |   Significant candidates and statistical power")
    bars.legend(handles=[Patch(color=COLORS["WT"],label="SpCas9"),Patch(color=COLORS["HiFi"],label="HiFi"),
                         Line2D([],[],ls="",marker="^",color="red",label="Positive control")],fontsize=8,frameon=False)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure",choices=[*FIGURES,"S13"])
    parser.add_argument("--tables",type=Path,help="Directory containing submitted SupplementaryTable2–5.xlsx")
    parser.add_argument("--figure-tables",type=Path,help="Original draft5 tables directory for the artwork version of S13")
    parser.add_argument("--results",type=Path,help="Use refitted estimates and BH q from the agent model directory")
    parser.add_argument("--out",type=Path,required=True,help="New output directory")
    args=parser.parse_args()
    chosen = [FIGURES[args.figure]] if args.figure in FIGURES else list(FIGURES.values())
    if args.figure_tables:
        if args.figure != "S13" or args.tables or args.results:
            parser.error("--figure-tables applies only to S13 and selects the original artwork inputs")
        rows, sources = load_summary_archive(args.figure_tables)
        estimate_source = "draft5 original figure tables"
    else:
        if not args.tables:
            parser.error("Supply --tables, or --figure-tables for original S13")
        rows, sources = load_panel(args, chosen)
        estimate_source = "agent refit" if args.results else "submitted workbook"
    fig = summary_plot(rows) if args.figure=="S13" else editing_plot(rows,args.figure)
    if args.out.exists():
        raise FileExistsError(args.out)
    if args.figure=="S13" and not args.figure_tables:
        fig.suptitle("Supplementary Figure 13   |   "+estimate_source+" estimates")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.out.parent) as tmp:
        dest = Path(tmp) / "result"
        dest.mkdir()
        stem = "Figure3b" if args.figure=="3b" else "SupplementaryFigure"+args.figure[1:]
        fig.savefig(dest/(stem+".pdf"))
        fig.savefig(dest/(stem+".png"),dpi=170)
        plt.close(fig)
        with (dest/"plotted_data.tsv").open("w") as f:
            writer=csv.DictWriter(f,rows[0],delimiter="\t",lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        receipt = dict(figure=args.figure,targets=79,rows=len(rows),input_sha256=sources,
                       estimates=estimate_source,
                       power="draft5 power tables" if args.figure_tables else "submitted candidate-specific simulations",coordinate_order="submitted Figure 3b")
        (dest/"provenance.json").write_text(json.dumps(receipt,indent=2)+"\n")
        dest.rename(args.out)


if __name__=="__main__":
    main()
