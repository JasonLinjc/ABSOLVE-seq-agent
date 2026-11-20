"""Submit and summarize the original paired-end ABSOLVE-seq workflow."""
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from panel_report.manuscript_reads import GITHUB_TRIM_STEPS, HISTORICAL_TRIM_STEPS, ROOT, STRATEGIES, digest, save, table
from analysis_report import complete_report
from settings import DATA, MANIFEST, PUBLIC, PUBLISHED, RUNS

ARM = {"WTCas9":"WT","HiFiCas9":"HiFi","noEP":"NoEP"}
PANEL = ROOT/"panel_report"


def compare(baseline, current, output):
    def index(path):
        rows=table(path)
        result={(r["target"],r["strategy"],r["editor"]):r for r in rows}
        if len(rows)!=632 or len(result)!=632:
            raise ValueError("Expected 632 distinct panel comparisons")
        return result
    old,new=index(baseline),index(current)
    if old.keys()!=new.keys():
        raise ValueError("Target panels differ")
    records=[]
    for key,r in sorted(new.items()):
        p=old[key]
        if any(abs(float(r[k])-float(p[k]))>1e-9 for k in ("manuscript_delta_pp","manuscript_q_bh")):
            raise ValueError("The baseline uses different manuscript estimates")
        before=p["hit_agrees"]=="TRUE"
        after=r["hit_agrees"]=="TRUE"
        status={ (True,True):"concordant",(False,True):"resolved",(False,False):"persistent",(True,False):"new"}[before,after]
        records.append(dict(target=key[0],strategy=key[1],editor=key[2],status=status,
            submitted_hit=r["submitted_hit"],fast_hit=p["agent_hit"],rerun_hit=r["agent_hit"],
            submitted_delta_pp=r["manuscript_delta_pp"],fast_delta_pp=p["delta_pp"],rerun_delta_pp=r["delta_pp"],
            submitted_q=r["manuscript_q_bh"],fast_q=p["q_bh"],rerun_q=r["q_bh"]))
    save(output,records)


def fit_panel(stage, refs, image):
    refs = refs.resolve(strict=True)
    margin = "equal-donor"
    command = ["Rscript", str(PANEL / "manuscript_model.R"), str(stage), str(refs), margin]
    if image:
        binds = ",".join(str(p) for p in sorted({stage.resolve(), refs.resolve(), PANEL.resolve()}))
        command = ["apptainer", "exec", "--cleanenv", "--bind", binds, image, *command]
    subprocess.run(command, check=True)
    return dict(donor_margin=margin, submitted_workbooks=json.loads((refs / "source_sha256.json").read_text()))


def refit_counts(source, output):
    """Fit the submitted model from frozen library counts with their existing filtering."""
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    refs = DATA / "submitted_model_reference"
    design = {(r["donor_id"], r["replicate"], ARM[r["cas_type"]])
              for r in table(PUBLIC / "NovaSeq3_sample_info.csv", ",")}
    targets = {r["OT_name"] for r in table(PUBLIC / "OT_guide_amplicon_seq.csv", ",")}
    expected = {(t, strategy, *sample) for t in targets for strategy in STRATEGIES for sample in design}
    with TemporaryDirectory(dir=output.parent) as tmp:
        stage = Path(tmp) / "results"
        stage.mkdir()
        counts = stage / "library_counts.tsv"
        shutil.copyfile(source, counts)
        rows = table(counts)
        keys = {(r["target"], r["strategy"], r["donor"], r["replicate"], r["arm"]) for r in rows}
        if len(rows) != 8216 or expected != keys or len(keys) != len(rows):
            raise ValueError("Expected all 8,216 distinct original panel count cells")
        receipt = fit_panel(stage, refs, os.environ.get("ABSOLVE_IMAGE"))
        receipt.update(mode="frozen-count-refit", count_sha256=digest(counts),
                       model_sha256=digest(PANEL / "manuscript_model.R"))
        (stage / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
        complete_report(stage, "frozen-count refit with equal donor weights", panel=True, image=os.environ.get("ABSOLVE_IMAGE"))
        stage.rename(output)
    print(output)
    print(output / "report.md", flush=True)


def finish(config):
    cfg=json.loads(config.read_text()); base=config.parent
    targets=sorted(r["OT_name"] for r in table(Path(cfg["reference"])/"OT_guide_amplicon_seq.csv",","))
    libs=cfg["libraries"]
    if len(targets)!=len(set(targets)) or len(targets)!=79:
        raise ValueError("Incomplete reference")
    rows=[]; paired={}
    for s in libs:
        for t in targets:
            record=json.loads((base/"work"/s["sample_name"]/"targets"/(t+".json")).read_text())
            if record.get("maxi_reference_sha256") != cfg.get("historical_reference", {}).get("reference_sha256"):
                raise ValueError("Target counts use a different Maxi reference")
            if record["target"]!=t: raise ValueError("Target mismatch")
            for strategy in STRATEGIES:
                n,d=record["counts"][strategy]
                a,b=record["before"][strategy]
                unit="Dedud+Dedup" if strategy.startswith("Dedud") else "Dedup"
                rows.append(dict(target=t,strategy=strategy,donor=s["donor"],replicate=s["replicate"],arm=s["arm"],
                                 indel=n,denominator=d,family_rows=record["counts"][unit][1],
                                 pre_step6_indel=a,pre_step6_denominator=b))
                if len(libs)==2: paired[t,strategy,s["arm"]]=n,d
    with TemporaryDirectory(dir=base) as tmp:
        stage=Path(tmp)/"results";stage.mkdir()
        save(stage/"library_counts.tsv",rows)
        inventory=[]
        for sample in libs:
            entry={k:sample[k] for k in ("sample_name","donor","replicate","arm","sha256_r1","sha256_r2")}
            entry["trim_steps"]=" ".join(sample.get("trim_steps", GITHUB_TRIM_STEPS))
            inventory.append(entry)
        save(stage/"sample_inventory.tsv",inventory)
        receipt=dict(read_mode="manuscript",crispresso_version="2.2.13",trimmomatic_version="0.36",
                     jar_sha256=digest(cfg["jar"]),adapters_sha256=digest(cfg["adapters"]),
                     libraries=len(libs),targets=len(targets),input_config_sha256=digest(config),
                     worker_sha256=digest(PANEL/"manuscript_reads.py"),model_sha256=digest(PANEL/"manuscript_model.R"))
        if cfg.get("image"):
            receipt["image_sha256"] = digest(cfg["image"])
        else:
            receipt["crispresso_executable_sha256"] = digest(cfg["crispresso"])
        receipt["executor"] = cfg.get("executor", "slurm")
        receipt["maxi_method"] = cfg.get("maxi_method", "manuscript")
        if "historical_reference" in cfg:
            receipt["historical_reference"] = cfg["historical_reference"]
        if len(libs)==26:
            refs=Path(cfg["refs"]).resolve(strict=True)
            receipt.update(fit_panel(stage, refs, cfg.get("stats_image")))
            if cfg.get("baseline"):
                compare(Path(cfg["baseline"]),stage/"panel_comparison.tsv",stage/"call_changes.tsv")
                receipt["baseline_sha256"]=digest(cfg["baseline"])
        elif len(libs)==2:
            rates=[]
            for t in targets:
                for strat in STRATEGIES:
                    r=dict(donor=libs[0]["donor"], replicate=libs[0]["replicate"],
                           target=t, strategy=strat, maxi_method=cfg.get("maxi_method", "manuscript"),
                           min_exact_reads=None, target_filter="manuscript")
                    for arm in ("NoEP","HiFi"):
                        n,d=paired[t,strat,arm]
                        r.update({arm+"_indel":n,arm+"_denominator":d,arm+"_rate_pct":100*n/d if d else None})
                    r["delta_HiFi_minus_NoEP_pp"]=(r["HiFi_rate_pct"]-r["NoEP_rate_pct"] if
                        r["HiFi_rate_pct"] is not None and r["NoEP_rate_pct"] is not None else None)
                    rates.append(r)
            save(stage/"paired_target_rates.tsv",rates)
        else:
            raise ValueError("Unexpected manuscript library design")
        (stage/"provenance.json").write_text(json.dumps(receipt,indent=2)+"\n")
        if len(libs) == 26:
            complete_report(stage, "manuscript reads with equal donor weights", panel=True, run_config="run.json", image=cfg.get("stats_image"))
        else:
            complete_report(stage, "manuscript", run_config="run.json", image=cfg.get("stats_image"))
        stage.rename(base/"results")
    print(base / "results/report.md", flush=True)


def submit(config, stages=("prepare","align","finish")):
    cfg=json.loads(config.read_text())
    workers=PANEL/"reads.sbatch"; dependency=None
    env=dict(os.environ,ABSOLVE_HOST_PYTHON=sys.executable)
    jobs={}
    for phase in stages:
        args=["sbatch","--parsable","--chdir",str(config.parent),"--output",str(config.parent/(phase+"-%A_%a.log"))]
        if dependency:
            args.append("--dependency="+("aftercorr:" if phase=="align" else "afterok:")+dependency)
        if phase=="finish":
            args.extend(["--cpus-per-task=1","--mem=4G","--time=01:00:00","--wrap",
                         shlex.join([sys.executable,str(ROOT/"paired.py"),"--finish-manuscript",str(config)])])
        else:
            args.extend(["--array",f"0-{len(cfg['libraries'])-1}",str(workers),str(PANEL/"manuscript_reads.py"),str(config),phase])
        dependency=subprocess.check_output(args,text=True,env=env).strip().split(";")[0]
        jobs[phase]=dependency
        (config.parent/"jobs.json").write_text(json.dumps(jobs,indent=2)+"\n")
        print(f"{phase} job {dependency}",flush=True)


def start(args, *, local, workers):
    manifest={}
    for r in table(MANIFEST):
        if r["status"]=="matched":
            if r["sample_name"] in manifest: raise ValueError("Duplicate sample: "+r["sample_name"])
            manifest[r["sample_name"]]=r
    libs=[]
    for r in table(PUBLIC/"NovaSeq3_sample_info.csv",","):
        d,rep,arm=r["donor_id"],r["replicate"],ARM[r["cas_type"]]
        if not args.manuscript_panel and (d!=args.donor or rep!=args.replicate or arm=="WT"): continue
        libs.append(dict(manifest[f"{arm}_pool_{d}_{rep}"],donor=d,replicate=rep,arm=arm,
                         trim_steps=HISTORICAL_TRIM_STEPS[arm]))
    if len(libs)!=(26 if args.manuscript_panel else 2) or len({s["sample_name"] for s in libs})!=len(libs):
        raise ValueError("Incomplete or duplicated manuscript design")
    from maxi_qc import PUBLISHED_SHA256
    if not args.maxi_reference and (PUBLISHED/"source.sha256").read_text().strip()!=PUBLISHED_SHA256:
        raise ValueError("Published Maxi annotation checksum mismatch")
    cfg=dict(libraries=libs,reference=str(PUBLIC),published=str(PUBLISHED),refs=str(DATA/"submitted_model_reference"),
             adapters=str(DATA/"upstream/test/data/adapters.fa"),stats_image=os.getenv("ABSOLVE_IMAGE"),
             baseline=str(Path(args.baseline_panel).resolve()) if args.baseline_panel else None)
    cfg["maxi_method"] = args.maxi_method
    if args.maxi_reference:
        from historical_maxi import reference
        cfg["historical_reference"] = reference(args.maxi_reference)
        del cfg["published"]
    jar = os.environ.get("ABSOLVE_TRIMMOMATIC_JAR")
    if not jar:
        raise ValueError("Set ABSOLVE_TRIMMOMATIC_JAR to the Trimmomatic 0.36 JAR")
    cfg["jar"] = str(Path(jar).expanduser().resolve(strict=True))
    cfg["executor"] = "local" if local else "slurm"
    image = os.environ.get("ABSOLVE_CRISPRESSO2213_IMAGE")
    if local:
        cfg["threads"] = workers
    if image:
        if local and not shutil.which("apptainer"):
            raise ValueError("Install Apptainer for the configured image, or unset ABSOLVE_CRISPRESSO2213_IMAGE and configure native CRISPResso2 2.2.13")
        cfg["image"] = str(Path(image).expanduser().resolve(strict=True))
    elif local:
        from panel_report.manuscript_reads import checked_version
        executable = shutil.which(os.getenv("ABSOLVE_CRISPRESSO2213_EXECUTABLE", "CRISPResso"))
        if not executable or not shutil.which("flash"):
            raise ValueError("Local manuscript mode requires CRISPResso2 2.2.13 and flash, or ABSOLVE_CRISPRESSO2213_IMAGE")
        cfg["crispresso"] = str(Path(executable).absolute())
        checked_version([cfg["crispresso"], "--version"], "2.2.13")
    else:
        raise ValueError("Set ABSOLVE_CRISPRESSO2213_IMAGE")
    inputs=[cfg["adapters"]]+[s[k] for s in libs for k in ("r1","r2")]
    if args.manuscript_panel:
        refs = Path(cfg["refs"])
        inputs.extend(refs / f"SupplementaryTable{i}.csv" for i in range(2, 6))
        inputs.append(refs / "source_sha256.json")
        if cfg["stats_image"]:
            cfg["stats_image"] = str(Path(cfg["stats_image"]).expanduser().resolve(strict=True))
    if cfg["baseline"]: inputs.append(cfg["baseline"])
    for path in inputs:
        if not Path(path).is_file():raise FileNotFoundError(path)
    if args.output:
        base=Path(args.output).resolve();base.mkdir(parents=True)
    else:
        import tempfile
        RUNS.mkdir(parents=True,exist_ok=True)
        base=Path(tempfile.mkdtemp(prefix="manuscript_reads_",dir=RUNS))
    config=base/"run.json";config.write_text(json.dumps(cfg,indent=2)+"\n")
    print("Manuscript run directory",base,flush=True)
    if local:
        run_local(config, len(libs), workers)
    else:
        submit(config)


def run_local(config, library_count, workers):
    from panel_report.manuscript_reads import worker
    for phase in ("prepare", "align"):
        for index in range(library_count):
            print(f"{phase} library {index + 1}/{library_count}", flush=True)
            worker(config, phase, index, threads=workers)
    finish(config)
    print(config.parent / "results", flush=True)
