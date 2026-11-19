"""Submit and summarize the original paired-end ABSOLVE-seq workflow."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from tempfile import TemporaryDirectory

from panel_report.manuscript_reads import ROOT, STRATEGIES, digest, save, table
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
        save(stage/"sample_inventory.tsv",[{k:s[k] for k in ("sample_name","donor","replicate","arm","sha256_r1","sha256_r2")} for s in libs])
        receipt=dict(read_mode="manuscript",crispresso_version="2.2.13",trimmomatic_version="0.36",
                     image_sha256=digest(cfg["image"]),jar_sha256=digest(cfg["jar"]),adapters_sha256=digest(cfg["adapters"]),
                     libraries=len(libs),targets=len(targets),input_config_sha256=digest(config),
                     worker_sha256=digest(PANEL/"manuscript_reads.py"),model_sha256=digest(PANEL/"manuscript_model.R"))
        if len(libs)==26:
            refs=Path(cfg["refs"])
            command=["Rscript",str(PANEL/"manuscript_model.R"),str(stage),str(refs)]
            if cfg.get("stats_image"):
                binds=",".join(str(p) for p in sorted({stage.resolve(),refs.resolve(),PANEL.resolve()}))
                command=["apptainer","exec","--cleanenv","--bind",binds,cfg["stats_image"],*command]
            subprocess.run(command,check=True)
            receipt["submitted_workbooks"]=json.loads((refs/"source_sha256.json").read_text())
            if cfg.get("baseline"):
                compare(Path(cfg["baseline"]),stage/"panel_comparison.tsv",stage/"call_changes.tsv")
                receipt["baseline_sha256"]=digest(cfg["baseline"])
        elif len(libs)==2:
            rates=[]
            for t in targets:
                for strat in STRATEGIES:
                    r=dict(donor=libs[0]["donor"], replicate=libs[0]["replicate"],
                           target=t, strategy=strat, maxi_method="manuscript",
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
        stage.rename(base/"results")


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


def start(args):
    manifest={}
    for r in table(MANIFEST):
        if r["status"]=="matched":
            if r["sample_name"] in manifest: raise ValueError("Duplicate sample: "+r["sample_name"])
            manifest[r["sample_name"]]=r
    libs=[]
    for r in table(PUBLIC/"NovaSeq3_sample_info.csv",","):
        d,rep,arm=r["donor_id"],r["replicate"],ARM[r["cas_type"]]
        if not args.manuscript_panel and (d!=args.donor or rep!=args.replicate or arm=="WT"): continue
        libs.append(dict(manifest[f"{arm}_pool_{d}_{rep}"],donor=d,replicate=rep,arm=arm))
    if len(libs)!=(26 if args.manuscript_panel else 2) or len({s["sample_name"] for s in libs})!=len(libs):
        raise ValueError("Incomplete or duplicated manuscript design")
    from maxi_qc import PUBLISHED_SHA256
    if (PUBLISHED/"source.sha256").read_text().strip()!=PUBLISHED_SHA256:
        raise ValueError("Published Maxi annotation checksum mismatch")
    cfg=dict(libraries=libs,reference=str(PUBLIC),published=str(PUBLISHED),refs=str(DATA/"submitted_model_reference"),
             adapters=str(DATA/"upstream/test/data/adapters.fa"),stats_image=os.getenv("ABSOLVE_IMAGE"),
             baseline=str(Path(args.baseline_panel).resolve()) if args.baseline_panel else None)
    for key,var in [("jar","ABSOLVE_TRIMMOMATIC_JAR"),("image","ABSOLVE_CRISPRESSO2213_IMAGE")]:
        if not os.getenv(var):raise ValueError("Set "+var)
        cfg[key]=str(Path(os.environ[var]).resolve(strict=True))
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
    submit(config)
