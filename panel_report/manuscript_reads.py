"""Paired-read manuscript processing using Trimmomatic 0.36 and CRISPResso2 2.2.13."""
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from historical_maxi import configured_reference, rebuilt_labels
from target_identity import TargetIdentity

STRATEGIES = ("Unfiltered", "Dedup", "Dedud", "Dedud+Dedup")
GITHUB_TRIM_STEPS = ("LEADING:3", "TRAILING:3", "AVGQUAL:20", "MINLEN:48")
HISTORICAL_TRIM_STEPS = {
    "WT": GITHUB_TRIM_STEPS,
    "NoEP": ("LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "AVGQUAL:20", "MINLEN:48", "CROP:105"),
}
HISTORICAL_TRIM_STEPS["HiFi"] = HISTORICAL_TRIM_STEPS["NoEP"]


def table(path, delimiter="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def save(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    with Path(path).open("rb") as handle:
        h = hashlib.sha256()
        for chunk in iter(lambda: handle.read(4 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def checked_version(command, version):
    text = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
    if not re.search(r"(?<![\d.])" + re.escape(version) + r"(?![\d.])", text):
        raise ValueError(f"Expected version {version} for {command[0]}")


def trim(cfg, sample, dest, workers):
    jar = cfg["jar"]
    checked_version(["java", "-jar", jar, "-version"], "0.36")
    for mate in ("r1", "r2"):
        if digest(sample[mate]) != sample["sha256_" + mate]:
            raise ValueError("Input checksum mismatch: " + sample[mate])
    command = ["java", "-Xmx3g", "-jar", jar, "PE", "-threads", str(workers), "-phred33",
               sample["r1"], sample["r2"], str(dest/"R1.gz"), str(dest/"unpaired1.gz"),
               str(dest/"R2.gz"), str(dest/"unpaired2.gz"),
               "ILLUMINACLIP:" + cfg["adapters"] + ":2:30:10",
               *sample.get("trim_steps", GITHUB_TRIM_STEPS)]
    with (dest/"trim.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True)
    for f in ("unpaired1.gz", "unpaired2.gz"):
        (dest/f).unlink()


def records(handle):
    for header in handle:
        seq, plus, quality = handle.readline().rstrip(), handle.readline(), handle.readline().rstrip()
        if not header.startswith(b"@") or not plus.startswith(b"+") or len(seq) != len(quality):
            raise ValueError("Invalid FASTQ record")
        yield header.split()[0], seq, quality


def demultiplex(cfg, sample_dir, dest):
    designs = table(Path(cfg["reference"])/"target_oligos_sequenes.csv", ",")
    barcodes = {r["Barcode"].encode(): r["OT"].replace("-", "_") for r in designs}
    constant = b"CCAACCTCATAGAACACTCATCC"
    counts = Counter()
    with ExitStack() as stack:
        inputs = [records(stack.enter_context(gzip.open(sample_dir/"trim"/f"R{i}.gz", "rb"))) for i in (1,2)]
        writers = {}
        for one, two in zip_longest(*inputs):
            if not one or not two or one[0].removesuffix(b"/1") != two[0].removesuffix(b"/2"):
                raise ValueError("Paired FASTQ records differ")
            seq = one[1]
            fixed = seq[16:39]
            target = barcodes.get(seq[39:48])
            if not target or (fixed != constant and
                    (fixed.count(b"N") > 5 or any(x != y and x != 78 for x,y in zip(fixed,constant)))):
                continue
            if target not in writers:
                writers[target] = [stack.enter_context((dest/(target+f".R{i}.fq")).open("wb")) for i in (1,2)]
            counts[target] += 1
            header = b"@" + str(counts[target]).encode() + b"|" + seq[:16] + b"\n"
            for handle, record, region in zip(writers[target], (one,two), (slice(16,None),slice(None,-16))):
                handle.write(b"".join((header,record[1][region],b"\n+\n",record[2][region],b"\n")))
    (dest/"pairs.json").write_text(json.dumps(counts,sort_keys=True))


def step(func, sample_dir, name, *args):
    final = sample_dir/name
    if final.is_dir():
        return
    sample_dir.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=sample_dir) as tmp:
        stage=Path(tmp)/name
        stage.mkdir()
        func(*args,stage)
        stage.rename(final)


def worker(config, action, index, threads=None):
    cfg=json.loads(config.read_text())
    sample=cfg["libraries"][index]
    root=config.parent/"work"/sample["sample_name"]
    if threads is None:
        threads = int(os.getenv("SLURM_CPUS_PER_TASK", "1"))
    if action=="prepare":
        step(lambda stage: trim(cfg,sample,stage,threads),root,"trim")
        step(demultiplex,root,"pairs",cfg,root)
    elif action=="align":
        align_library(cfg,sample,root,threads)
    else:
        raise ValueError(action)


def summarize(target, families, labels, identity):
    raw = [0,0]
    rep = [0,0]
    pre = [0,0]
    pre_rep = [0,0]
    final = [0,0]
    final_rep = [0,0]
    for barcode, alleles in families.items():
        # v2.2.13 sorts each family's allele table by abundance, aligned read, then reference.
        order = sorted(alleles, key=lambda a: (-alleles[a], a[1], a[0]))
        k = sum(n for a,n in alleles.items() if a[2])
        n = sum(alleles.values())
        raw[0] += k
        raw[1] += n
        rep[0] += int(order[0][2])
        rep[1] += 1
        if not labels.get(barcode,False):
            continue
        pre[0] += k
        pre[1] += n
        pre_rep[0] += int(order[0][2])
        pre_rep[1] += 1
        keep = [a for a in order if not a[2] or identity.retain(target,a)]
        if keep:
            final[0] += sum(alleles[a] for a in keep if a[2])
            final[1] += sum(alleles[a] for a in keep)
            final_rep[0] += int(keep[0][2])
            final_rep[1] += 1
    return dict(zip(STRATEGIES,[raw,rep,final,final_rep])), dict(zip(STRATEGIES,[raw,rep,pre,pre_rep]))


def crispresso_command(cfg, root=None):
    if not cfg.get("image"):
        return [cfg["crispresso"]]
    command = ["apptainer", "exec", "--cleanenv"]
    if root is not None:
        command.extend(["--bind", str(root)])
    return command + [cfg["image"], "CRISPResso"]


def read_aligned_families(results, target, pair_count):
    """Read paired-end CRISPResso annotations and reconcile them with its totals."""
    families = defaultdict(Counter)
    report = results/"CRISPResso_quantification_of_editing_frequency.txt"
    cache = {}
    with gzip.open(results/"CRISPResso_output.fastq.gz","rb") as f:
        for header in f:
            f.readline()
            ann = f.readline()
            if not f.readline():
                raise ValueError("Truncated CRISPResso FASTQ")
            if b" ALN=NA " in ann:
                continue
            # Repeated reads share the same short alignment annotation.
            allele = cache.get(ann)
            if allele is None:
                fields = dict(re.findall(rb"(?:^|\s)([A-Z_]+)=([^\s]*)",ann))
                if fields[b"ALN"] != b"Reference" or fields[b"CLASS"] not in (b"Reference_MODIFIED",b"Reference_UNMODIFIED"):
                    raise ValueError("Unexpected aligned read annotation")
                allele=(fields[b"ALN_REF"].decode(),fields[b"ALN_SEQ"].decode(),fields[b"CLASS"]==b"Reference_MODIFIED")
                cache[ann]=allele
            barcode=header.split()[0].split(b"|")[-1].decode()
            if len(barcode)!=16:
                raise ValueError("Lost original pBC during FLASH merging")
            families[barcode][allele]+=1
    q,=table(report)
    k=sum(n for vs in families.values() for a,n in vs.items() if a[2])
    total=sum(sum(vs.values()) for vs in families.values())
    if (k,total,pair_count)!=(int(q["Modified"]),int(q["Reads_aligned"]),int(q["Reads_in_input"])):
        raise ValueError(f"Read conservation failed for {target}: {k,total,pair_count} vs {q}")
    return families


def count_target(cfg, root, spec, pair_count):
    target = spec["OT_name"]
    final = root/"targets"/(target+".json")
    rebuilt = configured_reference(cfg)
    reference_hash = rebuilt.get("reference_sha256")
    if final.exists():
        if json.loads(final.read_text()).get("maxi_reference_sha256") != reference_hash:
            raise ValueError("Cached counts use a different Maxi reference")
        return
    families = defaultdict(Counter)
    with TemporaryDirectory(dir=root/"targets",prefix="."+target+"-") as tmp:
        work = Path(tmp)
        if pair_count:
            for i in (1,2):
                (work/f"R{i}.fq").symlink_to(root/"pairs"/f"{target}.R{i}.fq")
            command = crispresso_command(cfg, root) + ["--fastq_r1","R1.fq","--fastq_r2","R2.fq",
                       "--amplicon_seq",spec["amplicon_seq"],"--guide_seq",spec["guide_seq"],
                       "--ignore_substitutions","--plot_window_size","10","--min_frequency_alleles_around_cut_to_plot","0","--name","family","--flash_command","flash -t 1",
                       "--write_detailed_allele_table","--fastq_output","--suppress_plots","--suppress_report"]
            with (work/"crispresso.log").open("w") as log:
                proc = subprocess.run(command,cwd=work,stdout=log,stderr=log)
            results = work/"CRISPResso_on_family"
            report = results/"CRISPResso_quantification_of_editing_frequency.txt"
            if proc.returncode or not report.is_file():
                log = (work/"crispresso.log").read_text()
                # Upstream Batch --skip_failed omits families with no merged or aligned reads.
                if "No reads" not in log or not re.search("No reads (were aligned|aligned|survived)",log,re.I):
                    raise RuntimeError(target+" CRISPResso failed: "+log[-3000:])
            else:
                families = read_aligned_families(results, target, pair_count)
        if rebuilt:
            labels = rebuilt_labels(rebuilt, target)
        else:
            labels={}
            with (Path(cfg["published"])/(target+".tsv")).open() as handle:
                for ln in handle:
                    b,flag=ln.rstrip().split("\t");labels.setdefault(b,flag=="1")
        identity=TargetIdentity(Path(cfg["reference"])/"target_oligos_sequenes.csv")
        counts, before=summarize(target,families,labels,identity)
        if target in ("1617_OT_0041_ALT","1450_OT_0274_ALT"):
            with gzip.open(work/"focal.tsv.gz","wt") as handle:
                handle.write("pbc_hash\tindel\tdenominator\n")
                for b,vs in families.items():
                    if not labels.get(b):continue
                    retained={a:n for a,n in vs.items() if identity.retain(target,a)}
                    if retained:
                        tag=hashlib.sha256(b.encode()).hexdigest()[:16]
                        handle.write(f"{tag}\t{sum(n for a,n in retained.items() if a[2])}\t{sum(retained.values())}\n")
            (work/"focal.tsv.gz").replace(final.with_suffix(".families.tsv.gz"))
        receipt=dict(target=target,paired_input=pair_count,counts=counts,before=before,
                     maxi_reference_sha256=reference_hash)
        (work/"result").write_text(json.dumps(receipt))
        (work/"result").replace(final)


def align_library(cfg,sample,root,threads):
    checked_version(crispresso_command(cfg) + ["--version"], "2.2.13")
    (root/"targets").mkdir(exist_ok=True)
    specs=table(Path(cfg["reference"])/"OT_guide_amplicon_seq.csv",",")
    admitted=json.loads((root/"pairs/pairs.json").read_text())
    if len({s["OT_name"] for s in specs}) != 79:
        raise ValueError("Incomplete target design")
    specs.sort(key=lambda s:-admitted.get(s["OT_name"],0))
    with ThreadPoolExecutor(max_workers=threads) as pool:
        list(pool.map(lambda s:count_target(cfg,root,s,admitted.get(s["OT_name"],0)),specs))


if __name__=="__main__":
    worker(Path(sys.argv[1]).resolve(),sys.argv[2],int(sys.argv[3]))
