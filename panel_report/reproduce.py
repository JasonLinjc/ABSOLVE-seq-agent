"""Recount all original reporter libraries and compare the manuscript treatment margins."""
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

from maxi_qc import PUBLISHED_SHA256, sha256
from pilot import PUBLIC
from target_identity import TargetIdentity, filtered_counts

AGENT = Path(__file__).resolve().parents[1]
from settings import DEPLOYMENT, MANIFEST, DATA, PUBLISHED
PANEL = AGENT/'panel_report'
ARMS = {'WTCas9':'WT', 'noEP':'NoEP', 'HiFiCas9':'HiFi'}
STRATEGIES = ('Unfiltered','Dedup','Dedud','Dedud+Dedup')
REFERENCE = PUBLIC/'target_oligos_sequenes.csv'


def table(path, delimiter='\t'):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write(path, rows):
    with path.open('x', newline='') as f:
        writer = csv.DictWriter(f, rows[0], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def align_library(sample, donor, replicate, arm, output, workers):
    output.mkdir(parents=True, exist_ok=False)
    config = dict(donor=donor, replicate=replicate, conditions=[arm], samples={arm:sample},
                  maxi_method='manuscript', target_filter='pbc-only',
                  published=str(PUBLISHED))
    (output/'pair.json').write_text(json.dumps(config,indent=2))
    def task(script, index):
        with (output/f'{script}_{index}.log').open('x') as log:
            subprocess.run([sys.executable,str(AGENT/'workers'/(script+'.py')),str(index)],
                           cwd=output,stdout=log,stderr=subprocess.STDOUT,check=True)
    task('prepare',0)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(lambda i:task('target',i),range(79)))


def libraries(stage, workers):
    design = table(PUBLIC/'NovaSeq3_sample_info.csv', ',')
    manifest = {r['sample_name']:r for r in table(MANIFEST) if r['status']=='matched'}
    cohort = DEPLOYMENT/'wt_comparison_261119'
    history = cohort/'results/provenance.json'
    archive = json.loads(history.read_text())['runs'] if history.is_file() else []
    prior = {(r['donor'],r['replicate']):r for r in archive}
    result = []
    for row in design:
        donor, rep, arm = row['donor_id'], row['replicate'], ARMS[row['cas_type']]
        sample = manifest[f'{arm}_pool_{donor}_{rep}']
        if (donor,rep) in prior:
            run = prior[donor,rep]
            if run['samples'][arm] != sample or run['published_sha256'] != PUBLISHED_SHA256:
                raise ValueError('Archived sample or published annotation has changed.')
            root = cohort/'runs'/f'{donor}_{rep}' if arm=='WT' else Path(run['previous_pair'])
        else:
            root = DEPLOYMENT/'paired_agent/panel_report/reproduction_inputs_26'/arm
            archived = root/'pair.json'
            if not archived.exists() or json.loads(archived.read_text())['samples'].get(arm) != sample:
                root = stage/'alignments'/sample['sample_name']
                align_library(sample, donor, rep, arm, root, workers)
        folder = root/arm
        prep = json.loads((folder/'prepared/summary.json').read_text())
        expected_lib = sample['sample_name'] if arm=='WT' else row['LV_ID']
        if not (prep['scope']=='full' and prep['status']=='complete' and prep['counts_conserved'] and prep['sample_id']==expected_lib):
            raise ValueError('Incomplete or mismatched input library: '+sample['sample_name'])
        inputs, = json.loads((folder/'inputs.json').read_text())
        if any(inputs[m] != sample[m] for m in ('r1','r2')):
            raise ValueError('Source FASTQ path mismatch: '+sample['sample_name'])
        result.append(dict(donor=donor,replicate=rep,arm=arm,sample=sample['sample_name'],
                           sha256_r1=sample['sha256_r1'],sha256_r2=sample['sha256_r2'],folder=str(folder)))
    if len(result)!=26 or len({r['sample'] for r in result})!=26:
        raise ValueError('Expected 26 distinct manuscript libraries.')
    return result


def read_good_families(folder, target, good):
    counts = json.loads((folder/'counts.json').read_text())
    if counts['target']!=target or counts.get('target_filter','pbc-only')!='pbc-only':
        raise ValueError('Recount requires archived pre-step-6 counts: '+str(folder))
    families = defaultdict(Counter)
    if counts['prepared']:
        with gzip.open(folder/'CRISPResso_on_target/CRISPResso_output.fastq.gz','rt') as f:
            for header in f:
                f.readline()
                ann = f.readline()
                if not f.readline():
                    raise ValueError('Truncated FASTQ: '+str(folder))
                pbc = header.partition('|pbc=')[2].split('|',1)[0].strip()
                if not good.get(pbc):
                    continue
                fields = dict(re.findall(r'(?:^|\s)([A-Z_]+)=([^\s]*)',ann))
                if fields['ALN']=='NA':
                    continue
                if fields['ALN']!=target:
                    raise ValueError('Unexpected target: '+str(folder))
                edited = fields['CLASS']==target+'_MODIFIED'
                families[pbc][fields['ALN_REF'],fields['ALN_SEQ'],edited] += 1
    return families, counts


def recount_target(task):
    target, samples, published = task
    with (published/(target+'.tsv')).open() as f:
        good = {}
        for line in f:
            pbc, flag = line.rstrip().split('\t')
            good.setdefault(pbc,flag=='1')
    identity = TargetIdentity(REFERENCE)
    rows = []
    for sample in samples:
        families, record = read_good_families(Path(sample['folder'])/'targets'/target,target,good)
        before = record['counts']
        raw, rep, _ = filtered_counts(target,families,good)
        if raw!=before['Dedud'] or rep!=before['Dedud+Dedup'] or record['condition']!=sample['arm']:
            raise ValueError('Pre-filter counts do not match archive: '+str((sample['sample'],target)))
        raw, rep, _ = filtered_counts(target,families,good,identity)
        after = dict(before, Dedud=raw, **{'Dedud+Dedup':rep})
        for strategy in STRATEGIES:
            family_stage = 'Dedud+Dedup' if strategy in ('Dedud','Dedud+Dedup') else 'Dedup'
            k,n = after[strategy]
            rows.append(dict(target=target, strategy=strategy, donor=sample['donor'],replicate=sample['replicate'],
                             arm=sample['arm'], indel=k, denominator=n, family_rows=after[family_stage][1],
                             pre_step6_indel=before[strategy][0],pre_step6_denominator=before[strategy][1]))
    return rows


def run_panel(output):
    workers = int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
    output = Path(output).resolve() if output else AGENT/'runs'/('manuscript_'+os.environ['SLURM_JOB_ID'])
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    published = PUBLISHED
    if (published/'source.sha256').read_text().strip()!=PUBLISHED_SHA256:
        raise ValueError('Published reference SHA-256 mismatch')
    refs = DATA/'submitted_model_reference'
    targets = sorted(TargetIdentity(REFERENCE).designs)
    with TemporaryDirectory(dir=output.parent,prefix='.manuscript-') as tmp:
        stage = Path(tmp)/'result'
        stage.mkdir()
        samples = libraries(stage,workers)
        tasks = [(t,samples,published) for t in targets]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            counts = [r for block in pool.map(recount_target,tasks) for r in block]
        if len(counts)!=26*79*4:
            raise ValueError('Incomplete cohort')
        write(stage/'library_counts.tsv',counts)
        for s in samples:
            path = Path(s['folder'])
            if stage in path.parents:
                s['folder'] = str(output/path.relative_to(stage))
        write(stage/'sample_inventory.tsv',samples)
        subprocess.run(['Rscript',str(PANEL/'manuscript_model.R'),str(stage),str(refs)],check=True)
        receipt = dict(job=os.environ['SLURM_JOB_ID'],libraries=len(samples),targets=len(targets),
                       published_sha256=PUBLISHED_SHA256,oligo_sha256=sha256(REFERENCE),
                       manifest_sha256=sha256(MANIFEST),
                       submitted_workbooks=json.loads((refs/'source_sha256.json').read_text()),
                       implementation_sha256={p.name:sha256(p) for p in (Path(__file__),PANEL/'manuscript_model.R',AGENT/'target_identity.py')})
        (stage/'provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
        stage.rename(output)
    print(output,flush=True)
