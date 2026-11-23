"""Align one sample/target and count all four manuscript strategies."""
import csv
import json
import itertools
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

root=Path(__file__).resolve().parent
sys.path.insert(0,str(root.parent))
from tools.pilot import CUT, PUBLIC, fastq, one
from absolve.target_identity import TargetIdentity, filtered_counts, family_counts
from absolve.maxi_qc import sha256, PUBLISHED_SHA256, REQUIRED, DEFAULT_METHOD
from absolve.historical_maxi import configured_reference, rebuilt_labels
from absolve.settings import RUNTIME
sys.path.insert(0,str(RUNTIME/'src'))
import pipeline as cut

published=root/"published"
if (Path.cwd()/"pair.json").is_file():
    root=Path.cwd()
pair_file = root/'pair.json'
pair = json.loads(pair_file.read_text()) if pair_file.exists() else {}
method = pair.get('maxi_method', DEFAULT_METHOD)
source = Path(pair.get('alignment_source', root))
idx=int(sys.argv[1])
condition=pair.get('conditions', ('NoEP','HiFi'))[idx//79]
target=sorted(json.loads((CUT/'target_guide_mapping.json').read_text()),key=lambda x:x['target_id'])[idx%79]
tid=target['target_id']
work=root/condition
prep=source/condition/'prepared'
out=work/'targets'/tid
out.mkdir(parents=True)
with (prep/'targets.tsv').open() as h:
    meta=one(r for r in csv.DictReader(h,delimiter='\t') if r['target_id']==tid)
assert meta['amplicon_seq']==target['amplicon_seq'] and meta['tbc']==target['tbc']
total=int(meta['prepared_reads'])
assert int(meta['original_reads'])==total
families=defaultdict(Counter)
if total:
    reads=Path(meta['fastq_r1'])
    job=one(j for j in json.loads((prep/'crispresso_jobs.json').read_text()) if j['target_id']==tid)
    assert sha256(reads)==job['fastq_sha256']
    cmd=['CRISPResso','--fastq_r1',str(reads),'--amplicon_seq',target['amplicon_seq'],
         '--amplicon_name',tid,'--name','target','--guide_seq',target['guide_seq'],
         '--quantification_window_center','-3','--quantification_window_size','1',
         '--exclude_bp_from_left','0','--exclude_bp_from_right','0','--needleman_wunsch_gap_incentive','1',
         '--ignore_substitutions','--write_detailed_allele_table','--fastq_output','--suppress_plots',
         '--suppress_report','--output_folder',str(out),'--n_processes','1']
    if source == root:
        with (out/'alignment.log').open('x') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    results=source/condition/'targets'/tid/'CRISPResso_on_target'
    info=cut.unwrap(json.loads((results/'CRISPResso2_info.json').read_text()))
    ref=info['results']['refs'][tid]
    assert ref['include_idxs']==[51,52] and ref['sgRNA_cut_points']==[51] and ref['sgRNA_sequences']==[target['guide_seq']]
    q=one(cut.read_tsv(results/'CRISPResso_quantification_of_editing_frequency.txt'))
    aligned=int(q['Reads_aligned'])
    indels=cut.indel_count(q)
    assert int(q['Reads_in_input'])==total and int(q['Discarded'])==0 and int(q['Modified'])==indels
    count=0
    for seq,ann in itertools.zip_longest(fastq(reads),fastq(results/'CRISPResso_output.fastq.gz')):
        assert seq and ann and (seq[0],seq[1],seq[3])==(ann[0],ann[1],ann[3])
        count+=1
        fields=dict(re.findall(r'(?:^|\s)([A-Z_]+)=([^\s]*)',ann[2]))
        if fields['ALN']=='NA':
            continue
        tags=dict(v.split('=',1) for v in ann[0].split('|')[1:])
        assert fields['ALN']==tid and tags['tbc']==target['tbc'] and tags['structure']=='original'
        d,i,s=map(int,re.fullmatch(r'D(\d+);I(\d+);S(\d+)',fields['MODS']).groups())
        edited=bool(d or i)
        assert edited==(fields['CLASS']==tid+'_MODIFIED')
        families[tags['pbc']][fields['ALN_REF'],fields['ALN_SEQ'],edited]+=1
    assert count==total
else:
    aligned=indels=0

rebuilt = configured_reference(pair)
if rebuilt:
    labels = rebuilt_labels(rebuilt, tid)
elif method == DEFAULT_METHOD:
    # The published table can repeat a target/pBC. Repository step 5 takes the first row.
    published = Path(pair.get('published', published))
    assert (published/'source.sha256').read_text().strip()==PUBLISHED_SHA256
    labels={}
    with (published/(tid+'.tsv')).open() as h:
        for line in h:
            pbc,good=line.rstrip().split('\t')
            if pbc in families:
                labels.setdefault(pbc,good=='1')
else:
    # Each file holds the Good pBCs for one target in one independently evaluated reference.
    votes = Counter()
    for ref in pair['reference_dirs']:
        with (Path(ref)/(tid+'.tsv')).open() as h:
            votes.update(families.keys() & {line.rstrip() for line in h})
    labels = {pbc: n >= REQUIRED[method] for pbc, n in votes.items()}
methods=['Unfiltered','Dedup','Dedud','Dedud+Dedup']
counts={m:[0,0] for m in methods}
ties=Counter()
for alleles in families.values():
    edited, total_reads, representative, tie = family_counts(alleles)
    counts['Unfiltered'][0] += edited
    counts['Unfiltered'][1] += total_reads
    counts['Dedup'][0] += representative
    counts['Dedup'][1] += 1
    ties['Dedup'] += tie
identity_mode = pair.get('target_filter', 'pbc-only')  # legacy run configurations
identity = TargetIdentity(PUBLIC/'target_oligos_sequenes.csv') if identity_mode == 'manuscript' else None
counts['Dedud'], counts['Dedud+Dedup'], ties['Dedud+Dedup'] = filtered_counts(tid, families, labels, identity)
assert counts['Unfiltered']==[indels,aligned] and counts['Dedup'][1]==len(families)
record=dict(target_filter=identity_mode,target=tid,condition=condition,prepared=total,counts=counts,mixed_top_ties=dict(ties))
(out/'counts.json').write_text(json.dumps(record,indent=2)+'\n')
