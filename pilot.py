"""Run a bounded original-library smoke test inside a Luminary Slurm allocation."""
import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import maxi_qc

from settings import PUBLIC, REFERENCES as CUT, RUNTIME as FROZEN


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def one(rows):
    item, = rows
    return item


def fastq(path):
    with gzip.open(path, 'rt') as handle:
        for header in handle:
            yield header.rstrip(), handle.readline().rstrip(), handle.readline().rstrip(), handle.readline().rstrip()


def run(args):
    import prepare_crispresso as prep
    import barcode_core as bc
    import pipeline as cut
    from CRISPResso2 import CRISPRessoShared
    root = args.output.resolve()
    config = FROZEN / 'config'
    with args.manifest.open() as handle:
        sample = one(r for r in csv.DictReader(handle, delimiter='\t') if r['sample_name'] == args.sample and r['status'] == 'matched')
    if sample['dataset'] not in ('HiFi', 'NoEP'):
        raise ValueError('This pilot supports original LV_S-numbered HiFi and NoEP inputs.')
    lv = re.match(r'(LV_S\d+)_', Path(sample['r1']).name).group(1)
    metadata = one(r for r in table(PUBLIC / 'NovaSeq3_sample_info.csv') if r['LV_ID'] == lv)
    condition = {'noEP':'NoEP', 'HiFiCas9':'HiFi'}[metadata['cas_type']]
    if sample['sample_name'] != '_'.join([condition,'pool',metadata['donor_id'],metadata['replicate']]):
        raise ValueError('Repository and SRA sample identities disagree.')
    root.mkdir(parents=True, exist_ok=False)
    record = {'pool_id': lv}
    for mate in ('r1','r2'):
        path = sample[mate]
        if bc.sha256(path) != sample['sha256_'+mate]:
            raise ValueError('FASTQ checksum mismatch: '+path)
        record.update({mate:path, mate+'_identity':bc.identity(path)})
    manifest = root / 'inputs.json'
    manifest.write_text(json.dumps([record], indent=2))
    prepared = root / 'prepared'
    stats = prep.prepare(manifest, lv, config/'barcode_policy.json', config/'crispresso_prep_policy.json',
                         CUT/'targets.json', prepared, args.pairs)
    job = one(j for j in json.loads((prepared/'crispresso_jobs.json').read_text()) if j['target_id']==args.target)
    target = one(r for r in json.loads((CUT/'target_guide_mapping.json').read_text()) if r['target_id']==args.target)
    ref = one(r for r in table(PUBLIC/'OT_guide_amplicon_seq.csv') if r['OT-Name'].replace('-','_')==args.target)
    barcode = one(r for r in table(PUBLIC/'target_oligos_sequenes.csv') if r['OT'].replace('-','_')==args.target)
    old = job['args']
    reads = Path(old[old.index('--fastq_r1')+1])
    if (target['amplicon_seq'], target['guide_seq'], target['tbc']) != (ref['amplicon_seq'].upper(), ref['guide_seq'].upper(), barcode['Barcode']):
        raise ValueError('Published and HPC target definitions differ.')
    if target['amplicon_seq'] != old[old.index('--amplicon_seq')+1] or bc.sha256(reads)!=job['fastq_sha256']:
        raise ValueError('Prepared input or reference mismatch.')
    summary = dict(sample=sample, repository_sample=metadata, target=args.target, pairs=stats['pairs'],
                   prepared=stats['prepared_pairs'], targets_with_reads=stats['targets_with_reads'],
                   target_reads=job['prepared_reads'], version=subprocess.check_output(['CRISPResso','--version'],text=True).strip(),
                   job=os.environ['SLURM_JOB_ID'], scope='prefix pilot, raw cut-site counts', runs={})
    previous = None
    for mode in (['lean','report'] if args.compare_reports else ['lean']):
        dest = root/mode
        cmd = ['CRISPResso','--fastq_r1',str(reads),'--amplicon_seq',target['amplicon_seq'],'--amplicon_name',args.target,
               '--name',lv,'--guide_seq',target['guide_seq'],'--quantification_window_center','-3','--quantification_window_size','1',
               '--exclude_bp_from_left','0','--exclude_bp_from_right','0','--needleman_wunsch_gap_incentive','1',
               '--ignore_substitutions','--write_detailed_allele_table','--fastq_output','--output_folder',str(dest),'--n_processes','1']
        if mode == 'lean':
            cmd.extend(['--suppress_plots','--suppress_report'])
        start=time.monotonic()
        with (root/(mode+'.log')).open('x') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        elapsed=time.monotonic()-start
        folder=dest/('CRISPResso_on_'+CRISPRessoShared.slugify(lv))
        info=cut.unwrap(json.loads((folder/'CRISPResso2_info.json').read_text()))
        active=info['results']['refs'][args.target]
        assert active['include_idxs']==[51,52] and active['sgRNA_cut_points']==[51] and active['sgRNA_sequences']==[target['guide_seq']]
        quant=one(cut.read_tsv(folder/'CRISPResso_quantification_of_editing_frequency.txt'))
        counts=Counter()
        for a,b in itertools.zip_longest(fastq(reads),fastq(folder/'CRISPResso_output.fastq.gz')):
            assert a and b and (a[0],a[1],a[3])==(b[0],b[1],b[3])
            tags=dict(v.split('=',1) for v in a[0].split('|')[1:])
            ann=dict(re.findall(r'(?:^|\s)([A-Z_]+)=([^\s]*)',b[2]))
            assert tags['tbc']==target['tbc'] and tags['structure']=='original'
            counts[tags['pbc'],ann['ALN'],ann.get('CLASS',''),ann.get('ALN_REF',''),ann.get('ALN_SEQ','')]+=1
        aligned=sum(n for k,n in counts.items() if k[1]==args.target)
        assert aligned==int(quant['Reads_aligned']) and aligned>0
        assert sum(counts.values())==job['prepared_reads']==int(quant['Reads_in_input'])
        if previous is not None and previous!=(quant,counts):
            raise ValueError('Reporting modes disagree on counts or pBC-resolved alignments.')
        previous=quant,counts
        if mode=='lean':
            with (root/'pbc_alleles.tsv').open('x') as handle:
                writer=csv.writer(handle,delimiter='\t');writer.writerow(['pbc','alignment','class','aligned_reference','aligned_read','reads'])
                for key,n in sorted(counts.items()):
                    writer.writerow([*key,n])
        summary['runs'][mode]=dict(seconds=elapsed,aligned=aligned,cutsite_indels=cut.indel_count(quant),
                                   pbc_count=len({k[0] for k in counts}),command=cmd)
    summary['source_sha256']={str(p):bc.sha256(p) for p in (Path(__file__),args.manifest,Path(prep.__file__),Path(bc.__file__),
        Path(cut.__file__),CUT/'target_guide_mapping.json',CUT/'targets.json',config/'barcode_policy.json',config/'crispresso_prep_policy.json',
        PUBLIC/'OT_guide_amplicon_seq.csv',PUBLIC/'NovaSeq3_sample_info.csv',PUBLIC/'target_oligos_sequenes.csv')}
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--sample',default='NoEP_pool_D3734_T1')
    parser.add_argument('--target',default='1617_OT_0041_ALT')
    parser.add_argument('--pairs',type=int,default=50000)
    parser.add_argument('--compare-reports',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    maxi_qc.add_options(parser)
    parser.add_argument('--skip-maxi-qc', action='store_true', help='Run the raw alignment pilot only.')
    args=parser.parse_args()
    if args.skip_maxi_qc and args.maxi_method:
        parser.error('--skip-maxi-qc conflicts with --maxi-method.')
    methods = args.maxi_method or [maxi_qc.DEFAULT_METHOD]
    if args.min_exact_reads < 1 or (not any(m in maxi_qc.REQUIRED for m in methods) and args.min_exact_reads != 1):
        parser.error("A custom positive exact-read threshold requires --maxi-method.")
    if not os.environ.get('SLURM_JOB_ID') or not 1<=args.pairs<=250000:
        parser.error('A Slurm allocation and a prefix between 1 and 250000 pairs are required.')
    sys.path[:0]=[str(FROZEN/'src')]
    run(args)
    if not args.skip_maxi_qc:
        maxi_qc.classify(args.output, methods, args.min_exact_reads, args.output/"maxi_qc")
