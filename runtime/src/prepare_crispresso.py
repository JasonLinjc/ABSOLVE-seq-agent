#!/usr/bin/env python3
"""Prepare single-end CRISPResso2 inputs while retaining original pBC metadata."""
import argparse,collections,csv,gzip,json,os,re,shlex,time
from pathlib import Path
import barcode_core as b

VERSION='1.0.0'


def transform(one,two,recognizer,policy,number,sample):
    h,s,q=one
    found,reason=recognizer.recognize(s)
    if found is None:return None,reason
    name,prefix,y,pbc,tbc=found
    bp=recognizer.p
    if set(pbc)-set('ACGT'):return None,'ambiguous_pbc'
    if not b.quality_ok(q[len(prefix):y],bp['min_pbc_base_q'],bp['min_pbc_mean_q']):return None,'pbc_low_quality'
    if not b.quality_ok(q[y+23:y+32],bp['min_tbc_base_q'],bp['min_tbc_mean_q']):return None,'tbc_low_quality'
    start=y+32;adapter=bp['adapter'];hit=s.find(adapter,start)
    if hit<0:return None,'missing_exact_adapter'
    if s.find(adapter,hit+1)>=0:return None,'ambiguous_adapter_boundary'
    target=s[start:hit];tq=q[start:hit]
    if len(target)>policy['max_target_length']:return None,'target_length_outside_policy'
    if set(target)-set('ACGT'):return None,'ambiguous_target_bases'
    if target and not b.quality_ok(tq,policy['min_target_base_q'],policy['min_target_mean_q']):return None,'target_low_quality'
    # Keep four observed adapter bases as non-biological right alignment context.
    # They are excluded from the reference's quantification window.
    sequence=s[y:hit+4];qual=q[y:hit+4]
    if len(sequence)!=len(qual):raise AssertionError('Trimmed sequence/quality mismatch')
    original=h.split()[0][1:]
    if any(x in original for x in ['|','\t']):raise ValueError('Reserved separator in original FASTQ ID')
    header=f'@{original}|sample={sample}|record={number}|pbc={pbc}|tbc={tbc}|structure={name}'
    r2_tested=int(bool(target) and len(two[1])>=len(target))
    r2_agree=int(bool(r2_tested) and two[1][:len(target)]==b.rc(target))
    return {'header':header,'sequence':sequence,'quality':qual,'pbc':pbc,'tbc':tbc,'structure':name,
            'target_length':len(target),'r2_tested':r2_tested,'r2_agree':r2_agree},None


def prepare(manifest,sample_id,barcode_policy,policy_path,targets_path,output,max_pairs=None):
    policy=json.loads(Path(policy_path).read_text());bp,targets=b.load_config(barcode_policy,targets_path)
    if max_pairs is not None and max_pairs<1:raise ValueError('max_pairs must be positive')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',sample_id):raise ValueError('Unsafe sample ID')
    for t in targets.values():
        if not re.fullmatch(r'[A-Za-z0-9_-]+',t['id']):raise ValueError('Unsafe target ID')
    manifest_rows=json.loads(Path(manifest).read_text());matches=[s for s in manifest_rows if s['pool_id']==sample_id]
    if len(matches)!=1 or len({s['pool_id'] for s in manifest_rows})!=len(manifest_rows):raise ValueError('Sample IDs not unique')
    sample=matches[0];inputs=[b.identity(sample[f'r{m}']) for m in (1,2)]
    if inputs[0]['path']==inputs[1]['path']:raise ValueError('R1 and R2 must differ')
    for m,identity in zip((1,2),inputs):
        if sample.get(f'r{m}_identity')!=identity:raise ValueError('Input differs from frozen manifest')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'fastq').mkdir()
    rec=b.Recognizer(bp,targets);handles={};counts=collections.defaultdict(collections.Counter)
    dispositions=collections.Counter();structures=collections.Counter();lengths=collections.Counter();examples={};n=0
    started=time.monotonic()
    run={'version':VERSION,'sample_id':sample_id,'sample':sample,'inputs':inputs,'scope':'pilot' if max_pairs else 'full',
         'max_pairs':max_pairs,'policy':policy,'barcode_policy':bp,'targets':list(targets.values()),
         'started_utc':b.utcnow(),'slurm_job_id':os.getenv('SLURM_JOB_ID'),
         'source_sha256':{str(p):b.sha256(p) for p in [Path(__file__),Path(b.__file__)]},
         'input_config_sha256':{str(p):b.sha256(p) for p in [manifest,barcode_policy,policy_path,targets_path]},
         'layout':'constant23 + tBC9 + variable-length target + observed AGAT4',
         'reference_target_coordinates':'32-63, zero-based inclusive; whole target, not a nuclease-specific window',
         'pbc_filter':'structure/quality only; no plasmid-pool whitelist applied; no barcode correction or deduplication',
         'read_metadata':'Original first-token cluster ID followed by sample,record,pbc,tbc,structure in FASTQ header',
         'status':'running'}
    try:
        for one,two in b.pairs(sample['r1'],sample['r2'],max_pairs):
            n+=1;r,reason=transform(one,two,rec,policy,n,sample_id)
            dispositions[reason or 'prepared']+=1
            if reason:
                examples.setdefault(reason,{'record':n,'header':one[0],'r1_prefix':one[1][:80]});continue
            tbc=r['tbc'];target=targets[tbc];tid=target['id']
            if tbc not in handles:
                handles[tbc]=gzip.open(out/'fastq'/(tid+'.R1.fastq.gz.partial'),'wt',compresslevel=1,newline='')
            handles[tbc].write(r['header']+'\n'+r['sequence']+'\n+\n'+r['quality']+'\n')
            counts[tbc]['reads']+=1;counts[tbc]['r2_target_tested']+=r['r2_tested'];counts[tbc]['r2_target_agree']+=r['r2_agree']
            counts[tbc]['structure_'+r['structure']]+=1;structures[r['structure']]+=1;lengths[r['target_length']]+=1
        for handle in handles.values():handle.close()
        if n==0:raise ValueError('No input pairs')
        if max_pairs is not None and n!=max_pairs:raise ValueError('Input shorter than requested pilot')
        if max_pairs is None and sample.get('expected_pairs') is not None and n!=sample['expected_pairs']:raise ValueError('Read count differs from completed QC')
        if inputs!=[b.identity(sample[f'r{m}']) for m in (1,2)]:raise ValueError('Inputs changed during preparation')
        if sum(dispositions.values())!=n or sum(x['reads'] for x in counts.values())!=dispositions['prepared']:raise AssertionError('Read conservation failed')
        jobs=[];file_digests={}
        with open(out/'targets.tsv','x',newline='') as f,open(out/'amplicons.fasta','x') as fasta:
            w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','target_id','tbc','amplicon_seq','quantification_window_coordinates','prepared_reads','r2_target_tested','r2_target_agree','original_reads','CACC_reads','STUFFER_CACC_reads','status','fastq_r1'])
            for tbc,t in targets.items():
                amp=bp['constant']+tbc+t['target']+bp['adapter'][:4]
                assert len(amp)==68 and amp[32:64]==t['target']
                fasta.write('>'+t['id']+'\n'+amp+'\n')
                stat=counts[tbc];path=out/'fastq'/(t['id']+'.R1.fastq.gz')
                if stat['reads']:
                    partial=Path(str(path)+'.partial');os.link(partial,path);partial.unlink();file_digests[str(path.relative_to(out))]=b.sha256(path)
                    name=sample_id+'__'+t['id']
                    args=['--fastq_r1',str(path),'--amplicon_seq',amp,'--amplicon_name',t['id'],'--name',name,
                          '--quantification_window_coordinates','32-63','--exclude_bp_from_left','0','--exclude_bp_from_right','0',
                          '--needleman_wunsch_gap_incentive','0','--write_detailed_allele_table','--fastq_output',
                          '--output_folder',str(out/'crispresso_results')]
                    jobs.append({'sample_id':sample_id,'target_id':t['id'],'scope':run['scope'],'prepared_reads':stat['reads'],'fastq_sha256':file_digests[str(path.relative_to(out))],'args':args})
                w.writerow([sample_id,t['id'],tbc,amp,'32-63',stat['reads'],stat['r2_target_tested'],stat['r2_target_agree'],
                            stat['structure_original'],stat['structure_CACC'],stat['structure_STUFFER_CACC'],
                            'ready' if stat['reads'] else 'no_eligible_reads',str(path) if stat['reads'] else ''])
        (out/'crispresso_jobs.json').write_text(json.dumps(jobs,indent=2))
        (out/'output_sha256.json').write_text(json.dumps(file_digests,indent=2))
        run.update({'status':'complete','counts_conserved':True,'pairs':n,'prepared_pairs':dispositions['prepared'],
                    'dispositions':dict(dispositions),'structures':dict(structures),'target_length_counts':dict(lengths),
                    'targets_with_reads':len(jobs),'excluded_examples':examples,'finished_utc':b.utcnow(),
                    'elapsed_seconds':time.monotonic()-started,'crispresso_execution':'not_run'})
        temp=out/'summary.json.partial';temp.write_text(json.dumps(run,indent=2));os.link(temp,out/'summary.json');temp.unlink()
        return run
    except BaseException as exc:
        for handle in handles.values():
            if not handle.closed:handle.close()
        (out/'FAILED.json').write_text(json.dumps({'status':'failed','pairs_seen':n,'error':str(exc)}));raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['manifest','sample','barcode-policy','policy','targets','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--max-pairs',type=int);a=p.parse_args()
    r=prepare(a.manifest,a.sample,a.barcode_policy,a.policy,a.targets,a.output,a.max_pairs)
    print(json.dumps({k:r[k] for k in ['sample_id','status','scope','pairs','prepared_pairs','targets_with_reads','dispositions']}))

if __name__=='__main__':main()
