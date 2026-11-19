"""Audit representative alleles for Dedud indel families using archived CRISPResso annotations."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory

def write(path,rows):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

@lru_cache(maxsize=8192)
def events(ref,read):
    # Zero-based reference coordinates; boundary is after reference position 51.
    pos=0; i=0; changes=[]; payloads=[]
    while i<len(ref):
        if ref[i]=='-':
            j=i+1
            while j<len(ref) and ref[j]=='-':j+=1
            if pos-1 in (51,52) and pos in (51,52):
                changes.append(f'ins{j-i}@{pos-1}|{pos}');payloads.append(read[i:j])
        elif read[i]=='-':
            j=i+1
            while j<len(ref) and read[j]=='-':j+=1
            if pos<=52 and pos+j-i-1>=51:changes.append(f'del{j-i}@{pos}..{pos+j-i-1}')
            pos+=j-i
        else:
            j=i+1;pos+=1
        i=j
    kind=','.join(changes) or 'none'
    return kind,hashlib.sha256((kind+'|'+','.join(payloads)).encode()).hexdigest()[:16] if changes else ''

def scan(folder,target,allowed=None):
    counts=json.loads((folder/'counts.json').read_text())
    assert counts['target']==target
    families=defaultdict(Counter)
    if counts['prepared']:
        with gzip.open(folder/'CRISPResso_on_target/CRISPResso_output.fastq.gz','rt') as f:
            for header in f:
                f.readline();ann=f.readline().strip();f.readline()
                fields=dict(re.findall(r'(?:^|\s)([A-Z_]+)=([^\s]*)',ann))
                if fields['ALN']=='NA':continue
                assert fields['ALN']==target
                pbc=dict(s.split('=',1) for s in header.rstrip().split('|')[1:])['pbc']
                if allowed is not None and pbc not in allowed:continue
                d,i,_=map(int,re.fullmatch(r'D(\d+);I(\d+);S(\d+)',fields['MODS']).groups())
                edited=bool(d or i)
                families[pbc][(fields['ALN_REF'],fields['ALN_SEQ'],edited)]+=1
    return families,counts

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cohort',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--target',default='1450_OT_0071_REF')
    args=p.parse_args()
    cohort=args.cohort;target=args.target
    history=json.loads((cohort/'results/provenance.json').read_text())['runs']
    assert len(history)==8
    pub=Path(history[0]['published']); assert (pub/'source.sha256').read_text().strip()==history[0]['published_sha256']
    good={}
    with (pub/(target+'.tsv')).open() as f:
        for l in f:
            key,value=l.rstrip().split('\t');good.setdefault(key,value=='1')
    locations={}
    selected=set()
    for conf in history:
        d,r=conf['donor'],conf['replicate']
        for arm in ('NoEP','HiFi','WT'):
            src=cohort/'runs'/(d+'_'+r) if arm=='WT' else Path(conf['previous_pair'])
            folder=src/arm/'targets'/target
            locations[d,r,arm]=folder
            if r=='T2' and arm=='WT' and d in ('D3734','D3742'):
                families,counts=scan(folder,target)
                retained={k:v for k,v in families.items() if good.get(k)}
                edited={k for k,v in retained.items() if max(v,key=v.get)[2]}
                assert len(edited)==counts['counts']['Dedud+Dedup'][0] and len(retained)==counts['counts']['Dedud+Dedup'][1]
                assert sum(n for v in retained.values() for (_,_,e),n in v.items() if e)==counts['counts']['Dedud'][0]
                selected.update(edited)
    assert selected
    presence=[];alleles=[]
    for (d,r,arm),folder in sorted(locations.items()):
        families,counts=scan(folder,target,selected)
        for pbc in sorted(selected):
            label='pBC_'+hashlib.sha256(pbc.encode()).hexdigest()[:16]
            v=families.get(pbc,Counter());top=max(v,key=v.get) if v else None
            keys=dict(donor=d,replicate=r,condition=arm,target=target,family=label)
            n=sum(v.values());k=sum(n for (_,_,e),n in v.items() if e)
            presence.append(dict(**keys,aligned_reads=n,cut_indel_reads=k,distinct_aligned_alleles=len(v),
                top_support=v[top] if top else 0,top_ties=sum(a==v[top] for a in v.values()) if top else 0,
                top_cut_indel=int(top[2]) if top else '',top_event=events(*top[:2])[0] if top else '',top_event_id=events(*top[:2])[1] if top else '',
                sample_dedud_indels=counts['counts']['Dedud'][0],sample_dedud_families=counts['counts']['Dedud+Dedup'][1]))
            for (ref,read,e),num in v.most_common():
                event,event_id=events(ref,read)
                alleles.append(dict(**keys,allele=hashlib.sha256((ref+'|'+read).encode()).hexdigest()[:16],reads=num,
                      cut_indel=int(e),window_event=event,window_event_id=event_id,is_representative=(ref,read,e)==top))
                if (ref,read,e)==top and e: assert event!='none'
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as tmp:
        dest=Path(tmp)/'trace';dest.mkdir()
        write(dest/'family_presence.tsv',presence);write(dest/'family_alleles.tsv',alleles)
        dest.rename(args.output)
if __name__=='__main__':main()
