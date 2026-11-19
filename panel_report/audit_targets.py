"""Trace WT read- and representative-family support and replay the repository target-identity check."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory
from trace import scan, events, write

def distance(a,b):
    row=list(range(len(b)+1))
    for i,x in enumerate(a,1):
        prev=row;row=[i]
        for j,y in enumerate(b,1):
            row.append(min(prev[j]+1,row[-1]+1,prev[j-1]+(x!=y)))
    return row[-1]

def label(p):
    return 'pBC_'+hashlib.sha256(p.encode()).hexdigest()[:16]

def summarize(families):
    return dict(reads=sum(v.total() for v in families.values()),
                indel_reads=sum(n for v in families.values() for a,n in v.items() if a[2]),
                families=len(families),any_indel_families=sum(any(a[2] for a in v) for v in families.values()),
                indel_families=sum(max(v,key=v.get)[2] for v in families.values()))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cohort',type=Path);p.add_argument('panel_reference',type=Path);p.add_argument('out',type=Path)
    p.add_argument('targets',nargs='+')
    args=p.parse_args()
    history=json.loads((args.cohort/'results/provenance.json').read_text())['runs']
    assert len(history)==8
    with args.panel_reference.open() as f:
        panel={r['OT'].replace('-','_'):r['OT_sequence'].upper() for r in csv.DictReader(f)}
    assert len(panel)==79 and set(args.targets)<=panel.keys()
    counts=[];presence=[];alleles=[];identities=[];qc=[]
    @lru_cache(maxsize=10000)
    def match(target,seq):
        scored={k:distance(seq,v) for k,v in panel.items()}
        best=min(scored.values());ids=sorted(k for k,v in scored.items() if v==best)
        outcome='closer_other' if scored[target]!=best else ('tied' if len(ids)>1 else 'unique_expected')
        return outcome,','.join(ids),scored[target],best
    for t in args.targets:
        pub=Path(history[0]['published']);assert (pub/'source.sha256').read_text().strip()==history[0]['published_sha256']
        good={}
        with (pub/(t+'.tsv')).open() as f:
            for line in f:
                k,x=line.rstrip().split('\t');good.setdefault(k,x=='1')
        saved={};selected=set()
        for cfg in history:
            d,r=cfg['donor'],cfg['replicate']
            for c in ('WT','NoEP','HiFi'):
                folder=(args.cohort/'runs'/(d+'_'+r) if c=='WT' else Path(cfg['previous_pair']))/c/'targets'/t
                allfam,record=scan(folder,t)
                assert record['condition']==c
                fam={pbc:v for pbc,v in allfam.items() if good.get(pbc)}
                saved[d,r,c]=fam
                base=dict(target=t,donor=d,replicate=r,condition=c)
                totals=summarize(fam)
                assert [totals['indel_reads'],totals['reads']]==record['counts']['Dedud']
                assert [totals['indel_families'],totals['families']]==record['counts']['Dedud+Dedup']
                counts.append(dict(**base,rule='published_pBC',**totals))
                if c=='WT':selected.update(k for k,v in fam.items() if any(a[2] for a in v))
                for mode in ('step6_as_aligned','step6_target_only'):
                    filtered={}
                    scores=Counter()
                    for k,v in fam.items():
                        keep=Counter()
                        for a,n in v.items():
                            ref,read,e=a
                            if not e:
                                keep[a]=n
                                continue
                            ref=ref.replace('-','');seq=read.replace('-','')
                            assert ref[32:-4]==panel[t] and len(ref)==68
                            # Current prepared R1 ends in four adapter bases; step 6 slices the ungapped allele at 32.
                            assert seq.endswith(ref[-4:])
                            q=seq[32:] if mode=='step6_as_aligned' else seq[32:-4]
                            outcome,best,de,db=match(t,q)
                            scores[outcome,best,de,db]+=n
                            if outcome=='unique_expected':keep[a]=n
                        if keep:filtered[k]=keep
                    counts.append(dict(**base,rule=mode,**summarize(filtered)))
                    for (outcome,best,de,db),n in scores.items():
                        identities.append(dict(**base,rule=mode,outcome=outcome,nearest_targets=best,
                            distance_expected=de,distance_nearest=db,indel_reads=n))
        for (d,r,c),fam in sorted(saved.items()):
            for k in sorted(selected):
                v=fam.get(k,Counter())
                top=max(v,key=v.get) if v else None
                base=dict(target=t,donor=d,replicate=r,condition=c,family=label(k))
                presence.append(dict(**base,aligned_reads=v.total(),indel_reads=sum(n for a,n in v.items() if a[2]),
                       distinct_alleles=len(v),top_support=v[top] if v else 0,
                       top_ties=sum(n==v[top] for n in v.values()) if top else 0,representative_indel=int(top[2]) if top else '',
                       event=events(*top[:2])[0] if top else '',event_id=events(*top[:2])[1] if top else ''))
                for a,n in v.items():
                    ref,read,e=a;kind,event=events(ref,read)
                    seq=read.replace('-','')
                    res=match(t,seq[32:-4]) if e else ('unfiltered_unedited','',None,None)
                    literal=match(t,seq[32:]) if e else ('unfiltered_unedited','',None,None)
                    alleles.append(dict(**base,allele_id=hashlib.sha256((ref+'|'+read).encode()).hexdigest()[:16],
                           reads=n,indel=int(e),representative=int(a==top),event=kind,event_id=event,
                           identity_target_only=res[0],nearest_target_only=res[1],identity_as_aligned=literal[0],nearest_as_aligned=literal[1]))
        qc.append(dict(target=t,WT_distinct_pBCs_any_edit=len(selected),libraries=len(saved)))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=args.out.parent) as tmp:
        out=Path(tmp)/'audit';out.mkdir()
        for name,rows in [('sample_counts',counts),('family_presence',presence),('family_alleles',alleles),
                          ('target_identity',identities),('trace_summary',qc)]:write(out/(name+'.tsv'),rows)
        out.rename(args.out)

if __name__=='__main__':main()
