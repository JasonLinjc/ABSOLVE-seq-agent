"""Validate complete target coverage and summarize paired rates as percentages."""
import csv
import json
import statistics
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

root=Path(__file__).resolve().parent
pair={'donor':'D3734','replicate':'T1'}
if (Path.cwd()/'pair.json').is_file():
    root=Path.cwd()
    pair=json.loads((root/'pair.json').read_text())
methods=('Unfiltered' ,'Dedup','Dedud','Dedud+Dedup')
conditions=('NoEP','HiFi')
def write(name,rows):
    with (out/name).open('x',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
def rate(n,d):
    return 100*n/d if d else None
def difference(a,b):
    return a-b if a is not None and b is not None else None

counts={}
qc={}
for c in conditions:
    source = Path(pair.get('alignment_source', root))
    prep=json.loads((source/c/'prepared/summary.json').read_text())
    assert prep['status']=='complete' and prep['scope']=='full' and prep['counts_conserved']
    counts[c]={}
    for path in sorted((root/c/'targets').glob('*/counts.json')):
        j=json.loads(path.read_text())
        assert j['target']==path.parent.name and j['condition']==c and set(j['counts'])==set(methods)
        for n,d in j['counts'].values():
            assert 0<=n<=d
        for a,b in [('Dedud','Unfiltered'),('Dedud+Dedup','Dedup'),('Dedup','Unfiltered'),('Dedud+Dedup','Dedud')]:
            assert all(x<=y for x,y in zip(j['counts'][a],j['counts'][b]))
        counts[c][j['target']]=j
    assert len(counts[c])==79 and sum(j['prepared'] for j in counts[c].values())==prep['prepared_pairs']
    qc[c]={k:prep[k] for k in ('sample_id','pairs','prepared_pairs','dispositions','structures','input_config_sha256','source_sha256')}
assert counts['NoEP'].keys()==counts['HiFi'].keys()
rows=[]
for t in sorted(counts['NoEP']):
    for m in methods:
        row=dict(donor=pair['donor'],replicate=pair['replicate'],target=t,strategy=m)
        if 'maxi_method' in pair:
            row['maxi_method'] = pair['maxi_method']
            row['min_exact_reads'] = pair['min_exact_reads']
            row['target_filter'] = pair.get('target_filter', 'pbc-only')
        for c in conditions:
            j=counts[c][t]
            n,d=j['counts'][m]
            row.update({c+'_indel':n,c+'_denominator':d,c+'_rate_pct':rate(n,d),
                        c+'_mixed_top_ties':j['mixed_top_ties'].get(m,0) if 'Dedup' in m else None})
        row['delta_HiFi_minus_NoEP_pp']=difference(row['HiFi_rate_pct'],row['NoEP_rate_pct'])
        rows.append(row)
summary=[]
for m in methods:
    subset=[r for r in rows if r['strategy']==m]
    delta=[r['delta_HiFi_minus_NoEP_pp'] for r in subset if r['delta_HiFi_minus_NoEP_pp'] is not None]
    row=dict(strategy=m,targets=len(subset),paired_rates=len(delta))
    for c in conditions:
        n=sum(r[c+'_indel'] for r in subset)
        d=sum(r[c+'_denominator'] for r in subset)
        row.update({c+'_indel':n,c+'_denominator':d,c+'_pooled_rate_pct':rate(n,d)})
    row['pooled_delta_pp']=difference(row['HiFi_pooled_rate_pct'],row['NoEP_pooled_rate_pct'])
    row['median_target_delta_pp']=statistics.median(delta) if delta else None
    summary.append(row)
with TemporaryDirectory(dir=root) as temporary:
    out=Path(temporary)/'results'
    out.mkdir()
    write('paired_target_rates.tsv',rows)
    write('strategy_summary.tsv',summary)
    (out/'qc.json').write_text(json.dumps(qc,indent=2)+'\n')
    out.rename(root/'results')
print(len(rows),'sample-paired target/strategy rows')
