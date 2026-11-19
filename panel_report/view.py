"""Make the matched-panel evidence tables browsable in a self-contained offline report."""
import argparse
import csv
import json
from pathlib import Path

def tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def create(folder):
    data={'panel':tsv(folder/'panel_evidence.tsv'),'donor':tsv(folder/'donor_evidence.tsv'),'trace':[]}
    trace=folder/'family_trace/family_presence.tsv'
    if trace.exists():data['trace']=tsv(trace)
    payload=json.dumps(data,separators=(',',':')).replace('<',r'\u003c')
    html=Path(__file__).with_name('view.html').read_text().replace('/*DATA*/',payload)
    (folder/'index.html').write_text(html)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path)
    create(p.parse_args().folder)
