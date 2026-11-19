"""Extract the named panel-comparison columns from the four manuscript workbooks."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from openpyxl import load_workbook

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('destination',type=Path)
    parser.add_argument('--full-model', action='store_true', help='Include all submitted model columns for manuscript reproduction.')
    args=parser.parse_args()
    hashes={}
    args.destination.mkdir(parents=True,exist_ok=True)
    fields=['barcode','predicted-NoEP','predicted-HiFi','predicted-WTCas9','asymPval-treatHiFi','asymPval-treatWTCas9']
    for num in range(2,6):
        name=f'SupplementaryTable{num}'
        wb=load_workbook(args.source/(name+'.xlsx'),read_only=True,data_only=True)
        data=iter(wb.active.values);next(data);head=next(data)
        if args.full_model:
            fields=[k for k in head if k=='barcode' or k.startswith(('estim-','zScore-','asymPval-','predicted-','stdError-','confHigh-','confLow-'))]
        indices=[head.index(k) for k in fields]
        hashes[name+'.xlsx']=hashlib.sha256((args.source/(name+'.xlsx')).read_bytes()).hexdigest()
        rows=[[row[i] for i in indices] for row in data if row and row[0]]
        wb.close()
        assert len(rows)==len({r[0] for r in rows})==79
        assert all(0<=x<=1 for row in rows for k,x in zip(fields,row) if k.startswith(('predicted-','asymPval-','confHigh-','confLow-')))
        with (args.destination/(name+'.csv')).open('x',newline='') as f:
            w=csv.writer(f);w.writerow(fields);w.writerows(rows)
    if args.full_model:
        (args.destination/'source_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
if __name__=='__main__':main()
