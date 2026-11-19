"""Join SRA submission rows to verified HPC mates; retain unresolved pairs."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from openpyxl import load_workbook


def resolve(inventory, workbook, destination):
    sheet = load_workbook(workbook, read_only=True, data_only=True)['SRA_data']
    values = iter(sheet.values)
    header = next(values)
    by_pair = defaultdict(list)
    for number, values_row in enumerate(values, 2):
        row = dict(zip(header, values_row))
        if row['filename'] and row['filename2']:
            by_pair[(row['filename'], row['filename2'])].append((number, row))
    groups = defaultdict(dict)
    with open(inventory) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['hpc_status'] != 'sha256_verified':
                raise ValueError('Inventory contains an unverified file')
            group = groups[row['dataset'], row['pair_key']]
            if row['mate'] in group:
                raise ValueError('Duplicate mate in inventory')
            group[row['mate']] = row
    with open(destination, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['dataset','pair_key','r1','r2','sha256_r1','sha256_r2','sra_row','sample_name','library_ID','status'])
        for (dataset, key), mates in sorted(groups.items()):
            r1, r2 = (mates[m]['hpc_fastq_path'] for m in ('R1','R2'))
            hits = by_pair[Path(r1).name, Path(r2).name]
            number, row = hits[0] if len(hits) == 1 else ('', {})
            status = 'matched' if len(hits) == 1 and row.get('sample_name') else 'unresolved'
            writer.writerow([dataset,key,r1,r2,mates['R1']['sha256'],mates['R2']['sha256'],number,
                             row.get('sample_name',''),row.get('library_ID',''),status])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ['inventory','workbook','output']:
        parser.add_argument(arg)
    args = parser.parse_args()
    resolve(args.inventory, args.workbook, args.output)
