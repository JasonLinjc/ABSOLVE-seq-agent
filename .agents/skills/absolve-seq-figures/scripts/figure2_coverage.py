"""Export Figure 2a coverage from the archived plot_results.ipynb allele-table layout."""
import argparse
import csv
from itertools import groupby
import sys
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    counts = defaultdict(lambda: [0, 0])
    folders = [('target_umi_barcode_table_WTCas9_NoEP_filtered', {'WTCas9', 'NoEP'}),
               ('target_umi_barcode_table_HiFi_NoEP', {'HiFi'})]
    for folder, groups in folders:
        for method in ['raw', 'dedudS_filtered_v3']:
            files = list((args.archive/folder/method).glob('*annot.txt'))
            if not files:
                raise ValueError('Missing archived allele tables in ' + str(args.archive/folder/method))
            def target_id(path):
                return '_'.join(path.name.split('_')[:4])
            for name, paths in groupby(sorted(files), key=target_id):
                unique_barcodes = defaultdict(set)
                for path in paths:
                    if not any('_'+group+'_' in path.name for group in groups):
                        continue
                    with path.open() as stream:
                        it = csv.reader(stream, delimiter='\t')
                        columns = next(it)
                        target, arm, umi, reads = [columns.index(k) for k in ['OT_name','group','UMI','#Reads']]
                        for row in it:
                            if row[arm] in groups:
                                if row[target].replace('-', '_') != name:
                                    raise ValueError('Target mismatch in ' + str(path))
                                cell = counts[method,name,row[arm]]
                                cell[0] += int(row[reads])
                                unique_barcodes[row[arm]].add(row[umi])
                for arm, barcodes in unique_barcodes.items():
                    counts[method,name,arm][1] = len(barcodes)
    if len(counts) != 2*79*3:
        raise ValueError('Incomplete historical coverage table')
    writer=csv.writer(sys.stdout,delimiter='\t',lineterminator='\n')
    writer.writerow(['target','arm','strategy','coverage'])
    for (method,target,arm),(reads,observations) in sorted(counts.items()):
        labels=('Unfiltered','Dedup') if method=='raw' else ('Dedud','Dedud+Dedup')
        writer.writerow([target,arm,labels[0],reads])
        writer.writerow([target,arm,labels[1],observations])


if __name__=='__main__':
    main()
