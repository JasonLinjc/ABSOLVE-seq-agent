"""Read CRISPResso2 result tables and metadata."""
import csv
from pathlib import Path


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def unwrap(x):
    if isinstance(x, dict) and '_type' in x and 'value' in x:
        return unwrap(x['value'])
    if isinstance(x, dict):
        return {k:unwrap(v) for k,v in x.items()}
    if isinstance(x, list):
        return [unwrap(v) for v in x]
    return x


def indel_count(quant):
    return sum(int(quant[k]) for k in ['Only Insertions','Only Deletions',
        'Insertions and Deletions','Insertions and Substitutions',
        'Deletions and Substitutions','Insertions Deletions and Substitutions'])
