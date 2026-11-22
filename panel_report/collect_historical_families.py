"""Collect one historical per-pBC CRISPRessoBatch into its two allele tables."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import BadZipFile

import pandas as pd


WINDOW = '_Allele_frequency_table_withUMI.txt'
FULL = '_Allele_frequency_table_full_withUMI.txt'


def completed_families(batch):
    path = batch / 'CRISPResso2Batch_info.json'
    names = json.loads(path.read_text())['results']['completed_batch_arr']
    if not isinstance(names, list) or not names or not all(isinstance(n, str) and n.strip() for n in names):
        raise ValueError(f'{path}: expected a nonempty completed_batch_arr of pBC names')
    if len(names) != len(set(names)):
        raise ValueError(f'{path}: duplicate completed pBC')
    return set(names)


def check_membership(expected, tables):
    differences = []
    for name, table in tables.items():
        observed = set(table['UMI'])
        for label, ids in [('missing', expected - observed), ('unexpected', observed - expected)]:
            if ids:
                differences.append(f'{name}: {label} pBC: {", ".join(sorted(ids))}')
    if differences:
        raise ValueError('Historical pBC collection failed before counting.\n' + '\n'.join(differences))


def read_alleles(path):
    try:
        return pd.read_csv(path, sep='\t')
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def collect(batch, library_id, group, output):
    if output.exists():
        raise ValueError(f'Choose a new output directory: {output}')
    sample_id = batch.name.removeprefix('CRISPRessoBatch_on_').removesuffix('_input')
    if batch.name != f'CRISPRessoBatch_on_{sample_id}_input':
        raise ValueError(f'Expected a CRISPRessoBatch_on_<sample>_input directory: {batch}')
    _, donor, replicate, target = sample_id.split('_', 3)
    library = library_id.replace('_', '-')
    prefix = f'{target}_{library}_{donor}_{group}_{replicate}'
    if Path(prefix).name != prefix:
        raise ValueError('Historical table names cannot contain a directory separator')
    window_name, full_name = prefix + WINDOW, prefix + FULL
    frames = {window_name: [], full_name: []}
    expected = completed_families(batch)
    # Preserve the archived collector's window-driven discovery and row order.
    seen = set()
    for window in batch.glob('CRISPResso_on_*/Alleles_frequency_table_around_sgRNA*'):
        family = window.parent
        pbc = family.name.removeprefix('CRISPResso_on_')
        if pbc in seen:
            raise ValueError(f'{window_name}: multiple guide-window files for pBC {pbc}')
        seen.add(pbc)
        guide = read_alleles(window)
        guide['UMI'] = pbc
        frames[window_name].append(guide)
        archive = family / 'Alleles_frequency_table.zip'
        if archive.is_file():
            full = read_alleles(archive)
            info = json.loads((family / 'CRISPResso2_info.json').read_text())
            full['UMI'] = pbc
            full['sgRNA_cut_points'] = info['results']['refs']['Reference']['sgRNA_cut_points'][0]
            frames[full_name].append(full)

    tables = {name: pd.concat(parts).reset_index() if parts else pd.DataFrame(columns=['UMI'])
              for name, parts in frames.items()}
    check_membership(expected, tables)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Publish the two completed files together, including when a write fails partway.
    with TemporaryDirectory(dir=output.parent) as temp:
        staged = Path(temp) / 'tables'
        staged.mkdir()
        for name, table in tables.items():
            table.rename(columns={'%Reads': '%Reads_UMI'}, inplace=True)
            table['%Reads'] = table['#Reads'] / table['#Reads'].sum() * 100
            table['sample_id'] = sample_id
            table['donor'] = donor
            table['replicate'] = replicate
            table['group'] = group
            table['OT_name'] = target
            table['LV_id'] = library
            table.to_csv(staged / name, sep='\t', index=False)
        staged.rename(output)
    return len(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True, help='Original CRISPRessoBatch_on_<sample>_input directory')
    parser.add_argument('--library-id', required=True, help='Historical LV identifier, e.g. LV-S24')
    parser.add_argument('--group', required=True, help='Group label from the historical metadata, e.g. WTCas9')
    parser.add_argument('--output', type=Path, required=True, help='New directory for both verified allele tables')
    args = parser.parse_args()
    try:
        count = collect(args.batch, args.library_id, args.group, args.output)
    except (ValueError, OSError, KeyError, IndexError, BadZipFile) as error:
        parser.exit(1, f'{error}\n')
    print(f'PASS: all {count} successful pBCs are present in both tables in {args.output}')


if __name__ == '__main__':
    main()
