"""Analyze a matched NoEP/HiFi pair or reproduce Supplementary Tables 2–5 from the full study panel."""
import argparse
import csv
import os
from pathlib import Path
import re

from maxi_qc import DEFAULT_METHOD
from settings import MANIFEST, use_local

CONDITIONS = ('NoEP', 'HiFi')


def discover():
    groups = {}
    with MANIFEST.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            condition = row['dataset']
            if row['status'] != 'matched' or condition not in CONDITIONS:
                continue
            name = re.fullmatch(condition + r'_pool_(D\d+)_(T\d+)', row['sample_name'])
            if not name:
                raise ValueError('Unexpected original library name: ' + row['sample_name'])
            pair = groups.setdefault(name.groups(), {})
            if condition in pair:
                raise ValueError('Ambiguous library: ' + row['sample_name'])
            pair[condition] = row
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Show manifest pairs and any missing conditions.')
    parser.add_argument('--manuscript-panel', action='store_true', help='Reproduce Supplementary Tables 2–5 using historical processing for the original 26 libraries.')
    parser.add_argument('--target-filter', choices=['manuscript'], help='Cross-target identity check after Maxi filtering (default manuscript).')
    parser.add_argument('--donor')
    parser.add_argument('--replicate')
    parser.add_argument('--maxi-method', choices=[DEFAULT_METHOD, 'historical'])
    parser.add_argument('--maxi-reference', type=Path, help='Saved Maxi FASTQ rebuild for the selected rule (manuscript or historical).')
    parser.add_argument('--output', type=Path, help='New output directory (default: a unique directory under runs/).')
    parser.add_argument('--fresh', action='store_true', help='Start from FASTQs (the default for all read analyses).')
    parser.add_argument('--read-mode', choices=['manuscript', 'historical'], help='Use uniform paired-end processing for routine analysis or historical processing for the submitted panel.')
    parser.add_argument('--baseline-panel', type=Path, help='Earlier panel_comparison.tsv for tracking resolved, persistent and new disagreements.')
    parser.add_argument('--refit-counts', type=Path, help='Fit the equal-donor manuscript panel from a complete existing library_counts.tsv.')
    parser.add_argument('--finish-manuscript', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--executor', choices=['auto', 'local'], default='auto',
                        help='Use the local machine when Slurm is absent, or force local execution.')
    parser.add_argument('--threads', type=int, help='Local worker limit (default 1).')
    args = parser.parse_args()
    if args.read_mode is None:
        args.read_mode = 'historical' if args.manuscript_panel and not args.refit_counts else 'manuscript'
    if args.read_mode == 'historical':
        if not args.manuscript_panel or args.list or args.refit_counts or args.finish_manuscript:
            parser.error('Historical reads require --manuscript-panel and raw FASTQs')
        if args.maxi_method not in (None, 'historical'):
            parser.error('Historical reads require the historical pBC rule')
        args.maxi_method = 'historical'
        args.maxi_reference = args.maxi_reference or os.getenv('ABSOLVE_HISTORICAL_MAXI_REFERENCE')
        if not args.maxi_reference:
            parser.error('Supply --maxi-reference or ABSOLVE_HISTORICAL_MAXI_REFERENCE for the historical six-Maxi rebuild')
    args.maxi_method = args.maxi_method or DEFAULT_METHOD
    local = use_local(args.executor)
    if args.threads is not None and (args.threads < 1 or not local):
        parser.error('--threads requires local execution and a positive integer')
    workers = (args.threads or 1) if local else int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    if args.refit_counts:
        if not args.manuscript_panel or not args.output or any((args.list, args.donor, args.replicate, args.maxi_reference,
                args.maxi_method != DEFAULT_METHOD, args.fresh, args.target_filter,
                args.baseline_panel, args.finish_manuscript)):
            parser.error('Use --manuscript-panel --refit-counts FILE --output DIR to refit counts with their existing filtering')
        from panel_report.manuscript_run import refit_counts
        refit_counts(args.refit_counts, args.output)
        return
    if args.finish_manuscript:
        from panel_report.manuscript_run import finish
        finish(args.finish_manuscript.resolve())
        return
    if args.maxi_method == 'historical' and not args.maxi_reference:
        parser.error('--maxi-method historical requires --maxi-reference')
    if args.baseline_panel and not args.manuscript_panel:
        parser.error('--baseline-panel requires --manuscript-panel')
    if args.list:
        if args.manuscript_panel:
            parser.error('Use --list to list matched donor/replicate pairs.')
        print('donor\treplicate\tNoEP\tHiFi')
        for key, group in sorted(discover().items()):
            print(*key, *(group.get(c, {}).get('sample_name', 'missing') for c in CONDITIONS), sep='\t')
        return
    if args.manuscript_panel and (args.donor or args.replicate):
        parser.error('The manuscript panel covers the entire submitted design.')
    if not args.manuscript_panel and not (args.donor and args.replicate):
        parser.error('Choose --donor and --replicate, or --manuscript-panel.')
    from panel_report.manuscript_run import start
    start(args, local=local, workers=workers)


if __name__ == '__main__':
    main()
