"""Export all four ABSOLVE-seq strategies for a matched original-library NoEP/HiFi pair."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import csv
import gzip
import json
import os
from pathlib import Path
import re
from uuid import uuid4
import subprocess
import sys
import tempfile

from maxi_qc import DEFAULT_METHOD, MERGED, POOLS, REQUIRED, reconstructed_call, sha256, verified_reference
from historical_maxi import reference as historical_reference
from pilot import CUT, PUBLIC, table

AGENT = Path(__file__).resolve().parent
SCRIPTS = AGENT / 'workers'
from settings import DATA, MANIFEST, DEPLOYMENT, PUBLISHED, RUNS, use_local
ARCHIVE = DEPLOYMENT / 'paired_cohort_261119/cohort_results/jobs.json'
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


def previous_alignment(donor, replicate, samples):
    if not ARCHIVE.is_file():
        return None
    entry = json.loads(ARCHIVE.read_text())['pairs'].get(f'{donor}_{replicate}')
    if not entry or entry['samples'] != samples:
        return None
    work = Path(entry['work'])
    if not (work/'results/qc.json').is_file():
        return None
    qc = json.loads((work/'results/qc.json').read_text())
    for c, sample in samples.items():
        inputs, = json.loads((work/c/'inputs.json').read_text())
        if qc[c]['sample_id'] != inputs['pool_id'] or any(inputs[k] != sample[k] for k in ('r1','r2')):
            raise ValueError('Archived alignment has inconsistent sample identity: ' + str(work))
        for mate in ('r1', 'r2'):
            if sha256(Path(sample[mate])) != sample['sha256_'+mate]:
                raise ValueError('Input checksum mismatch: ' + sample[mate])
    return work


def partition_reference(folder, dest, targets, minimum):
    source = verified_reference(folder)
    dest.mkdir()
    with ExitStack() as stack:
        streams = {t: stack.enter_context((dest/(t+'.tsv')).open('x')) for t in targets}
        with gzip.open(source['path'], 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                tid = row['target_matched'].replace('-', '_')
                if tid not in streams:
                    continue
                if row['tbc'] != targets[tid]:
                    raise ValueError('Reference tBC differs from target map: ' + tid)
                if reconstructed_call(row, source, minimum)['label'] == 'Good':
                    streams[tid].write(row['plasmid_barcode']+'\n')
    return source


def execute(args, samples, workers, job_id):
    from maxi_qc import PUBLISHED_SHA256
    targets = {r['target_id']: r['tbc'] for r in json.loads((CUT/'target_guide_mapping.json').read_text())}
    published_targets = {r['OT'].replace('-', '_'): r['Barcode'] for r in table(PUBLIC/'target_oligos_sequenes.csv')}
    if len(targets) != 79 or targets != published_targets:
        raise ValueError('Target map disagrees with the published 79-target panel.')
    work = args.output.resolve() if args.output else RUNS / f"paired_{job_id or ('local_' + uuid4().hex)}"
    work.mkdir(parents=True, exist_ok=False)
    previous = None if args.fresh else previous_alignment(args.donor, args.replicate, samples)
    config = dict(donor=args.donor, replicate=args.replicate, samples=samples, maxi_method=args.maxi_method, target_filter=args.target_filter,
                  oligo_sha256=sha256(PUBLIC/'target_oligos_sequenes.csv'),
                  min_exact_reads=None if args.maxi_method in (DEFAULT_METHOD, "historical") else args.min_exact_reads,
                  manifest_sha256=sha256(MANIFEST), target_map_sha256=sha256(CUT/'target_guide_mapping.json'),
                  job_id=job_id, executor='slurm' if job_id else 'local', workers=workers)
    if previous is not None:
        config['alignment_source'] = str(previous)
    print('Working directory:', work, flush=True)

    def run_tasks(script, indices):
        def run(index):
            command = [sys.executable, str(SCRIPTS/script), str(index)]
            with (work/f'{Path(script).stem}_{index}.log').open('x') as log:
                subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True,
                               env=dict(os.environ, ABSOLVE_DATA=str(DATA), ABSOLVE_MANIFEST=str(MANIFEST)))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(run, indices))

    with tempfile.TemporaryDirectory(dir=work, prefix='maxi_') as tmp:
        if args.maxi_method == DEFAULT_METHOD:
            published = PUBLISHED
            if (published/'source.sha256').read_text().strip() != PUBLISHED_SHA256:
                raise ValueError('Manuscript reference checksum mismatch.')
            config.update(published=str(published), references=[dict(path=str(published), source_sha256=PUBLISHED_SHA256)])
        elif args.maxi_method == 'historical':
            config['historical_reference'] = historical_reference(args.maxi_reference)
            config['references'] = [config['historical_reference']]
        else:
            folders = [MERGED] if args.maxi_method == 'merged' else POOLS
            config['reference_dirs'] = [str(Path(tmp)/f.name) for f in folders]
            config['references'] = [partition_reference(f, Path(dest), targets, args.min_exact_reads)
                                    for f, dest in zip(folders, config['reference_dirs'])]
        (work/'pair.json').write_text(json.dumps(config, indent=2)+'\n')
        if previous is None:
            run_tasks('prepare.py', range(2))
        run_tasks('target.py', range(2*len(targets)))
        run_tasks('summarize.py', [0])
    print(work/'results/paired_target_rates.tsv', flush=True)
    print(work/'results/report.md', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Show manifest pairs and any missing conditions.')
    parser.add_argument('--manuscript-panel', action='store_true', help='Process the full submitted sample design and compare Supplementary Tables 2–5.')
    parser.add_argument('--target-filter', choices=['manuscript','pbc-only'], help='Cross-target identity check after Maxi filtering (default manuscript).')
    parser.add_argument('--donor')
    parser.add_argument('--replicate')
    parser.add_argument('--maxi-method', choices=[DEFAULT_METHOD, 'historical', *REQUIRED], default=DEFAULT_METHOD)
    parser.add_argument('--maxi-reference', type=Path, help='Verified reference directory built with setup_references.py --historical-maxi.')
    parser.add_argument('--min-exact-reads', type=int, help='Reconstructed Maxi methods only (default 1).')
    parser.add_argument('--output', type=Path, help='New output directory (default: a unique directory under runs/).')
    parser.add_argument('--fresh', action='store_true', help='Prepare and align full FASTQs again.')
    parser.add_argument('--read-mode', choices=['fast','manuscript'], help='Read processing mode. Pairs default to fast and manuscript panels default to manuscript.')
    parser.add_argument('--baseline-panel', type=Path, help='Earlier panel_comparison.tsv for tracking resolved, persistent and new disagreements.')
    parser.add_argument('--refit-counts', type=Path, help='Fit the equal-donor manuscript panel from a complete existing library_counts.tsv.')
    parser.add_argument('--finish-manuscript', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--executor', choices=['auto', 'local'], default='auto',
                        help='Use the local machine when Slurm is absent, or force local execution.')
    parser.add_argument('--threads', type=int, help='Local worker limit (default 1).')
    args = parser.parse_args()
    local = use_local(args.executor)
    if args.threads is not None and (args.threads < 1 or not local):
        parser.error('--threads requires local execution and a positive integer')
    workers = (args.threads or 1) if local else int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    job_id = None if local else os.environ.get('SLURM_JOB_ID')
    if args.refit_counts:
        if not args.manuscript_panel or not args.output or any((args.list, args.donor, args.replicate, args.maxi_reference,
                args.maxi_method != DEFAULT_METHOD, args.min_exact_reads is not None, args.fresh, args.target_filter,
                args.read_mode == 'fast', args.baseline_panel, args.finish_manuscript)):
            parser.error('Use --manuscript-panel --refit-counts FILE --output DIR to refit counts with their existing filtering')
        from panel_report.manuscript_run import refit_counts
        refit_counts(args.refit_counts, args.output)
        return
    if args.finish_manuscript:
        from panel_report.manuscript_run import finish
        finish(args.finish_manuscript.resolve())
        return
    if (args.maxi_method == 'historical') != bool(args.maxi_reference):
        parser.error('--maxi-method historical and --maxi-reference must be supplied together')
    if args.maxi_method == 'historical' and args.min_exact_reads is not None:
        parser.error('The historical rule has no adjustable read threshold')
    mode = args.read_mode or ('manuscript' if args.manuscript_panel else 'fast')
    if args.baseline_panel and not (args.manuscript_panel and mode=='manuscript'):
        parser.error('--baseline-panel requires a manuscript read-mode panel')
    if mode == 'manuscript' and not args.list:
        if args.maxi_method not in (DEFAULT_METHOD, "historical") or args.min_exact_reads is not None or args.target_filter == 'pbc-only':
            parser.error('Manuscript reads use manuscript or historical Maxi membership and cross-target filtering.')
        if args.manuscript_panel and (args.donor or args.replicate):
            parser.error('The manuscript panel covers the entire submitted design.')
        if not args.manuscript_panel and not (args.donor and args.replicate):
            parser.error('Choose --donor and --replicate, or --manuscript-panel.')
        from panel_report.manuscript_run import start
        start(args, local=local, workers=workers)
        return
    if args.manuscript_panel:
        if args.list or args.donor or args.replicate or args.maxi_method not in (DEFAULT_METHOD, "historical") or args.min_exact_reads is not None or args.fresh or args.target_filter == 'pbc-only':
            parser.error('--manuscript-panel uses the full original design, manuscript or historical Maxi membership and manuscript target filter.')
        if not local and not job_id:
            parser.error('Submit the panel with sbatch paired.sbatch --manuscript-panel.')
        from panel_report.reproduce import run_panel
        run_panel(args.output, historical_reference(args.maxi_reference) if args.maxi_reference else None, workers=workers, job_id=job_id)
        return
    args.target_filter = args.target_filter or 'manuscript'
    groups = discover()
    if args.list:
        print('donor\treplicate\tNoEP\tHiFi')
        for key, group in sorted(groups.items()):
            print(*key, *(group.get(c, {}).get('sample_name', 'missing') for c in CONDITIONS), sep='\t')
        return
    samples = groups.get((args.donor, args.replicate), {})
    if set(samples) != set(CONDITIONS):
        parser.error('Select a complete NoEP/HiFi pair using --donor and --replicate. See --list.')
    if args.min_exact_reads is not None and (args.min_exact_reads < 1 or args.maxi_method == DEFAULT_METHOD):
        parser.error('--min-exact-reads requires a reconstructed Maxi method and a positive integer.')
    args.min_exact_reads = args.min_exact_reads or 1
    if not local and not job_id:
        parser.error('Run with sbatch paired.sbatch or use --list.')
    execute(args, samples, workers, job_id)


if __name__ == '__main__':
    main()
