"""Classify original-library pilot pBC–tBC pairs using frozen Maxi evidence."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import tempfile

from settings import DATA, location

MERGED = location('ABSOLVE_MAXI_MERGED', DATA/'maxi'/'Maxi_merged')
INDIVIDUAL = location('ABSOLVE_MAXI_POOLS', DATA/'maxi'/'pool_qc')
POOLS = [INDIVIDUAL / f'LV_Maxi_{i}' for i in (1, 2, 3)]
MANUSCRIPT = location('ABSOLVE_ANNOTATION', DATA/'upstream'/'test/data/plasmid_barcode_category.tsv.gz')
PUBLISHED_SHA256 = '2a0eb6d1a58652f234e49e47cf9f72807875d3287c40b5c8fc70a36c3b4effad'
DEFAULT_METHOD = 'manuscript'
REQUIRED = {'merged': 1, 'any-pool': 1, 'two-pools': 2, 'all-pools': 3}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verified_reference(folder):
    path = folder / 'plasmid_barcode_category.tsv.gz'
    meta = json.loads((folder / 'summary.json').read_text())
    policy = meta['annotation_policy']
    digest = sha256(path)
    if (meta['status'] != 'complete' or not meta['counts_conserved'] or
            policy['min_exact_fraction'] != 1 or digest != meta['annotation_sha256']):
        raise ValueError(f'Unverified Maxi reference: {folder}')
    return dict(path=str(path), sha256=digest, entropy_threshold=meta['entropy_threshold'],
                max_eligible_tbc=policy['max_eligible_tbc'])


def reconstructed_call(row, source, minimum):
    exact = int(row['target_exact_reads'])
    evaluable = int(row['target_evaluable_reads'])
    reasons = []
    if evaluable == 0:
        reasons.append('no_evaluable_target')
    if exact != evaluable:
        reasons.append('target_mismatch')
    if exact < minimum:
        reasons.append('insufficient_exact_reads')
    if float(row['entropy']) < source['entropy_threshold']:
        reasons.append('low_entropy')
    if int(row['#barcode']) > source['max_eligible_tbc']:
        reasons.append('too_many_tbc')
    return dict(label='Dud' if reasons else 'Good', reasons=reasons, exact=exact, evaluable=evaluable)


def load_evidence(folder, wanted, target, minimum):
    source = verified_reference(folder)
    path = Path(source['path'])
    found = {}
    with gzip.open(path, 'rt') as handle:
        columns = handle.readline().rstrip().split('\t')
        offset = columns.index('plasmid_barcode')
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            pbc = fields[offset]
            if pbc not in wanted:
                continue
            row = dict(zip(columns, fields))
            if row['tbc'] != wanted[pbc]:
                continue
            if pbc in found or row['target_matched'].replace('-', '_') != target:
                raise ValueError(f'Duplicate or conflicting reference record in {folder}')
            found[pbc] = reconstructed_call(row, source, minimum)
    return found, source


def published_evidence(path, wanted, target):
    digest = sha256(path)
    if digest != PUBLISHED_SHA256:
        raise ValueError(f'Published Dedud annotation checksum mismatch: {path}')
    found = {}
    with gzip.open(path, 'rt') as handle:
        columns = handle.readline().rstrip().split('\t')
        ipbc, itarget, igood, ientropy, inum = (columns.index(k) for k in
            ('plasmid_barcode', 'target_matched', 'if_good', 'entropy', '#barcode'))
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            pbc = fields[ipbc]
            # Repository step 5 filters by target, then retains the first row per pBC.
            if pbc not in wanted or fields[itarget] != target or pbc in found:
                continue
            flag, entropy, n = fields[igood], float(fields[ientropy]), int(fields[inum])
            if flag not in ('True','False'):
                raise ValueError(f'Unexpected if_good value: {flag}')
            reasons = []
            if flag != 'True':
                reasons.append('published_if_good_false')
            if not entropy > 1.477:
                reasons.append('low_entropy')
            if n not in (1, 2):
                reasons.append('tbc_count_outside_1_2')
            found[pbc] = dict(label='Dud' if reasons else 'Good', reasons=reasons,
                              if_good=flag == 'True', entropy=entropy, tbc_count=n)
    return found, dict(path=str(path), sha256=digest,
                       commit='271355dfe40cf46b952e82b88dd4521f23f9feb4',
                       rule='if_good and entropy > 1.477 and tbc_count in (1, 2)')


def decide(labels, required):
    good = labels.count('Good')
    missing = labels.count('Unseen')
    if good >= required:
        return 'Good'
    if missing == len(labels):
        return 'Unseen'
    if good + missing < required:
        return 'Dud'
    return 'Undetermined'


def classify(run, methods, minimum, output):
    if not methods or not set(methods) <= {DEFAULT_METHOD, *REQUIRED} or minimum < 1:
        raise ValueError('Select a listed method and a positive exact-read threshold.')
    if output.exists():
        raise FileExistsError(output)
    methods = list(dict.fromkeys(methods))
    summary = json.loads((run/'summary.json').read_text())
    if summary['sample']['dataset'] not in ('HiFi','NoEP'):
        raise ValueError('Expected an original-library pilot.')
    target = summary['target']
    with (run/'prepared/targets.tsv').open() as handle:
        record, = (r for r in csv.DictReader(handle, delimiter='\t') if r['target_id'] == target)
    aligned = Counter()
    with (run/'pbc_alleles.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['alignment'] == target:
                aligned[row['pbc']] += int(row['reads'])
    if sum(aligned.values()) != summary['runs']['lean']['aligned']:
        raise ValueError('Aligned counts disagree with the pilot receipt.')
    wanted = dict.fromkeys(aligned, record['tbc'])
    references = {DEFAULT_METHOD: [MANUSCRIPT], 'merged': [MERGED]}
    refs = {m: references.get(m, POOLS) for m in methods}
    evidence, sources = {}, {}
    for folder in dict.fromkeys(f for folders in refs.values() for f in folders):
        if folder == MANUSCRIPT:
            evidence[folder], sources[folder.name] = published_evidence(folder, wanted, target)
            continue
        evidence[folder], sources[folder.name] = load_evidence(folder, wanted, target, minimum)
    receipt = dict(sample=summary['sample']['sample_name'], target=target, min_exact_reads=minimum if any(m in REQUIRED for m in methods) else None,
                   sources=sources, pilot_receipt_sha256=sha256(run/'summary.json'),
                   classifier_sha256=sha256(Path(__file__)), methods={})
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        stage = Path(temporary)/'result'
        stage.mkdir()
        with (stage/'pbc_labels.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['method','pbc','tbc','label','aligned_reads','pool_evidence_json'])
            for method, folders in refs.items():
                pairs, reads = Counter(), Counter()
                for pbc, count in sorted(aligned.items()):
                    votes = {f.name:evidence[f].get(pbc, {'label':'Unseen'}) for f in folders}
                    if method == DEFAULT_METHOD:
                        label = votes[MANUSCRIPT.name]['label']
                    else:
                        label = decide([v['label'] for v in votes.values()], REQUIRED[method])
                    pairs[label] += 1
                    reads[label] += count
                    writer.writerow([method,pbc,wanted[pbc],label,count,json.dumps(votes,sort_keys=True)])
                receipt['methods'][method] = dict(required_good_pools=REQUIRED.get(method),
                                                  pairs=dict(pairs),aligned_reads=dict(reads))
        (stage/'summary.json').write_text(json.dumps(receipt,indent=2)+'\n')
        stage.rename(output)
    return receipt


def add_options(parser):
    parser.add_argument('--maxi-method', choices=(DEFAULT_METHOD, *REQUIRED), action='append',
                        help='Maxi evidence rule (default: manuscript); repeat to compare methods on identical alignments.')
    parser.add_argument('--min-exact-reads', type=int, default=1,
                        help='Minimum exact reads for reconstructed Maxi references only (default 1).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    add_options(parser)
    args = parser.parse_args()
    methods = args.maxi_method or [DEFAULT_METHOD]
    if args.min_exact_reads != 1 and not any(m in REQUIRED for m in methods):
        parser.error('--min-exact-reads applies to reconstructed Maxi references only.')
    classify(args.run, methods, args.min_exact_reads, args.output)
