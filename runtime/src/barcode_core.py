#!/usr/bin/env python3
"""Read-only ABSOLVE-seq pBC annotation. Python >=3.9; standard library only."""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sqlite3
import sys
import time

VERSION = '1.0.0'
COUNTERS = ('reads', 'target_evaluable', 'target_exact', 'target_variant',
            'target_unresolved', 'pbc_qsum', 'pbc_qbases', 'r2_tbc_tested',
            'r2_tbc_agree', 'r2_pbc_tested', 'r2_pbc_agree',
            'r2_target_tested', 'r2_target_agree')
TRANS = str.maketrans('ACGTN', 'TGCAN')


def rc(seq):
    return seq.translate(TRANS)[::-1]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def identity(path):
    p = Path(path).resolve(strict=True)
    s = p.stat()
    return {'path': str(p), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def hamming(a, b, limit):
    if len(a) != len(b):
        return limit + 1
    n = 0
    for x, y in zip(a, b):
        n += x != y
        if n > limit:
            break
    return n


def quality_ok(q, min_base, min_mean):
    values = [ord(x) - 33 for x in q]
    return bool(values) and min(values) >= min_base and sum(values) >= min_mean * len(values)


def load_config(policy_path, targets_path):
    policy = json.loads(Path(policy_path).read_text())
    targets = json.loads(Path(targets_path).read_text())
    if policy.get('schema_version') != 1 or policy.get('pbc_length') != 16:
        raise ValueError('Only schema 1 and N16 pBC are supported')
    if policy['constant'] != 'CCAACCTCATAGAACACTCATCC':
        raise ValueError('This implementation requires the verified 23nt constant')
    if not 0 <= policy['constant_max_mismatches'] <= 2:
        raise ValueError('Constant mismatch tolerance must be 0..2')
    if not 0 <= policy['entropy_percentile'] <= 100:
        raise ValueError('Invalid entropy percentile')
    if not 0 < policy['min_exact_fraction'] <= 1 or policy['min_exact_support'] < 1:
        raise ValueError('Invalid target support policy')
    if policy['max_eligible_tbc'] < 1:
        raise ValueError('max_eligible_tbc must be positive')
    if not 8 <= len(policy['adapter']) <= 32 or set(policy['adapter']) - set('ACGT'):
        raise ValueError('Invalid adapter motif')
    for k in ['pbc', 'tbc', 'target']:
        if not 0 <= policy['min_' + k + '_base_q'] <= policy['min_' + k + '_mean_q'] <= 93:
            raise ValueError('Invalid quality thresholds')
    if policy['structures'].get('original') != '':
        raise ValueError('The original structure must have an empty prefix')
    for prefix in policy['structures'].values():
        if set(prefix) - set('ACGT'):
            raise ValueError('Invalid canonical prefix')
    if len(set(policy['structures'].values())) != len(policy['structures']):
        raise ValueError('Duplicate canonical prefixes')
    by_tbc, ids = {}, set()
    for t in targets:
        if len(t['tbc']) != 9 or len(t['target']) != 32:
            raise ValueError('Reference must contain 9nt tBC and 32nt target')
        if set(t['tbc'] + t['target']) - set('ACGT'):
            raise ValueError('Ambiguous reference sequence')
        if t['id'] in ids or t['tbc'] in by_tbc:
            raise ValueError('Duplicate target ID or tBC')
        if policy['adapter'] in t['target']:
            raise ValueError('Adapter motif occurs inside a target')
        ids.add(t['id'])
        by_tbc[t['tbc']] = t
    if not by_tbc:
        raise ValueError('Empty reference')
    return policy, by_tbc


def read_record(handle, source, number):
    h = handle.readline()
    if not h:
        return None
    lines = [handle.readline() for _ in range(3)]
    if any(not x for x in lines):
        raise ValueError(f'{source}: incomplete FASTQ record {number}')
    h = h.rstrip('\r\n')
    seq, plus, q = [x.rstrip('\r\n') for x in lines]
    if not h.startswith('@') or not h[1:].split() or not plus.startswith('+'):
        raise ValueError(f'{source}: malformed FASTQ record {number}')
    if not seq or len(seq) != len(q) or set(seq) - set('ACGTN'):
        raise ValueError(f'{source}: invalid sequence/quality length or bases at record {number}')
    if any(ord(c) < 33 or ord(c) > 126 for c in q):
        raise ValueError(f'{source}: invalid Phred+33 encoding at record {number}')
    return h, seq, q


def read_id(header, mate):
    fields = header.split()
    token = fields[0][1:]
    if token.endswith(('/1', '/2')):
        if token[-1] != str(mate):
            raise ValueError('Wrong mate suffix in FASTQ header')
        token = token[:-2]
    if len(fields) > 1 and fields[1].startswith(('1:', '2:')) and fields[1][0] != str(mate):
        raise ValueError('Wrong Illumina mate field in FASTQ header')
    return token


def pairs(r1, r2, max_pairs=None):
    def opener(path):
        return gzip.open(path, 'rt', encoding='ascii', newline='') if str(path).endswith('.gz') else open(path, encoding='ascii', newline='')
    with opener(r1) as a, opener(r2) as b:
        n = 0
        while max_pairs is None or n < max_pairs:
            x, y = read_record(a, r1, n + 1), read_record(b, r2, n + 1)
            if x is None or y is None:
                if x is not None or y is not None:
                    raise ValueError(f'Mate counts differ at record {n + 1}')
                return
            if read_id(x[0], 1) != read_id(y[0], 2):
                raise ValueError(f'Mate IDs differ at record {n + 1}: {x[0]} / {y[0]}')
            n += 1
            yield x, y


class Recognizer:
    def __init__(self, policy, targets):
        self.p = policy
        self.targets = targets
        self.constant = policy['constant']
        self.structures = list(policy['structures'].items())

    def recognize(self, s):
        found, unknown_tbc = [], False
        for name, prefix in self.structures:
            y = len(prefix) + 16
            if not s.startswith(prefix) or len(s) < y + 32:
                continue
            constant = s[y:y + 23]
            if 'N' in constant or hamming(constant, self.constant, self.p['constant_max_mismatches']) > self.p['constant_max_mismatches']:
                continue
            tbc = s[y + 23:y + 32]
            if tbc not in self.targets:
                unknown_tbc = True
                continue
            found.append((name, prefix, y, s[len(prefix):y], tbc))
        if len(found) > 1:
            return None, 'ambiguous_structure'
        if found:
            return found[0], None
        if unknown_tbc:
            return None, 'unknown_tbc'
        # Diagnostic only: shifted exact constants never silently redefine a pBC.
        y = s.find(self.constant)
        if y >= 0:
            if y == 40 and not s.startswith(self.p['structures'].get('STUFFER_CACC', '\0')):
                return None, 'noncanonical_prefix'
            matching = [(len(prefix), name) for name, prefix in self.structures if s.startswith(prefix)]
            prefix_len, _ = max(matching)
            if y - prefix_len != 16:
                return None, 'non16_or_unsupported_prefix'
            if len(s) < y + 32:
                return None, 'partial_barcode_region'
        return None, 'unrecognized_structure'

    def target_region(self, s, q, start):
        end = start + 32
        adapter = self.p['adapter']
        hit = s.find(adapter, start)
        if hit == end:
            boundary = 'exact_adapter'
        elif len(s) >= end + len(adapter) and hamming(s[end:end + len(adapter)], adapter, 2) <= 2:
            # Approximate only a known adapter at the expected boundary; never correct target bases.
            hit, boundary = end, 'approx_adapter_at_expected_boundary'
        elif hit >= start:
            boundary = 'exact_adapter_length_variant'
        elif len(s) == end:
            hit, boundary = end, 'read_end'
        else:
            return None, None, 'partial_target' if len(s) < end else 'unresolved_adapter_boundary'
        target = s[start:hit]
        qual = q[start:hit]
        if not target or 'N' in target:
            return None, None, 'ambiguous_or_empty_target'
        if not quality_ok(qual, self.p['min_target_base_q'], self.p['min_target_mean_q']):
            return None, None, 'target_low_quality'
        return target, boundary, None

    def process(self, one, two):
        _, s, q = one
        _, t, _ = two
        found, reason = self.recognize(s)
        if found is None:
            return None, reason, None
        name, prefix, y, pbc, tbc = found
        if set(pbc) - set('ACGT'):
            return None, 'ambiguous_pbc', None
        pq, bq = q[len(prefix):y], q[y + 23:y + 32]
        if not quality_ok(pq, self.p['min_pbc_base_q'], self.p['min_pbc_mean_q']):
            return None, 'pbc_low_quality', None
        if not quality_ok(bq, self.p['min_tbc_base_q'], self.p['min_tbc_mean_q']):
            return None, 'tbc_low_quality', None
        target, boundary, issue = self.target_region(s, q, y + 32)
        values = [0] * len(COUNTERS)
        values[0] = 1
        values[5], values[6] = sum(ord(x) - 33 for x in pq), len(pq)
        if target is None:
            values[4] = 1
        else:
            values[1] = 1
            values[2 if target == self.targets[tbc]['target'] else 3] = 1
            n = len(target)
            for offset, start, length, expected in [(7, n, 9, rc(tbc)), (9, n + 32, 16, rc(pbc)), (11, 0, n, rc(target))]:
                if len(t) >= start + length:
                    values[offset] = 1
                    values[offset + 1] = int(t[start:start + length] == expected)
        return (pbc, tbc, name, values), 'barcode_eligible', issue or boundary


def open_counts(path):
    db = sqlite3.connect(path)
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('PRAGMA synchronous=NORMAL')
    db.execute('PRAGMA temp_store=FILE')
    db.execute('PRAGMA cache_size=-32768')
    columns = ','.join(k + ' INTEGER NOT NULL' for k in COUNTERS)
    db.execute('CREATE TABLE counts(pbc TEXT NOT NULL,tbc TEXT NOT NULL,structure TEXT NOT NULL,' + columns + ', PRIMARY KEY(pbc,tbc,structure)) WITHOUT ROWID')
    db.execute('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)')
    return db


def set_metadata(db, key, value):
    db.execute('INSERT OR REPLACE INTO metadata VALUES(?,?)', (key, json.dumps(value, sort_keys=True)))


def flush(db, cache):
    placeholders = ','.join('?' for _ in range(3 + len(COUNTERS)))
    updates = ','.join(k + '=' + k + '+excluded.' + k for k in COUNTERS)
    db.executemany('INSERT INTO counts VALUES(' + placeholders + ') ON CONFLICT(pbc,tbc,structure) DO UPDATE SET ' + updates,
                   ((*key, *counts) for key, counts in cache.items()))
    db.commit()
    cache.clear()


def extract(manifest, pool, policy_path, targets_path, output, max_pairs=None, cache_keys=100000):
    if max_pairs is not None and max_pairs < 1 or cache_keys < 1:
        raise ValueError('max_pairs and cache_keys must be positive')
    policy, targets = load_config(policy_path, targets_path)
    samples = json.loads(Path(manifest).read_text())
    ids = [x['pool_id'] for x in samples]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate pool IDs')
    matches = [x for x in samples if x['pool_id'] == pool]
    if len(matches) != 1:
        raise ValueError('Pool must identify exactly one R1/R2 pair')
    sample = matches[0]
    inputs = [identity(sample[k]) for k in ['r1', 'r2']]
    if inputs[0]['path'] == inputs[1]['path']:
        raise ValueError('R1 and R2 must be different files')
    for key, observed in zip(['r1_identity', 'r2_identity'], inputs):
        if key in sample and sample[key] != observed:
            raise ValueError('Input does not match frozen manifest: ' + key)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(output) + '.partial')
    if output.exists():
        raise FileExistsError(output)
    with open(tmp, 'xb'):
        pass
    db = open_counts(tmp)
    rec = Recognizer(policy, targets)
    dispositions, diagnostics, structures = collections.Counter(), collections.Counter(), collections.Counter()
    examples, cache = {}, {}
    n = 0
    start = time.monotonic()
    meta = {'version': VERSION, 'source_sha256': sha256(__file__), 'python': platform.python_version(),
            'pool_id': pool, 'sample': sample, 'inputs': inputs, 'policy': policy,
            'policy_sha256': sha256(policy_path), 'targets': list(targets.values()),
            'targets_sha256': sha256(targets_path), 'manifest_sha256': sha256(manifest),
            'scope': 'pilot' if max_pairs else 'full', 'max_pairs': max_pairs,
            'slurm_job_id': os.getenv('SLURM_JOB_ID'), 'started_utc': utcnow()}
    set_metadata(db, 'run', meta)
    set_metadata(db, 'status', 'running')
    db.commit()
    try:
        for one, two in pairs(inputs[0]['path'], inputs[1]['path'], max_pairs):
            n += 1
            row, disposition, diagnostic = rec.process(one, two)
            dispositions[disposition] += 1
            if diagnostic:
                diagnostics[diagnostic] += 1
            if row:
                pbc, tbc, name, values = row
                structures[name] += 1
                key = (pbc, tbc, name)
                if key in cache:
                    counts = cache[key]
                    for i, value in enumerate(values):
                        counts[i] += value
                else:
                    cache[key] = values
                if len(cache) >= cache_keys:
                    flush(db, cache)
            elif disposition not in examples:
                examples[disposition] = {'record_number': n, 'id': one[0], 'r1_prefix': one[1][:80]}
        flush(db, cache)
        if not n:
            raise ValueError('No read pairs')
        if inputs != [identity(sample[k]) for k in ['r1', 'r2']]:
            raise ValueError('Inputs changed while processing')
        if sum(dispositions.values()) != n:
            raise AssertionError('Disposition counts do not conserve input reads')
        totals = db.execute('SELECT ' + ','.join('COALESCE(SUM(' + k + '),0)' for k in COUNTERS) + ' FROM counts').fetchone()
        if totals[0] != dispositions['barcode_eligible'] or totals[0] != totals[1] + totals[4] or totals[1] != totals[2] + totals[3]:
            raise AssertionError('Counts do not conserve eligible/evaluable reads')
        meta.update({'finished_utc': utcnow(), 'elapsed_seconds': time.monotonic() - start,
                     'pairs': n, 'dispositions': dict(dispositions), 'diagnostics': dict(diagnostics),
                     'structures': dict(structures), 'excluded_examples': examples,
                     'totals': dict(zip(COUNTERS, totals))})
        set_metadata(db, 'run', meta)
        set_metadata(db, 'status', 'complete')
        db.commit()
        db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        db.execute('PRAGMA journal_mode=DELETE')
        db.close()
        # Atomic and no-clobber on the same filesystem.
        os.link(tmp, output)
        tmp.unlink()
        return meta
    except BaseException as exc:
        try:
            set_metadata(db, 'status', 'failed')
            set_metadata(db, 'error', {'message': str(exc), 'pairs_seen': n})
            db.commit()
            db.close()
        except sqlite3.Error:
            pass
        raise


def entropy(seq):
    counts = collections.Counter(seq)
    return -math.fsum((n / len(seq)) * math.log2(n / len(seq)) for n in counts.values())


def percentile(histogram, percentile_value):
    n = sum(histogram.values())
    if n == 0:
        return None
    rank = (n - 1) * percentile_value / 100
    lower, upper = math.floor(rank), math.ceil(rank)
    ranks, seen = {}, 0
    for value, count in sorted(histogram.items()):
        for r in [lower, upper]:
            if seen <= r < seen + count:
                ranks[r] = value
        seen += count
    return ranks[lower] + (rank - lower) * (ranks[upper] - ranks[lower])


def groups(db):
    cols = ','.join('SUM(' + k + ')' for k in COUNTERS)
    cursor = db.execute('SELECT pbc,tbc,' + cols + ',GROUP_CONCAT(structure) FROM counts GROUP BY pbc,tbc ORDER BY pbc,tbc')
    pbc, records = None, []
    for row in cursor:
        if pbc is not None and row[0] != pbc:
            yield pbc, records
            records = []
        pbc = row[0]
        records.append({'tbc': row[1], **dict(zip(COUNTERS, row[2:-1])), 'structures': sorted(row[-1].split(','))})
    if pbc is not None:
        yield pbc, records


def annotate(database, output_dir, annotation_policy=None):
    import csv
    db = sqlite3.connect('file:' + str(Path(database).resolve()) + '?mode=ro', uri=True)
    db.execute('PRAGMA temp_store=FILE')
    metadata = {k: json.loads(v) for k, v in db.execute('SELECT key,value FROM metadata')}
    if metadata.get('status') != 'complete':
        raise ValueError('Refusing incomplete extraction')
    meta = metadata['run']
    policy = dict(meta['policy'])
    if annotation_policy:
        change = json.loads(Path(annotation_policy).read_text())
        allowed = {'min_exact_support', 'min_exact_fraction', 'entropy_percentile', 'max_eligible_tbc'}
        if set(change) - allowed:
            raise ValueError('Only annotation thresholds can change without re-extracting reads')
        policy.update(change)
    if not 0 <= policy['entropy_percentile'] <= 100 or not 0 < policy['min_exact_fraction'] <= 1 or policy['min_exact_support'] < 1 or policy['max_eligible_tbc'] < 1:
        raise ValueError('Invalid annotation thresholds')
    targets = {t['tbc']: t for t in meta['targets']}
    hist = collections.Counter()
    for pbc, in db.execute('SELECT DISTINCT pbc FROM counts ORDER BY pbc'):
        hist[entropy(pbc)] += 1
    cutoff = percentile(hist, policy['entropy_percentile'])
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=False)
    pair_columns = ['pool_id', 'plasmid_barcode', 'target_matched', 'tbc', 'structures',
                    'pair_reads', 'target_evaluable_reads', 'target_exact_reads', 'target_variant_reads',
                    'target_unresolved_reads', 'exact_fraction', 'evaluable_fraction', 'pbc_total_reads',
                    '#barcode', 'n_supported_tbc', 'entropy', 'entropy_threshold', 'mean_pbc_q',
                    'r2_tbc_tested', 'r2_tbc_agree', 'r2_pbc_tested', 'r2_pbc_agree',
                    'r2_target_tested', 'r2_target_agree', 'if_good', 'reason_codes', 'scope']
    reasons = collections.Counter()
    target_summary = {t['id']: collections.Counter() for t in targets.values()}
    summary = {'version': VERSION, 'source_sha256': sha256(__file__), 'extraction': meta,
               'annotation_policy': policy, 'entropy_threshold': cutoff,
               'entropy_quantile_method': 'linear interpolation of unique quality-eligible pBC entropies; strictly below cutoff excluded; ties retained',
               'unique_pbc': sum(hist.values()), 'pbc_tbc_pairs': 0, 'good_pairs': 0, 'good_pbc': 0,
               'counts_conserved': False, 'status': 'running'}
    totals = collections.Counter()
    try:
        with gzip.open(out/'plasmid_barcode_category.tsv.gz.partial', 'wt', newline='') as f, gzip.open(out/'pbc_summary.tsv.gz.partial', 'wt', newline='') as g:
            writer = csv.DictWriter(f, pair_columns, delimiter='\t')
            writer.writeheader()
            family_columns = ['pool_id', 'plasmid_barcode', 'reads', 'n_tbc', 'n_supported_tbc', 'entropy', 'good_target_pairs', 'scope']
            fw = csv.DictWriter(g, family_columns, delimiter='\t')
            fw.writeheader()
            for pbc, records in groups(db):
                h = entropy(pbc)
                n_tbc = len(records)
                n_supported = sum(r['target_exact'] >= policy['min_exact_support'] for r in records)
                reads = sum(r['reads'] for r in records)
                good_family = 0
                for r in records:
                    why = []
                    fraction = r['target_exact'] / r['target_evaluable'] if r['target_evaluable'] else 0.0
                    if cutoff is not None and h < cutoff:
                        why.append('low_entropy')
                    if n_tbc > policy['max_eligible_tbc']:
                        why.append('too_many_eligible_tbc')
                    if r['target_exact'] < policy['min_exact_support']:
                        why.append('insufficient_exact_support')
                    if fraction < policy['min_exact_fraction']:
                        why.append('target_inconsistent' if r['target_evaluable'] else 'no_evaluable_target')
                    good = not why
                    good_family += good
                    reasons.update(why)
                    summary['pbc_tbc_pairs'] += 1
                    summary['good_pairs'] += good
                    totals.update({k: r[k] for k in COUNTERS})
                    row = {'pool_id': meta['pool_id'], 'plasmid_barcode': pbc, 'target_matched': targets[r['tbc']]['id'],
                           'tbc': r['tbc'], 'structures': ','.join(r['structures']), 'pair_reads': r['reads'],
                           'target_evaluable_reads': r['target_evaluable'], 'target_exact_reads': r['target_exact'],
                           'target_variant_reads': r['target_variant'], 'target_unresolved_reads': r['target_unresolved'],
                           'exact_fraction': fraction, 'evaluable_fraction': r['target_evaluable'] / r['reads'],
                           'pbc_total_reads': reads, '#barcode': n_tbc, 'n_supported_tbc': n_supported,
                           'entropy': h, 'entropy_threshold': cutoff, 'mean_pbc_q': r['pbc_qsum'] / r['pbc_qbases'],
                           'if_good': str(good), 'reason_codes': ';'.join(why) if why else 'passes_configured_policy', 'scope': meta['scope']}
                    row.update({k: r[k] for k in COUNTERS if k.startswith('r2_')})
                    writer.writerow(row)
                    tc = target_summary[row['target_matched']]
                    tc['pbc_tbc_pairs'] += 1
                    tc['good_pairs'] += good
                    for field in ['reads', 'target_evaluable', 'target_exact', 'target_variant', 'target_unresolved']:
                        tc[field] += r[field]
                summary['good_pbc'] += good_family > 0
                fw.writerow({'pool_id': meta['pool_id'], 'plasmid_barcode': pbc, 'reads': reads, 'n_tbc': n_tbc,
                             'n_supported_tbc': n_supported, 'entropy': h, 'good_target_pairs': good_family, 'scope': meta['scope']})
        if any(totals[k] != meta['totals'][k] for k in COUNTERS):
            raise AssertionError('Annotation counts differ from extraction')
        summary.update({'counts_conserved': True, 'reason_counts': dict(reasons), 'totals': dict(totals),
                        'status': 'complete', 'finished_utc': utcnow(), 'slurm_job_id': os.getenv('SLURM_JOB_ID')})
        for name in ['plasmid_barcode_category.tsv.gz', 'pbc_summary.tsv.gz']:
            os.link(out/(name+'.partial'), out/name)
            (out/(name+'.partial')).unlink()
        with open(out/'entropy_histogram.tsv', 'x') as f:
            f.write('entropy\tunique_pbc\n')
            for h, n in sorted(hist.items()):
                f.write(f'{h}\t{n}\n')
        with open(out/'target_summary.tsv', 'x', newline='') as f:
            fields = ['target_id', 'pbc_tbc_pairs', 'good_pairs', 'reads', 'target_evaluable', 'target_exact', 'target_variant', 'target_unresolved']
            writer = csv.DictWriter(f, fields, delimiter='\t')
            writer.writeheader()
            for target_id, counts in target_summary.items():
                writer.writerow({'target_id': target_id, **{k: counts[k] for k in fields[1:]}})
        db.close()
        # Completion marker is published last; consumers must require this file.
        with open(out/'summary.json.partial', 'x') as f:
            json.dump(summary, f, indent=2)
        os.link(out/'summary.json.partial', out/'summary.json')
        (out/'summary.json.partial').unlink()
        return summary
    except BaseException as exc:
        with open(out/'FAILED.json', 'w') as f:
            json.dump({'error': str(exc), 'status': 'failed'}, f)
        db.close()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    sub = parser.add_subparsers(dest='command', required=True)
    e = sub.add_parser('extract')
    for name in ['manifest', 'pool', 'policy', 'targets', 'output']:
        e.add_argument('--'+name, required=True)
    e.add_argument('--max-pairs', type=int)
    e.add_argument('--cache-keys', type=int, default=100000)
    a = sub.add_parser('annotate')
    a.add_argument('--database', required=True)
    a.add_argument('--output-dir', required=True)
    a.add_argument('--annotation-policy')
    args = parser.parse_args()
    if args.command == 'extract':
        meta = extract(args.manifest, args.pool, args.policy, args.targets, args.output, args.max_pairs, args.cache_keys)
        print(json.dumps({'pool_id': meta['pool_id'], 'pairs': meta['pairs'], 'scope': meta['scope'], 'elapsed_seconds': meta['elapsed_seconds'], 'dispositions': meta['dispositions']}))
    else:
        result = annotate(args.database, args.output_dir, args.annotation_policy)
        print(json.dumps({k: result[k] for k in ['status', 'unique_pbc', 'pbc_tbc_pairs', 'good_pairs', 'good_pbc', 'counts_conserved']}))


if __name__ == '__main__':
    main()
