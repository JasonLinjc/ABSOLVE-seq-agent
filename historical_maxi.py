"""Reconstruct the original three-pool Maxi membership and verify the archived panel."""
import csv
import hashlib
import gzip
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from itertools import groupby
from pathlib import Path
from tempfile import TemporaryDirectory

from maxi_qc import sha256

ASSETS = Path(__file__).resolve().parent / 'runtime' / 'maxi'
MANIFEST = ASSETS / 'reference_sha256.tsv'


def targets():
    with MANIFEST.open() as handle:
        return {row['target']: row for row in csv.DictReader(handle, delimiter='\t')}


def load_labels(directory, target, expected):
    data = (Path(directory) / (target + '.tsv')).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError('Historical Maxi membership checksum mismatch for ' + target)
    return dict.fromkeys(data.decode('ascii').splitlines(), True)


def historical_labels(config, target):
    source = config['historical_reference']
    return load_labels(source['path'], target, targets()[target]['archive_sha256'])


def reference(directory):
    directory = Path(directory).resolve(strict=True)
    receipt = json.loads((directory / 'provenance.json').read_text())
    manifest_hash = sha256(MANIFEST)
    if receipt['method'] != 'historical' or receipt['reference_sha256'] != manifest_hash:
        raise ValueError('Invalid historical Maxi reference receipt')
    return dict(path=str(directory), reference_sha256=manifest_hash, provenance_sha256=sha256(directory/'provenance.json'))


def count_pool(job):
    pool, mates, work, jar, threads = job
    output = work / str(pool)
    output.mkdir()
    trimmed = output / 'paired.fq'
    adapter = ASSETS / 'adapters.fa'
    with (work / f'trim-{pool}.log').open('w') as log:
        subprocess.run(['java', '-Xmx4g', '-jar', jar, 'PE', '-threads', str(threads), '-phred33',
                        *mates, trimmed, os.devnull, os.devnull, os.devnull,
                        f'ILLUMINACLIP:{adapter}:2:30:10', 'LEADING:3', 'TRAILING:3',
                        'SLIDINGWINDOW:4:20', 'AVGQUAL:20', 'MINLEN:48', 'CROP:105'],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    script = ('awk -f "$1" "$2" "$3" | sort --parallel="$4" -S 2G -T "$5" | uniq -c '
              """| awk 'BEGIN{OFS="\\t"}{print $2,$3,$4,$1}'""")
    counts = output / 'counts.tsv'
    with counts.open('w') as handle:
        subprocess.run(['bash', '-o', 'pipefail', '-c', script, 'count', ASSETS/'aggregate.awk',
                        ASSETS/'design.tsv', trimmed, str(threads), output], stdout=handle,
                       env=dict(os.environ, LC_ALL='C'), check=True)
    trimmed.unlink()
    return counts


def build_lists(lines, dest):
    with (ASSETS/'design.tsv').open() as f:
        design = {tbc: name for tbc, _, name in csv.reader(f, delimiter='\t')}
    if set(design.values()) != set(targets()):
        raise ValueError('Target design differs from historical panel')
    with ExitStack() as stack:
        outputs = {name: stack.enter_context((dest/(name+'.tsv')).open('w')) for name in design.values()}
        for barcode, observations in groupby(csv.reader(lines, delimiter='\t'), key=lambda row: row[0]):
            eligible = []
            ntbc = 0
            for tbc, rows in groupby(observations, key=lambda row: row[1]):
                ntbc += 1
                if all(row[2] == '1' for row in rows) and tbc in design:
                    eligible.append(design[tbc])
            if ntbc > 2 or not eligible:
                continue
            frequencies = [barcode.count(base) / len(barcode) for base in 'ACGT']
            if ntbc == 1 and -sum(p * math.log2(p) for p in frequencies if p) <= 1.477:
                continue
            for name in eligible:
                outputs[name].write(barcode + '\n')


def export_membership(directory):
    """Publish the combined table and measured per-target checks from the partitions."""
    with gzip.open(directory/'target_pbc.tsv.gz', 'wt', encoding='ascii') as output, (directory/'target_comparison.tsv').open('w') as report:
        output.write('target\tpbc\n')
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['target', 'archive', 'reconstruction', 'shared', 'archive_only',
                         'reconstruction_only', 'archive_sha256', 'reconstruction_sha256'])
        for name, row in targets().items():
            labels = load_labels(directory, name, row['archive_sha256'])
            count = len(labels)
            if count != int(row['archive']):
                raise ValueError('Historical membership count differs for ' + name)
            output.writelines(name + '\t' + pbc + '\n' for pbc in labels)
            writer.writerow([name, row['archive'], count, count, 0, 0, row['archive_sha256'],
                             sha256(directory/(name + '.tsv'))])


def build(raw, jar, output, threads):
    if threads < 1 or output.exists():
        raise ValueError('Choose a new output directory and at least one thread')
    jar = jar.resolve(strict=True)
    version = subprocess.check_output(['java', '-jar', jar, '-version'], text=True).strip()
    if version != '0.36':
        raise ValueError('Expected Trimmomatic 0.36, found ' + version)
    inputs = {}
    for entry in (ASSETS/'inputs.sha256').read_text().splitlines():
        digest, name = entry.split()
        path = raw / Path(name).name
        if sha256(path) != digest:
            raise ValueError('Original Maxi FASTQ checksum mismatch for ' + str(path))
        inputs[path.name] = digest
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix='.maxi-') as scratch:
        work = Path(scratch)
        dest = work/'reference'
        dest.mkdir()
        jobs = []
        concurrency = min(3, threads)
        for pool in (1, 2, 3):
            mates = sorted(raw/name for name in inputs if name.startswith(f'LV_Maxi_{pool}_'))
            jobs.append((pool, mates, work, jar, max(1, threads//concurrency)))
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            counts = list(executor.map(count_pool, jobs))
        with subprocess.Popen(['sort','-m','-t','\t','-k1,1','-k2,2','-k3,3',*counts],
                              env=dict(os.environ, LC_ALL='C'), stdout=subprocess.PIPE, text=True) as merge:
            build_lists(merge.stdout, dest)
            if merge.wait():
                raise RuntimeError('Could not merge the three pools')
        export_membership(dest)
        for log in work.glob('trim-*.log'):
            log.rename(dest/log.name)
        receipt = dict(method='historical', reference_sha256=sha256(MANIFEST), inputs=inputs,
                       trimmomatic_version=version, jar_sha256=sha256(jar), adapter_sha256=sha256(ASSETS/'adapters.fa'),
                       design_sha256=sha256(ASSETS/'design.tsv'), targets=len(targets()),
                       accepted_keys=sum(int(row['archive']) for row in targets().values()))
        (dest/'provenance.json').write_text(json.dumps(receipt, indent=2)+'\n')
        dest.rename(output)
    print(output)
