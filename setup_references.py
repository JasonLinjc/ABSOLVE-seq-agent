"""Build the published reporter-panel index or rebuild its pBC membership from original Maxi FASTQs."""
import argparse
import csv
import gzip
import hashlib
import json
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory

from maxi_qc import PUBLISHED_SHA256, sha256
from settings import DATA, REFERENCES

BLOBS = {
    'NovaSeq3_sample_info.csv': '303737304a0f2785cf7fa8fc2945a823412241a7',
    'OT_guide_amplicon_seq.csv': 'a94d103edc68ba9a0a40af03bacd7f311d783eaf',
    'target_oligos_sequenes.csv': 'f88d16bedf6695a52fc78cb78696877fc6510458',
}


def build_panel(source, output):
    target_info = output/'target_info'
    target_info.mkdir()
    for name, expected in BLOBS.items():
        content = (source/name).read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()
        if blob != expected:
            raise ValueError('Upstream file differs from the pinned revision: '+name)
        (target_info/name).write_bytes(content)
    with (target_info/'OT_guide_amplicon_seq.csv').open() as f:
        guides = {row['OT_name']:row for row in csv.DictReader(f)}
    with (target_info/'target_oligos_sequenes.csv').open() as f:
        oligos = list(csv.DictReader(f))
    targets, mapping = [], []
    for row in oligos:
        key = row['OT'].replace('-', '_')
        tbc, sequence = row['Barcode'], row['OT_sequence'].upper()
        amplicon, guide = (guides[key][field].upper() for field in ('amplicon_seq','guide_seq'))
        if amplicon[23:32] != tbc or amplicon[32:64] != sequence or amplicon.count(guide)!=1:
            raise ValueError('Inconsistent target definitions: '+key)
        if amplicon.index(guide) + len(guide) - 4 != 51:
            raise ValueError('Unexpected guide window: '+key)
        targets.append(dict(id=key,tbc=tbc,target=sequence))
        mapping.append(dict(target_id=key,tbc=tbc,amplicon_seq=amplicon,guide_seq=guide))
    if len(targets)!=79 or len({t['id'] for t in targets})!=79 or len(guides)!=79:
        raise ValueError('Expected 79 unique targets')
    for name, rows in [('targets',targets),('target_guide_mapping',mapping)]:
        (output/(name+'.json')).write_text(json.dumps(rows,indent=2)+'\n')
    return targets


def build(upstream):
    if REFERENCES.exists():
        raise FileExistsError(f'Reference directory already exists: {REFERENCES}')
    annotation = upstream/'test/data/plasmid_barcode_category.tsv.gz'
    if sha256(annotation)!=PUBLISHED_SHA256:
        raise ValueError('Published plasmid annotation checksum mismatch')
    REFERENCES.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=REFERENCES.parent) as tmp:
        stage = Path(tmp)/'references'
        stage.mkdir()
        targets = build_panel(upstream/'test/data/target_info', stage)
        published = stage/'published'
        published.mkdir()
        with ExitStack() as handles, gzip.open(annotation,'rt') as f:
            files = {t['id']:handles.enter_context((published/(t['id']+'.tsv')).open('x')) for t in targets}
            for row in csv.DictReader(f,delimiter='\t'):
                key = row['target_matched']
                if key in files:
                    good = row['if_good']=='True' and float(row['entropy'])>1.477 and int(row['#barcode']) in (1,2)
                    files[key].write(f"{row['plasmid_barcode']}\t{int(good)}\n")
        if any(p.stat().st_size==0 for p in published.iterdir()):
            raise ValueError('A target has no published pBC records')
        (published/'source.sha256').write_text(PUBLISHED_SHA256+'\n')
        stage.rename(REFERENCES)
    print(REFERENCES)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('upstream',type=Path,nargs='?',default=DATA/'upstream')
    rebuild = parser.add_mutually_exclusive_group()
    rebuild.add_argument('--published-maxi', type=Path, metavar='FASTQ_DIR', help='Rebuild qualifying pBCs from the six original Maxi FASTQs using the published entropy rule.')
    rebuild.add_argument('--historical-maxi', type=Path, metavar='FASTQ_DIR', help='Rebuild the historical reference from the six original Maxi FASTQs.')
    parser.add_argument('--trimmomatic-jar', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--threads', type=int, default=3)
    args = parser.parse_args()
    if args.historical_maxi or args.published_maxi:
        from historical_maxi import build as build_maxi
        if not args.trimmomatic_jar:
            parser.error('A Maxi FASTQ rebuild requires --trimmomatic-jar')
        method = 'manuscript' if args.published_maxi else 'historical'
        default_dir = DATA / ('published-maxi' if args.published_maxi else 'historical-maxi')
        build_maxi(args.published_maxi or args.historical_maxi, args.trimmomatic_jar,
                   args.output or default_dir, args.threads, method)
    elif args.trimmomatic_jar or args.output:
        parser.error('--trimmomatic-jar and --output require --published-maxi or --historical-maxi')
    else:
        build(args.upstream)
