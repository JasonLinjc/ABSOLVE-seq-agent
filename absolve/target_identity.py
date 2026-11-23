"""Repository step 6 and representative-family counting for a reporter panel."""
import csv
from collections import Counter
from functools import lru_cache
from pathlib import Path


class TargetIdentity:
    def __init__(self, oligo_file):
        with Path(oligo_file).open() as handle:
            rows = list(csv.DictReader(handle))
        self.designs = {r['OT'].replace('-', '_'): r['OT_sequence'].upper() for r in rows}
        self.offsets = {r['OT'].replace('-', '_'): len(r['Barcode']) + 23 for r in rows}
        if len(rows) != 79 or len(self.designs) != 79:
            raise ValueError('Expected 79 unique reporter targets')
        self.patterns = {}
        for design in set(self.designs.values()):
            masks = {}
            for pos, base in enumerate(design):
                masks[base] = masks.get(base, 0) | (1 << pos)
            self.patterns[design] = masks, 1 << (len(design) - 1)

    def distance(self, design, sequence):
        # Myers bit-vector algorithm for exact, unit-cost global Levenshtein distance.
        masks, final_bit = self.patterns[design]
        positive, negative, score = -1, 0, len(design)
        for base in sequence:
            eq = masks.get(base, 0)
            xv = eq | negative
            xh = (((eq & positive) + positive) ^ positive) | eq
            ph = negative | ~(xh | positive)
            mh = positive & xh
            score += bool(ph & final_bit) - bool(mh & final_bit)
            ph = (ph << 1) | 1
            positive = (mh << 1) | ~(xv | ph)
            negative = ph & xv
        return score

    @lru_cache(maxsize=32768)
    def unique_match(self, target, sequence):
        expected = self.designs[target]
        best = self.distance(expected, sequence)
        return all(self.distance(other, sequence) > best
                   for name, other in self.designs.items() if name != target)

    def retain(self, target, allele):
        ref, read, edited = allele
        if not edited:
            return True
        offset = self.offsets[target]
        if not ref.replace('-', '')[offset:].startswith(self.designs[target]):
            raise ValueError('Aligned reference disagrees with reporter design: ' + target)
        return self.unique_match(target, read.replace('-', '')[offset:])


def family_counts(alleles):
    top = max(alleles, key=alleles.get)
    tie = any(n == alleles[top] and a[2] != top[2] for a, n in alleles.items())
    return sum(n for a, n in alleles.items() if a[2]), sum(alleles.values()), int(top[2]), int(tie)


def filtered_counts(target, families, good, identity=None):
    reads, representatives, ties = [0, 0], [0, 0], 0
    for pbc, alleles in families.items():
        if not good.get(pbc):
            continue
        if identity is not None:
            alleles = {a:n for a, n in alleles.items() if identity.retain(target, a)}
        if alleles:
            edited, total, representative, tie = family_counts(alleles)
            reads[0] += edited
            reads[1] += total
            representatives[0] += representative
            representatives[1] += 1
            ties += tie
    return reads, representatives, ties
