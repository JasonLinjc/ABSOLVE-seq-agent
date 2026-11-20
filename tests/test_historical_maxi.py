"""Guard the historical two-tBC exception, all-read pooling and canonical file checks."""
import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from historical_maxi import build_lists, load_labels


class HistoricalMaxiTest(unittest.TestCase):
    def test_historical_gates(self):
        good_single = 'ACGT' * 4
        low_single = 'A' * 16
        low_double = 'C' * 16
        mixed = 'TGCA' * 4
        triple = 'GATC' * 4
        tbc = 'CAACAACGC'
        rows = [(low_single,tbc,1), (low_double,tbc,1), (low_double,'AAAAAAAAA',0),
                (mixed,tbc,1), (mixed,tbc,0), (good_single,tbc,1),
                (triple,tbc,1), (triple,'AAAAAAAAA',0), (triple,'CCCCCCCCC',0)]
        lines = io.StringIO(''.join(f'{b}\t{t}\t{v}\t1\n' for b,t,v in sorted(rows)))
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            build_lists(lines, dest)
            name = '1617_OT_0000_REF'
            expected = ''.join(b+'\n' for b in sorted([good_single,low_double])).encode()
            digest = hashlib.sha256(expected).hexdigest()
            self.assertEqual(set(load_labels(dest, name, digest)), {good_single,low_double})
            (dest/(name+'.tsv')).write_text(low_single+'\n')
            with self.assertRaises(ValueError):
                load_labels(dest, name, digest)


if __name__ == '__main__':
    unittest.main()
