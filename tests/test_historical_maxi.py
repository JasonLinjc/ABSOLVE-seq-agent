"""Guard the historical two-tBC exception, all-read pooling and canonical file checks."""
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from historical_maxi import RULE_MANIFESTS, build_lists, configured_reference, load_labels, reference
from maxi_qc import sha256


class HistoricalMaxiTest(unittest.TestCase):
    def test_historical_gates(self):
        good_single = 'ACGT' * 4
        low_single = 'A' * 16
        low_double = 'C' * 16
        good_double = 'AGCT' * 4
        mixed = 'TGCA' * 4
        triple = 'GATC' * 4
        tbc = 'CAACAACGC'
        rows = [(low_single,tbc,1), (low_double,tbc,1), (low_double,'AAAAAAAAA',0),
                (mixed,tbc,1), (mixed,tbc,0), (good_single,tbc,1),
                (good_double,tbc,1), (good_double,'AAAAAAAAA',0),
                (triple,tbc,1), (triple,'AAAAAAAAA',0), (triple,'CCCCCCCCC',0)]
        lines = io.StringIO(''.join(f'{b}\t{t}\t{v}\t1\n' for b,t,v in sorted(rows)))
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            build_lists(lines, dest)
            name = '1617_OT_0000_REF'
            expected = ''.join(b+'\n' for b in sorted([good_single,good_double,low_double])).encode()
            digest = hashlib.sha256(expected).hexdigest()
            self.assertEqual(set(load_labels(dest, name, digest)), {good_single,good_double,low_double})
            (dest/(name+'.tsv')).write_text(low_single+'\n')
            with self.assertRaises(ValueError):
                load_labels(dest, name, digest)
            lines.seek(0)
            build_lists(lines, dest, method='manuscript')
            self.assertEqual((dest/(name+'.tsv')).read_text().splitlines(), sorted([good_single, good_double]))

    def test_saved_reference_must_match_requested_rule(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for method, other in [('manuscript', 'historical'), ('historical', 'manuscript')]:
                receipt = dict(method=method, reference_sha256=sha256(RULE_MANIFESTS[method]))
                (folder/'provenance.json').write_text(json.dumps(receipt))
                self.assertEqual(reference(folder, method)['method'], method)
                with self.assertRaises(ValueError):
                    reference(folder, other)
                receipt['reference_sha256'] = sha256(RULE_MANIFESTS[other])
                (folder/'provenance.json').write_text(json.dumps(receipt))
                with self.assertRaises(ValueError):
                    reference(folder, method)

    def test_legacy_historical_config_and_published_source(self):
        source = {'path': '/saved'}
        self.assertEqual(configured_reference({'historical_reference': source}), source)
        self.assertEqual(configured_reference({'maxi_method': 'manuscript', 'rebuilt_reference': source}), source)
        self.assertFalse(configured_reference({'maxi_method': 'manuscript', 'published': '/published'}))


if __name__ == '__main__':
    unittest.main()
