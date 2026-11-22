import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile

import pandas as pd

from panel_report.collect_historical_families import collect, FULL, WINDOW


class HistoricalCollectionTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        self.batch = root / 'CRISPRessoBatch_on_WT_D1_T1_target_input'
        self.batch.mkdir()
        self.output = root / 'verified'
        self.prefix = 'target_LV-S1_D1_WTCas9_T1'
        self.info = self.batch / 'CRISPResso2Batch_info.json'
        self.info.write_text(json.dumps({'results': {'completed_batch_arr': ['AAA', 'CCC']}}))
        self.alleles = '#Reads\t%Reads\n3\t75\n1\t25\n'
        for pbc in ('AAA', 'CCC'):
            family = self.batch / ('CRISPResso_on_' + pbc)
            family.mkdir()
            (family / 'Alleles_frequency_table_around_sgRNA.txt').write_text(self.alleles)
            with ZipFile(family / 'Alleles_frequency_table.zip', 'w') as archive:
                archive.writestr('Alleles_frequency_table.txt', self.alleles)
            (family / 'CRISPResso2_info.json').write_text(
                json.dumps({'results': {'refs': {'Reference': {'sgRNA_cut_points': [10]}}}}))
        self.first = self.batch / 'CRISPResso_on_AAA'

    def collect(self):
        return collect(self.batch, 'LV_S1', 'WTCas9', self.output)

    def test_missing_window_then_recovery(self):
        window = self.first / 'Alleles_frequency_table_around_sgRNA.txt'
        window.unlink()
        with self.assertRaises(ValueError) as caught:
            self.collect()
        for suffix in (WINDOW, FULL):
            self.assertIn(f'{self.prefix}{suffix}: missing pBC: AAA', str(caught.exception))
        self.assertFalse(self.output.exists())
        window.write_text(self.alleles)
        self.assertEqual(self.collect(), 2)
        for suffix in (WINDOW, FULL):
            table = pd.read_csv(self.output / (self.prefix + suffix), sep='\t')
            self.assertEqual(table.groupby('UMI')['#Reads'].sum().to_dict(), {'AAA': 4, 'CCC': 4})
            self.assertEqual(table['%Reads'].sum(), 100)

    def test_full_table_missing_one_family(self):
        (self.first / 'Alleles_frequency_table.zip').unlink()
        with self.assertRaises(ValueError) as caught:
            self.collect()
        self.assertIn(self.prefix + FULL + ': missing pBC: AAA', str(caught.exception))
        self.assertNotIn(self.prefix + WINDOW, str(caught.exception))
        self.assertFalse(self.output.exists())

    def test_empty_window_and_equal_size_identity_swap(self):
        (self.first / 'Alleles_frequency_table_around_sgRNA.txt').write_text('#Reads\t%Reads\n')
        self.info.write_text(json.dumps({'results': {'completed_batch_arr': ['CCC', 'GGG']}}))
        with self.assertRaises(ValueError) as caught:
            self.collect()
        for suffix in (WINDOW, FULL):
            self.assertIn(self.prefix + suffix + ': missing pBC: GGG', str(caught.exception))
        self.assertIn(self.prefix + FULL + ': unexpected pBC: AAA', str(caught.exception))
        self.assertNotIn(self.prefix + WINDOW + ': unexpected', str(caught.exception))
        self.assertFalse(self.output.exists())
