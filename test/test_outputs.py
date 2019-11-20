"""
Test sanity of various output formats for a minimal real-world example.
"""

import pytest
import subprocess
import tempfile
import unittest
import os.path
import csv
import re

os.environ["PACBIO_TEST_ENV"] = "1"  # turns off --verbose

DATA_DIR = "/pbi/dept/secondary/siv/testdata/kineticsTools"
REF_DIR = "/pbi/dept/secondary/siv/references/Helicobacter_pylori_J99"
ALIGNMENTS = os.path.join(DATA_DIR, "Hpyl_1_5000.xml")
REFERENCE = os.path.join(REF_DIR, "sequence", "Helicobacter_pylori_J99.fasta")

skip_if_no_data = unittest.skipUnless(
    os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))


@skip_if_no_data
class TestOutputs(unittest.TestCase):

    @classmethod
    @skip_if_no_data
    def setUpClass(cls):
        # prefix = tempfile.NamedTemporaryFile().name  # not sure of the problem, but misused anyway
        prefix = 'TestOutputs'
        cls.csv_file = "{p}.csv".format(p=prefix)
        cls.gff_file = "{p}.gff".format(p=prefix)
        cls.bw_file = "{p}.bw".format(p=prefix)
        args = [
            "ipdSummary", "--log-level", "WARNING",
            "--csv", cls.csv_file,
            "--gff", cls.gff_file,
            #"--bigwig", cls.bw_file,
            "--numWorkers", "12",
            "--pvalue", "0.001",
            "--identify", "m6A,m4C",
            "--referenceStride", "100",
            "--referenceWindows", "gi|12057207|gb|AE001439.1|:0-200",
            "--reference", REFERENCE,
            "--useChemistry", "P6-C4",  # FIXME hacky workaround
            ALIGNMENTS
        ]
        print(" ".join(args))
        assert subprocess.call(args) == 0
        with open(cls.csv_file) as f:
            cls.csv_records = [l.split(",") for l in f.read().splitlines()][1:]

    @classmethod
    def tearDownClass(cls):
        return
        for fn in [cls.csv_file, cls.gff_file, cls.bw_file]:
            if os.path.exists(fn):
                os.remove(fn)

    def test_csv_output(self):
        self.assertEqual(len(self.csv_records), 400)
        self.assertEqual(self.csv_records[0][3], "A")
        self.assertEqual(self.csv_records[100][3], "T")

    @pytest.mark.skip(reason="missing bw_file, so i disabled this (cd)")
    def test_bigwig(self):
        import pyBigWig
        f = pyBigWig.open(self.bw_file)
        for i_rec, rec in enumerate(self.csv_records):
            seqid = re.sub('\"', "", rec[0])
            tpl = int(rec[1]) - 1
            s = int(f.values(seqid, tpl, tpl+1)[0])
            ipd_minus = (s % 65536) / 100.0
            ipd_plus = (s >> 16) / 100.0
            if rec[2] == "1":
                self.assertAlmostEqual(ipd_minus, float(rec[8]), places=1)
            else:
                self.assertAlmostEqual(ipd_plus, float(rec[8]), places=1)


if __name__ == "__main__":
    unittest.main()
