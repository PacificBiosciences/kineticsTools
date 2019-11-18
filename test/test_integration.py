
from collections import defaultdict
import unittest
import tempfile
import os.path as op

from pbcommand.testkit import PbIntegrationBase

DATA_DIR = "/pbi/dept/secondary/siv/testdata/kineticsTools"
REF_DIR = "/pbi/dept/secondary/siv/references"
skip_if_no_data = unittest.skipUnless(
    op.isdir(DATA_DIR) and op.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))

EXPECTED_GFF = """\
Escherichia_coli_K12\tkinModCall\tmodified_base\t6307\t6307\t34\t-\t.\tcoverage=56;context=CAGATTAGCACGCTGATGCGCATCAGCGACAAACTGGCGGG;IPDRatio=1.92
Escherichia_coli_K12\tkinModCall\tmodified_base\t6455\t6455\t32\t-\t.\tcoverage=53;context=CCTGCAAGGACTGGATATGCTGATTCTTATTTCACCTGCGA;IPDRatio=1.78
Escherichia_coli_K12\tkinModCall\tmodified_base\t6512\t6512\t32\t-\t.\tcoverage=55;context=GTAATCACAACTATCGATCAACTCATTCTCATTTTTTGCTA;IPDRatio=1.81
Escherichia_coli_K12\tkinModCall\tm6A\t6984\t6984\t34\t-\t.\tcoverage=56;context=TACTGGCGGGTAACGGCACAACCTACATGCCGCTGGAAGGT;IPDRatio=2.02;identificationQv=5
Escherichia_coli_K12\tkinModCall\tm4C\t7111\t7111\t24\t+\t.\tcoverage=19;context=GGAGGCCAGGACGCCGCTGCCGCTGCCGCGTTTGGCGTCGA;IPDRatio=2.45;identificationQv=4
Escherichia_coli_K12\tkinModCall\tmodified_base\t20051\t20051\t44\t+\t.\tcoverage=19;context=ACGTCAAAGGGTGACAGCAGGCTCATAAGACGCCCCAGCGT;IPDRatio=2.54
Escherichia_coli_K12\tkinModCall\tm4C\t20296\t20296\t24\t-\t.\tcoverage=17;context=GGATGCCGGGCAACAGCCCGCATTATGGGCGTTGGCCTCAA;IPDRatio=1.86;identificationQv=9
Escherichia_coli_K12\tkinModCall\tmodified_base\t20449\t20449\t31\t+\t.\tcoverage=18;context=TGTCCGGCGGTGCTTTTGCCGTTACGCACCACCCCGTCAGT;IPDRatio=2.49"""


class TestKineticsTools(PbIntegrationBase):

    def test_help(self):
        args = ["ipdSummary", "--help"]
        self._check_call(args)

    @skip_if_no_data
    def test_ecoli_first_50k(self):
        gff_out = tempfile.NamedTemporaryFile(suffix=".gff").name
        csv_out = tempfile.NamedTemporaryFile(suffix=".csv").name
        ref_set = op.join(REF_DIR, "ecoli_k12_MG1655_first50k",
                          "ecoli_k12_MG1655_first50k.referenceset.xml")
        alignments = op.join(DATA_DIR, "ecoli_first_50k.mapped.bam")
        args = [
            "ipdSummary",
            "--gff", gff_out,
            "--csv", csv_out,
            "--pvalue", "0.001",
            "--identify", "m6A,m4C",
            "--reference", ref_set,
            alignments
        ]
        self._check_call(args)
        self.assertTrue(op.isfile(gff_out))
        self.assertTrue(op.isfile(csv_out))
        records = []
        with open(gff_out, "rt") as gff:
            for line in gff.read().splitlines():
                if not line.startswith("##"):
                    records.append(line)
        for rec1, rec2 in zip(records, EXPECTED_GFF.splitlines()):
            self.assertEqual(rec1, rec2)
        csv_records = []
        with open(csv_out, "rt") as csv:
            for line in csv.read().splitlines():
                csv_records.append(line)
        self.assertEqual(len(csv_records), 3222)
        self.assertEqual(
            csv_records[1], "\"Escherichia_coli_K12\",6214,0,C,8,0.784,0.144,0.631,1.242,33")
