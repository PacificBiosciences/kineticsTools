
import unittest
import os.path as op

from pbcommand.testkit import PbIntegrationBase

TESTDATA = "/pbi/dept/secondary/siv/testdata/kineticsTools"
skip_unless_testdata = unittest.skipUnless(op.isdir(TESTDATA),
                                           "Testdata missing")

class TestSummarizeModifications(PbIntegrationBase):

    @skip_unless_testdata
    def test_integration(self):
        args = [
            "python", "-m", "kineticsTools.summarizeModifications",
            op.join(TESTDATA, "basemods.gff"),
            op.join(TESTDATA, "coverage.gff"),
            "basemods_summary.gff"
        ]
        self._check_call(args)
