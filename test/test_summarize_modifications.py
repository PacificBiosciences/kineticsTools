import os.path as op

import pytest

from pbcommand.testkit import PbIntegrationBase

TESTDATA = "/pbi/dept/secondary/siv/testdata/kineticsTools"


class TestSummarizeModifications(PbIntegrationBase):

    @pytest.mark.internal_data
    def test_integration(self):
        args = [
            "python", "-m", "kineticsTools.summarizeModifications",
            op.join(TESTDATA, "basemods.gff"),
            op.join(TESTDATA, "coverage.gff"),
            "basemods_summary.gff",
        ]
        self._check_call(args)
