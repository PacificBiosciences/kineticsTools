import os.path as op

import pytest

from pbcore.io import AlignmentSet

from kineticsTools.ipdModel import IpdModel
from kineticsTools import ReferenceUtils

from test_integration import REF_DIR, DATA_DIR


@pytest.mark.internal_data
class TestIpdModel:

    @classmethod
    def setup_class(cls):
        ref_set = op.join(REF_DIR, "ecoli_k12_MG1655_first50k",
                          "ecoli_k12_MG1655_first50k.referenceset.xml")
        alignments = op.join(DATA_DIR, "ecoli_first_50k.mapped.bam")
        ds_aln = AlignmentSet(alignments)
        contigs = ReferenceUtils.loadReferenceContigs(
            ref_set,
            alignmentSet=ds_aln)
        from kineticsTools.loader import getIpdModelFilename
        from kineticsTools.ipdSummary import _getResourcePathSpec
        model_file = getIpdModelFilename(
            None, "S/P2-C2", [_getResourcePathSpec()])
        cls._model = IpdModel(contigs, model_file)

    def test_refLength(self):
        assert self._model.refLength(0) == 50000

    def test_snippetFunc(self):
        f = self._model.snippetFunc(0, 5, 5)
        assert f(9, 0) == "TTTCATTCTGA"
        assert f(9, 1) == "TCAGAATGAAA"

    def test_getReferenceWindow(self):
        w1 = self._model.getReferenceWindow(0, 0, 4, 15)
        assert w1 == "TTTCATTCTGA"
        w2 = self._model.getReferenceWindow(0, 1, 4, 15)
        assert w2 == "GTCAGAATGAA"
