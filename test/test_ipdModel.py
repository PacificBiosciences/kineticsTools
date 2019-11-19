import unittest
import os.path as op

from pbcore.io import AlignmentSet

from kineticsTools.ipdModel import IpdModel
from kineticsTools import ReferenceUtils

from test_integration import REF_DIR, DATA_DIR, skip_if_no_data


class TestIpdModel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ref_set = op.join(REF_DIR, "ecoli_k12_MG1655_first50k",
                          "ecoli_k12_MG1655_first50k.referenceset.xml")
        alignments = op.join(DATA_DIR, "ecoli_first_50k.mapped.bam")
        ds_aln = AlignmentSet(alignments)
        contigs = ReferenceUtils.loadReferenceContigs(
            ref_set,
            alignmentSet=ds_aln)
        from kineticsTools.loader import getIpdModelFilename
        from kineticsTools.ipdSummary import _getResourcePathSpec
        model_file = getIpdModelFilename(None, "S/P2-C2", [_getResourcePathSpec()])
        cls._model = IpdModel(contigs, model_file)

    @skip_if_no_data
    def test_refLength(self):
        self.assertEqual(self._model.refLength(0), 50000)

    def test_snippetFunc(self):
        f = self._model.snippetFunc(0, 5, 5)
        self.assertEqual(f(9, 0), "TTTCATTCTGA")
        self.assertEqual(f(9, 1), "TCAGAATGAAA")

    def test_getReferenceWindow(self):
        w1 = self._model.getReferenceWindow(0, 0, 4, 15)
        self.assertEqual(w1, "TTTCATTCTGA")
        w2 = self._model.getReferenceWindow(0, 1, 4, 15)
        self.assertEqual(w2, "GTCAGAATGAA")
