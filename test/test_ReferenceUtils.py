
import logging
import unittest
import os.path as op

from kineticsTools import ReferenceUtils
from pbcore.io import AlignmentSet

big_data_dir = "/pbi/dept/secondary/siv/testdata/kineticsTools"
ref_dir = "/pbi/dept/secondary/siv/references"

logging.basicConfig()
log = logging.getLogger()


@unittest.skipUnless(op.isdir(big_data_dir), "Shared data folder missing")
class ReferenceUtilsTest (unittest.TestCase):

    def test_bam(self):
        bamFile = op.join(big_data_dir, "ecoli_first_50k.mapped.bam")
        refFile = op.join(ref_dir, "ecoli_k12_MG1655_first50k",
                          "ecoli_k12_MG1655_first50k.referenceset.xml")
        ds = AlignmentSet(bamFile, referenceFastaFname=refFile)
        contigs = ReferenceUtils.loadReferenceContigs(refFile, ds)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].alignmentID, 0)
        chemistry = ReferenceUtils.loadAlignmentChemistry(ds)
        self.assertEquals(chemistry, "S/P3-C3/5.0")

    def test_dataset(self):
        pass  # TODO

    def test_parseReferenceWindow(self):
        window = "gi|12057207|gb|AE001439.1|:1-5000"
        bamFile = op.join(big_data_dir, "Hpyl_1_5000.bam")
        refFile = op.join(ref_dir, "Helicobacter_pylori_J99", "sequence",
                          "Helicobacter_pylori_J99.fasta")
        alnFile = AlignmentSet(bamFile, referenceFastaFname=refFile)
        win = ReferenceUtils.parseReferenceWindow(window,
                                                  alnFile.referenceInfo)
        self.assertEquals([win.refId, win.start, win.end], [0, 1, 5000])

    def test_createReferenceWindows(self):
        bamFile = op.join(big_data_dir, "Hpyl_1_5000.bam")
        ds = AlignmentSet(bamFile, referenceFastaFname=None)
        refInfoTable = ds.referenceInfoTable
        windows = ReferenceUtils.createReferenceWindows(refInfoTable)
        self.assertEqual(len(windows), 1)
        w = windows[0]
        self.assertEqual(w.refId, 0)
        self.assertEqual(w.refName, 'gi|12057207|gb|AE001439.1|')
        self.assertEqual(w.start, 0)
        self.assertEqual(w.end, 1643831)

    def test_enumerateChunks(self):
        pass  # TODO


if __name__ == "__main__":
    unittest.main()
