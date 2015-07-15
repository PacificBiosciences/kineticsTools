
import logging
import unittest
import os.path

from kineticsTools.ReferenceUtils import ReferenceUtils
from pbcore.io import AlignmentSet

big_data_dir = "/mnt/secondary-siv/testdata/kineticsTools"
ref_dir = "/mnt/secondary-siv/references"

logging.basicConfig()
log = logging.getLogger()

@unittest.skipUnless(os.path.isdir(big_data_dir), "Shared data folder missing")
class ReferenceUtilsTest (unittest.TestCase):
    def setUp (self):
        pass

    def test_cmph5 (self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        dataDir = os.path.join(base_dir,'data')
        resourcesDir = os.path.join(base_dir, '../kineticsTools/resources')
        refFile = os.path.join(dataDir, 'lambda', 'sequence', 'lambda.fasta')
        cmpFile = os.path.join(dataDir, "p4-c2-lambda-mod-decode.cmp.h5")
        ds = AlignmentSet(cmpFile)
        contigs = ReferenceUtils.loadReferenceContigs(refFile, ds)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 1)
        chemistry = ReferenceUtils.loadAlignmentChemistry(ds)
        self.assertEquals(chemistry, "P4-C2")

    def test_bam (self):
        bamFile = os.path.join(big_data_dir, "Hpyl_1_5000.bam")
        refFile = os.path.join(ref_dir, "Helicobacter_pylori_J99", "sequence",
            "Helicobacter_pylori_J99.fasta")
        ds = AlignmentSet(bamFile)
        contigs = ReferenceUtils.loadReferenceContigs(refFile, ds)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 0)
        chemistry = ReferenceUtils.loadAlignmentChemistry(ds)
        self.assertEquals(chemistry, "P6-C4")

    def test_dataset (self):
        pass # TODO

    def test_parseReferenceWindow (self):
        window = "gi|12057207|gb|AE001439.1|:1-5000"
        bamFile = os.path.join(big_data_dir, "Hpyl_1_5000.bam")
        refFile = os.path.join(ref_dir, "Helicobacter_pylori_J99", "sequence",
            "Helicobacter_pylori_J99.fasta")
        alnFile = AlignmentSet(bamFile)
        win = ReferenceUtils.parseReferenceWindow(window,
            alnFile.referenceInfo)
        self.assertEquals([win.refId, win.start, win.end], [0, 1, 5000])

    def test_createReferenceWindows (self):
        bamFile = os.path.join(big_data_dir, "Hpyl_1_5000.bam")
        ds = AlignmentSet(bamFile)
        refInfoTable = ds.referenceInfoTable
        windows = ReferenceUtils.createReferenceWindows(refInfoTable)
        self.assertEqual(len(windows), 1)
        w = windows[0]
        self.assertEqual(w.refId, 0)
        self.assertEqual(w.refName, 'gi|12057207|gb|AE001439.1|')
        self.assertEqual(w.start, 0)
        self.assertEqual(w.end, 1643831)

    def test_enumerateChunks (self):
        pass # TODO

if __name__ == "__main__":
    unittest.main()
