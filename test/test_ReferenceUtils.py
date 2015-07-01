
import logging
import unittest
import os.path

from kineticsTools.ReferenceUtils import ReferenceUtils

big_data_dir = "/mnt/secondary-siv/testdata/kineticsTools"
ref_dir = "/mnt/secondary-siv/references"

logging.basicConfig()
log = logging.getLogger()

class ReferenceUtilsTest (unittest.TestCase):
    def setUp (self):
        pass

    def test_cmph5 (self):
        dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
        resourcesDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../kineticsTools/resources')
        refFile = os.path.join(dataDir, 'lambda', 'sequence', 'lambda.fasta')
        cmpFile = os.path.join(dataDir, "p4-c2-lambda-mod-decode.cmp.h5")
        contigs = ReferenceUtils.loadReferenceContigs(refFile, cmpFile)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 1)
        chemistry = ReferenceUtils.loadAlignmentChemistry(cmpFile)
        self.assertEquals(chemistry, "P4-C2")

    def test_bam (self):
        bamFile = os.path.join(big_data_dir, "Hpyl_1_5000.bam")
        refFile = os.path.join(ref_dir, "Helicobacter_pylori_J99", "sequence",
            "Helicobacter_pylori_J99.fasta")
        contigs = ReferenceUtils.loadReferenceContigs(refFile, bamFile)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 0)
        chemistry = ReferenceUtils.loadAlignmentChemistry(bamFile)
        self.assertEquals(chemistry, "P6-C4")

    def test_dataset (self):
        pass # TODO

    def test_parseReferenceWindow (self):
        pass # TODO

    def test_enumerateChunks (self):
        pass # TODO

if __name__ == "__main__":
    unittest.main()
