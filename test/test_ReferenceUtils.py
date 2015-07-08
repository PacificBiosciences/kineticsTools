
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
        table1, table2 = ReferenceUtils.loadAlignmentTables(cmpFile)
        self.assertEquals(list(table1[0]), [1, 1, 'lambda_NEB3011',
            'lambda_NEB3011', 48502, 'a1319ff90e994c8190a4fe6569d0822a',0,341])
        self.assertEquals(list(table2[0]), [1,
            'm130830_001214_ethan_c100541172550000001823090111241380_s1_p0',
            'standard', 'P4-C2', 75.0])
        contigs = ReferenceUtils.loadReferenceContigs(refFile, cmpFile)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 1)
        chemistry = ReferenceUtils.loadAlignmentChemistry(cmpFile)
        self.assertEquals(chemistry, "P4-C2")

    def test_bam (self):
        bamFile = os.path.join(big_data_dir, "Hpyl_1_5000.bam")
        refFile = os.path.join(ref_dir, "Helicobacter_pylori_J99", "sequence",
            "Helicobacter_pylori_J99.fasta")
        table1, table2 = ReferenceUtils.loadAlignmentTables(bamFile)
        self.assertEquals(list(table1[0]), [0, 0, 'gi|12057207|gb|AE001439.1|',             'gi|12057207|gb|AE001439.1|', 1643831,
            'cc055ae276c41fc1c6d7fbedd1c1a27f', 0, 0])
        self.assertEquals(list(table2[0]), [660570681,
            'm140907_161732_42142_c100676542550000001823129611271416_s1_p0',
            'SUBREAD', 'P6-C4', 75.0])
        contigs = ReferenceUtils.loadReferenceContigs(refFile, bamFile)
        self.assertEquals(len(contigs), 1)
        self.assertEquals(contigs[0].cmph5ID, 0)
        chemistry = ReferenceUtils.loadAlignmentChemistry(bamFile)
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
        (refInfoTable, _) = ReferenceUtils.loadAlignmentTables(bamFile)
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
