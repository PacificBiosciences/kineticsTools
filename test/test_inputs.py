
"""
Test for BAM file and AlignmentSet support.
"""

import logging
import os
import platform
import unittest
from pbcore.io import AlignmentSet
from kineticsTools.KineticWorker import KineticWorker
from kineticsTools.ipdModel import IpdModel
from kineticsTools.ReferenceUtils import ReferenceUtils, ReferenceWindow

logging.basicConfig()
log = logging.getLogger()

# FIXME
data_dir = "/mnt/secondary-siv/testdata/kineticsTools"

@unittest.skipUnless(os.path.isdir(data_dir), "Missing test data directory")
class TestBam(unittest.TestCase):

    def getOpts(self):
        """Derived tests can override this to customize behaviour"""
        return self.basicOpts()

    def basicOpts(self):
        """Mock up some options for the kinetic worker"""
        class opts:
            def __init__(self):
                self.mapQvThreshold = -1
                self.cap_percentile = 99.0
                self.minCoverage = 3
                self.subread_norm = True
                self.maxCoverage = 200
                self.identify = True
                self.methylFraction = False
                self.pvalue = 0.01
                self.modsToCall = ['H', 'J', 'K']
                # Bug 23546: need to set values for these two new flags:
                self.identifyMinCov = 5
                self.methylMinCov = 10
                self.useLDA = False
                self.maxAlignments = 1500
                self.randomSeed = None
        return opts()

    def getAlignments (self):
        return os.path.join(data_dir, "Hpyl_1_5000.bam")

    def getReference (self):
        refDir = "/mnt/secondary-siv/references"
        return os.path.join(refDir, "Helicobacter_pylori_J99", "sequence",
            "Helicobacter_pylori_J99.fasta")

    def setUp(self):
        self.cmpH5 = None
        resourcesDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../kineticsTools/resources')
        ref = self.getReference()
        alnFile = self.getAlignments()
        assert os.path.exists(alnFile) and os.path.exists(ref)

        self.contigs = ReferenceUtils.loadReferenceContigs(ref, alnFile)
        self.ipdModel = IpdModel(self.contigs, os.path.join(resourcesDir, "P6-C4.h5"))
        # Create a functional KineticWorker object that can be poked at
        self.kw = KineticWorker(self.ipdModel)
        self.cmpH5 = AlignmentSet(alnFile)
        self.cmpH5.addReference(ref)
        # Put in our cmp.h5 - this is normally supplied by the Worker
        self.kw.caseCmpH5 = self.cmpH5
        self.kw.controlCmpH5 = None

        self.kw.options = self.getOpts()

    def testSmallDecode (self):
        """Test for known modifications near the start of H. pylori genome"""
        # XXX should have mods on 60- (m4C), 89+ (m6A), 91- (m6A)
        start = 50
        end = 100
        referenceWindow = ReferenceWindow(0, "gi|12057207|gb|AE001439.1|",
            start, end)
        bounds = (start, end)

        self.kw._prepForReferenceWindow(referenceWindow)
        kinetics = self.kw._summarizeReferenceRegion(bounds, False, True)
        mods = self.kw._decodePositiveControl(kinetics, bounds)
        log.critical(mods)

        # Verify that we detect m6A mods at 14982 and 14991
        m6AMods = [x for x in mods if x['modification'] == 'm6A' and x['tpl'] in (88, 90) ]
        self.assertEqual(len(m6AMods), 2)
        m4CMods = [x for x in mods if x['modification'] == 'm4C' and x['tpl'] in (59,) ]
        self.assertEqual(len(m4CMods), 1)
        for x in mods:
            if x['strand'] == 0:
                self.assertEqual(x['tpl'], 88)
            else:
                self.assertTrue(x['tpl'] in [59,90])


@unittest.skipUnless(os.path.isdir(data_dir), "Missing test data directory")
class TestDataset (TestBam):
    def getAlignment (self):
        return os.path.join(data_dir, "Hpyl_1_5000.xml")

    def testSmallDecode (self):
        pass # TODO


if __name__ == '__main__':
    unittest.main()
