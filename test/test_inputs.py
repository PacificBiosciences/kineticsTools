
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

class _TestBase(object):
    """
    Common test functionality.  All input type tests should inherit from this,
    and yield identical results.
    """

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
        raise NotImplementedError()

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

        self.ds = AlignmentSet(alnFile, referenceFastaFname=ref)
        self.contigs = ReferenceUtils.loadReferenceContigs(ref, self.ds)
        self.ipdModel = IpdModel(self.contigs, os.path.join(resourcesDir, "P6-C4.h5"))
        # Create a functional KineticWorker object that can be poked at
        self.kw = KineticWorker(self.ipdModel)
        # Put in our cmp.h5 - this is normally supplied by the Worker
        self.kw.caseCmpH5 = self.ds
        self.kw.controlCmpH5 = None

        self.kw.options = self.getOpts()

    def test_private_api (self):
        start = 50
        end = 100
        REF_GROUP_ID = "gi|12057207|gb|AE001439.1|"
        referenceWindow = ReferenceWindow(0, REF_GROUP_ID, start, end)
        bounds = (start, end)
        rir = list(self.kw.caseCmpH5.readsInRange(referenceWindow.refName,
            referenceWindow.start, referenceWindow.end))
        self.assertEqual(len(rir), 301)
        chunks = self.kw._fetchChunks(REF_GROUP_ID, (start, end),
            self.kw.caseCmpH5)
        factor = 1.0 / self.ds.readGroupTable[0].FrameRate
        rawIpds = self.kw._loadRawIpds(rir, start, end, factor)
        logging.critical(len(rawIpds))
        # XXX note that this is very dependent on the exact order of reads
        # found by readsInRange(), which may be altered by changes to the
        # implementation of the dataset API.  It should, however, remain
        # consistent across equivalent input types.
        # XXX 2015-08-28 disabling this for now because it will change if the
        # dataset contains multiple .bam files
        #self.assertEqual("%.4f" % rawIpds[0][2], "0.2665")
        log.info(rawIpds)
        chunks = self.kw._chunkRawIpds(rawIpds)
        #log.critical(chunks)

    def test_small_decode (self):
        """Test for known modifications near the start of H. pylori genome"""
        # XXX should have mods on 60- (m4C), 89+ (m6A), 91- (m6A)
        start = 50
        end = 100
        REF_GROUP_ID = "gi|12057207|gb|AE001439.1|"
        referenceWindow = ReferenceWindow(0, REF_GROUP_ID, start, end)
        bounds = (start, end)

        self.kw._prepForReferenceWindow(referenceWindow)
        kinetics = self.kw._summarizeReferenceRegion(bounds, False, True)
        mods = self.kw._decodePositiveControl(kinetics, bounds)
        log.info(mods)

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
class TestBam(_TestBase, unittest.TestCase):
    def getAlignments (self):
        return os.path.join(data_dir, "Hpyl_1_5000.bam")


@unittest.skipUnless(os.path.isdir(data_dir), "Missing test data directory")
class TestDataset (TestBam, unittest.TestCase):
    def getAlignments (self):
        return os.path.join(data_dir, "Hpyl_1_5000.xml")


@unittest.skipUnless(os.path.isdir(data_dir), "Missing test data directory")
class TestSplitDataset(_TestBase, unittest.TestCase):
    def getAlignments (self):
        return os.path.join(data_dir, "Hpyl_1_5000_split.xml")


@unittest.skipUnless(os.path.isdir(data_dir), "Missing test data directory")
class TestChunkedDataset(_TestBase, unittest.TestCase):

    def getAlignments(self):
        return os.path.join(data_dir, "Hpyl_1_5000_chunk.xml")

    @unittest.skip
    def test_private_api(self):
        pass

    def test_small_decode(self):
        start = 985
        end = 1065
        REF_GROUP_ID = "gi|12057207|gb|AE001439.1|"
        referenceWindow = ReferenceWindow(0, REF_GROUP_ID, start, end)
        bounds = (start, end)

        self.kw._prepForReferenceWindow(referenceWindow)
        kinetics = self.kw._summarizeReferenceRegion(bounds, False, True)
        mods = self.kw._decodePositiveControl(kinetics, bounds)
        self.assertEqual(len(mods), 4)


if __name__ == '__main__':
    unittest.main()
