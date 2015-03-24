import logging
import os
import platform
import unittest
from pbcore.io import CmpH5Reader
from kineticsTools.KineticWorker import KineticWorker
from kineticsTools.ipdModel import IpdModel
from kineticsTools.ReferenceUtils import ReferenceUtils


class TestSetup(unittest.TestCase):

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

        return opts()

    def setUp(self):

        # Load the lambda genome from our sample data

        dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
        resourcesDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../kineticsTools/resources')
        ref = os.path.join(dataDir, 'lambda', 'sequence', 'lambda.fasta')
        cmpFile = os.path.join(dataDir, "p4-c2-lambda-mod-decode.cmp.h5")

        self.contigs = ReferenceUtils.loadReferenceContigs(ref, cmpFile)
        self.ipdModel = IpdModel(self.contigs, os.path.join(resourcesDir, "P4-C2.h5"))

        # Create a functional KineticWorker object that can be poked at manually.
        self.kw = KineticWorker(self.ipdModel)
        self.cmpH5 = CmpH5Reader(cmpFile)

        # Put in our cmp.h5 - this is normally supplied by the Worker superclass
        self.kw.caseCmpH5 = self.cmpH5
        self.kw.controlCmpH5 = None

        self.kw.options = self.getOpts()

    def runTest(self):
        pass


class TestBasic(TestSetup):

    def _testIpdModel(self):

        contig = self.contigs[0].sequence

        snippetFunc = self.ipdModel.snippetFunc(1, 10, 10)
        snip = snippetFunc(0, 0)

        print "Got snippet at pos 0: %s" % snip
        print "First 10 bases of lambda: %s" % (contig[0:10])

        lastPos = len(contig) - 1
        snip = snippetFunc(lastPos, 0)

        print "Got snippet at pos %d: %s" % (lastPos, snip)
        print "Last 10 bases of lambda: %s" % contig[-10:]

    def _testSpeed(self):

        contig = self.contigs[0].sequence
        snippetFunc = self.ipdModel.snippetFunc(1, 3, 8)
        ipdFunc = self.ipdModel.predictIpdFunc(1)

        snips = [snippetFunc(x, 0) for x in xrange(10000)]

        pFast = self.ipdModel.gbmModel.getPredictions(snips)
        pSlow = self.ipdModel.gbmModel.getPredictionsSlow(snips)

    def testCompareNullToGbm(self):
        """
        Check the null model against a few hard-coded contexts
        """

        contig = self.contigs[0].sequence
        snippetFunc = self.ipdModel.snippetFunc(1, 4, 10)
        ipdFunc = self.ipdModel.predictIpdFunc(1)
        ipdModelFunc = self.ipdModel.predictIpdFuncModel(1)

        for (pos, tplStrand) in [(3, 0), (10, 0), (20, 0), (30, 0), (31, 0), (32, 0), (33, 0), (34, 0)]:
            snip = snippetFunc(pos, tplStrand)

            print "Pos: %d, TplStrand: %d" % (pos, tplStrand)
            print "Got ctx: %s" % snip
            print "From lambda: %s" % (contig[(pos - 4):(pos + 11)])

            print "Lut prediction: %f" % ipdFunc(pos, tplStrand)

            gbmPred = self.ipdModel.gbmModel.getPredictionsSlow([snip])[0]
            print "Gbm prediction: %f" % gbmPred

            gbmPred = self.ipdModel.gbmModel.getPredictions([snip])[0]
            print "Gbm prediction fast: %f" % gbmPred

            gbmSnippetPred = ipdModelFunc(pos, tplStrand)
            print "Gbm pred via predictIpdFuncModel: %f" % gbmSnippetPred

            if snip[4] == 'A':
                snip2 = snip[0:4] + 'H' + snip[5:]
                snip3 = snip[0:9] + 'H' + snip[10:]
                gbmPred = self.ipdModel.gbmModel.getPredictionsSlow([snip2, snip3])
                print "Methylated prediction: %s ->  %f" % (snip2, gbmPred[0])
                print "Methylated prediction: %s ->  %f" % (snip3, gbmPred[1])

                gbmPred = self.ipdModel.gbmModel.getPredictions([snip2, snip3])
                print "Methylated prediction fast: %s ->  %f" % (snip2, gbmPred[0])
                print "Methylated prediction fast: %s ->  %f" % (snip3, gbmPred[1])

            if snip[4] == 'C':
                snip2 = snip[0:4] + 'J' + snip[5:]
                snip3 = snip[0:9] + 'J' + snip[10:]
                gbmPred = self.ipdModel.gbmModel.getPredictionsSlow([snip2, snip3])
                print "Methylated prediction: %s ->  %f" % (snip2, gbmPred[0])
                print "Methylated prediction: %s ->  %f" % (snip3, gbmPred[1])

                gbmPred = self.ipdModel.gbmModel.getPredictions([snip2, snip3])
                print "Methylated prediction fast: %s ->  %f" % (snip2, gbmPred[0])
                print "Methylated prediction fast: %s ->  %f" % (snip3, gbmPred[1])

            print ""

    def testSmallDecode(self):
        """
        Test a modification decode around a known modification in lambda
        """

        # First methlyated A in lambda:
        # strand            motif onTarget seqid   tpl
        #      0    GCACNNNNNNGTT       On     1 14983

        start = 14900
        end = 15100
        referenceWindow = (1, start, end)
        bounds = (start, end)

        self.kw._prepForReferenceWindow(referenceWindow)
        kinetics = self.kw._summarizeReferenceRegion(bounds, False, True)
        mods = self.kw._decodePositiveControl(kinetics, bounds)
        print mods

        # Verify that we detect m6A mods at 14982 and 14991
        m6AMods = [x for x in mods if x['modification'] == 'm6A' and x['tpl'] in (14982, 14991)]
        self.assertEqual(len(m6AMods), 2)

    def _testReferenceBoundary(self):
        start = 0
        end = 400
        referenceWindow = (1, start, end)
        bounds = (start, end)

        res = self.kw.onChunk(referenceWindow)
        # print res


if __name__ == '__main__':
    unittest.main()
