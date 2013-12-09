import logging
import os
import platform
import unittest
from pbcore.io import CmpH5Reader
from pbcore.deprecated import ReferenceEntry
from pbtools.kineticsTools.KineticWorker import KineticWorker
from pbtools.kineticsTools.ipdModel import IpdModel

from test import TestSetup


class TestSpeed(TestSetup):

    def testSpeed(self):

        contig = self.contigs[0].sequence
        snippetFunc = self.ipdModel.snippetFunc(1, 3, 9)
        ipdFunc = self.ipdModel.predictIpdFunc(1)

        snips = [snippetFunc(x, 0) for x in xrange(1000)]

        pFast = self.ipdModel.gbmModel.getPredictions(snips)
        #pSlow = self.ipdModel.gbmModel.getPredictionsSlow(snips)


if __name__ == '__main__':
    unittest.main()
