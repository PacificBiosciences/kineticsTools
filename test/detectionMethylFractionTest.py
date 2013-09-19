import logging
import os
import platform
import unittest
from pbcore.io import CmpH5Reader
from pbcore.deprecated import ReferenceEntry
from pbtools.kineticsTools.KineticWorker import KineticWorker
from pbtools.kineticsTools.ipdModel import IpdModel

from test import TestSetup


class TestDetectionMethylFraction(TestSetup):

    # We inherit the setup method for test.py.
    # If you need to customize your dataset, we should set up some different conventions
    # def setUp(self):

    def getOpts(self):
        opts = self.basicOpts()
        opts.identify = False
        opts.methylFraction = True
        return opts

    def testSmallDecode(self):
        """
        Test modified fraction estimation in detection mode around a known modification in lambda
        """

        # First methlyated A in lambda:
        # strand            motif onTarget seqid   tpl
        #      0    GCACNNNNNNGTT       On     1 14983

        start = 14900
        end = 15100
        referenceWindow = (1, start, end)
        bounds = (start, end)

        self.kw._prepForReferenceWindow(referenceWindow)
        kinetics = self.kw._summarizeReferenceRegion(bounds, True, False)
        

        # Verify that we detect m6A mods at 14982 and 14991
        m6AMods = [ {'frac': x['frac'], 'fracLow': x['fracLow'], 'fracUp': x['fracUp'], 'tpl': x['tpl'], 'strand': x['strand']} \
                     for x in kinetics if x.has_key('frac') and x['tpl'] in (14982, 14991)]
        print m6AMods

        for mod in m6AMods:
            self.assertGreater(mod["frac"], 0.5)



if __name__ == '__main__':
    unittest.main()
