from math import sqrt
import math
import scipy.stats as s
import array as a

from scipy.optimize import fminbound
from scipy.special import gammaln as gamln
from numpy import log, pi, log10, e, log1p, exp
import numpy as np


log10e = log10(e)
canonicalBaseMap = {'A': 'A', 'C': 'C', 'G': 'G',
                    'T': 'T', 'H': 'A', 'I': 'C', 'J': 'C', 'K': 'C'}
modNames = {'H': 'm6A', 'I': 'm5C', 'J': 'm4C', 'K': 'm5C'}
ModificationPeakMask = {
    'm6A': [0, -5], 'm4C': [0, -5], 'm5C': [2, 0, -1, -2, -4, -5, -6]}

# Labels for modified fraction:

FRAC = 'frac'
FRAClow = 'fracLow'
FRACup = 'fracUp'

# Try computing these only once

k1 = s.norm.ppf(0.025)
k2 = s.norm.ppf(0.975)


class MultiSiteCommon(object):

    def __init__(self, gbmModel, sequence, rawKinetics):
        """
        All indexes are 0-based into the the sequence.

        find a set of sites that _might_ have a modification - each modification type will include a list of
        'neighbor peaks' that can add the current site to the 'options' list.
        6mA and 4mC will use only the on-target peak
        5caC will use on target, -2 and -6.

        Only hits that make this list will be tested in the mod identification process

        Use the viterbi algorithm to find the optimal modifications to include, by measuring the per-site likelihood
        of the observed IPD, given the underlying sequence and methylation states.
        """

        log1p = math.log(0.05)
        self.modPriors = {'H': log1p, 'I': log1p, 'J': log1p, 'K': log1p}

        self.gbmModel = gbmModel
        self.sequence = sequence

        # These switch because we changing viewpoints
        self.pre = gbmModel.post
        self.post = gbmModel.pre

        self.lStart = self.pre
        self.lEnd = len(self.sequence) - self.post

        # Extents that we will use for likelihoods
        self.likelihoodRange = range(self.lStart, self.lEnd)
        self.alternateBases = dict(
            (x, set(sequence[x])) for x in range(len(sequence)))

        self.rawKinetics = rawKinetics

    def _possibleConfigs(self, start, end):

        if start == end:
            return self.alternateBases[start]
        else:
            r = []
            currentChars = self.alternateBases[start]
            for suffix in self._possibleConfigs(start + 1, end):
                for c in currentChars:
                    r.append(c + suffix)

            return r

    def getConfigs(self, centerIdx):
        start = centerIdx - self.pre
        end = centerIdx + self.post
        return self._possibleConfigs(start, end)

    def computeContextMeans(self):
        """Generate a hash of the mean ipd for all candidate contexts"""
        allContexts = list(
            set([cfg for pos in self.likelihoodRange for cfg in self.getConfigs(pos)]))
        predictions = self.gbmModel.getPredictions(allContexts)
        self.contextMeanTable = dict(zip(allContexts, predictions))

    # Log-t pdf - copied from scipy distributions.py line 3836
    def _logpdf(self, x, df):
        r = df * 1.0
        lPx = gamln((r + 1) / 2) - gamln(r / 2)
        lPx -= 0.5 * log(r * pi) + (r + 1) / 2 * log(1 + (x ** 2) / r)
        return lPx

    def singleScore(self, position, context):
        if position in self.rawKinetics:
            siteObs = self.rawKinetics[position]

            # mu of model, error in model
            um = self.contextMeanTable[context]

            # FIXME -- unify this with the error model used in KineticWorker.py
            # em = 0.06 * um + 0.12 * um**2.0
            em = 0.01 + 0.03 * um + 0.06 * um ** (1.7)

            uo = siteObs['tMean']
            eo = siteObs['tErr']

            t = -(uo - um) / sqrt(em ** 2 + eo ** 2)
            df = max(1, siteObs['coverage'] - 1)

            logLikelihood = self._logpdf(t, df).item()
            # logLikelihood = s.t.logpdf(t, df).item()
        else:
            logLikelihood = 0

        return logLikelihood

    def scorePosition(self, position, context):
        """ Compute the likelihood of the observed IPDs at position, given the context"""

        # Handle the prior for a modification at the current base here
        # unmodified bases get a prior of 0, modified bases get a prior less
        # than 0.
        prior = 0.0
        if context[self.pre] in self.modPriors:
            prior = self.modPriors[context[self.pre]]

        # Handle positions where we don't have enough coverage
        if position not in self.rawKinetics:
            return prior

        ll = self.singleScore(position, context)
        # return logLikelihood.item() + prior
        return ll + prior

    # Return expected IPDs for a portion [start, end] of the sequence.
    def getContextMeans(self, start, end, sequence):
        meanVector = []
        for pos in range(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring().decode("ascii")
            if ctx in self.contextMeanTable:
                meanVector.append(self.contextMeanTable[ctx])
            else:
                meanVector.append(self.gbmModel.getPredictions([ctx]))
        return meanVector
