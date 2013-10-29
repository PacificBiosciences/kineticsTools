#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

from math import sqrt
import math
import scipy.stats as s
import array as a

from scipy.optimize import fminbound
from scipy.special import gammaln as gamln
from numpy import log, pi, log10, e, log1p, exp
import numpy as np


log10e = log10(e)
canonicalBaseMap = { 'A': 'A', 'C':'C', 'G':'G', 'T':'T', 'H':'A', 'I':'C', 'J':'C', 'K':'C' }
modNames = { 'H':'m6A', 'I':'m5C', 'J':'m4C', 'K':'m5C' }
ModificationPeakMask = { 'm6A' : [0, -5], 'm4C': [0, -5], 'm5C': [2, 0, -1, -2, -4, -5, -6]  }

# Labels for modified fraction:

FRAC = 'frac'
FRAClow = 'fracLow'
FRACup = 'fracUp'

# Try computing these only once

k1 = s.norm.ppf(0.025)
k2 = s.norm.ppf(0.975)


class MultiSiteCommon(object):

    def __init__(self, gbmModel, sequence, rawKinetics ):

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
        self.modPriors = { 'H': log1p, 'I': log1p, 'J': log1p, 'K': log1p }

        self.gbmModel = gbmModel
        self.sequence = sequence


        # These switch because we changing viewpoints
        self.pre = gbmModel.post
        self.post = gbmModel.pre

        self.lStart = self.pre
        self.lEnd = len(self.sequence) - self.post

        # Extents that we will use for likelihoods
        self.likelihoodRange = xrange(self.lStart, self.lEnd)
        self.alternateBases = dict((x, set(sequence[x])) for x in xrange(len(sequence)))

        self.rawKinetics = rawKinetics



    def _possibleConfigs(self, start, end):

        if start == end:
            return self.alternateBases[start]
        else:
            r = []
            currentChars = self.alternateBases[start]
            for suffix in self._possibleConfigs(start+1, end):
                for c in currentChars:
                    r.append(c + suffix)

            return r


    def getConfigs(self, centerIdx):
        start = centerIdx - self.pre
        end = centerIdx + self.post
        return self._possibleConfigs(start, end)



    def computeContextMeans(self):
        """Generate a hash of the mean ipd for all candidate contexts"""
        allContexts = list(set([ cfg for pos in self.likelihoodRange for cfg in self.getConfigs(pos) ]))
        predictions = self.gbmModel.getPredictions(allContexts)
        self.contextMeanTable = dict(zip(allContexts, predictions))


    # Log-t pdf - copied from scipy distributions.py line 3836

    def _logpdf(self, x, df):
        r = df*1.0
        lPx = gamln((r+1)/2)-gamln(r/2)
        lPx -= 0.5*log(r*pi) + (r+1)/2*log(1+(x**2)/r)
        return lPx


    def singleScore(self, position, context):
        if self.rawKinetics.has_key(position):
            siteObs = self.rawKinetics[ position ]

            # mu of model, error in model
            um = self.contextMeanTable[context]

            # FIXME -- unify this with the error model used in KineticWorker.py
            # em = 0.06 * um + 0.12 * um**2.0
            em = 0.01 + 0.03*um + 0.06*um**(1.7)

            uo = siteObs['tMean']
            eo = siteObs['tErr']

            t = -(uo - um) / sqrt(em**2 + eo**2)
            df = max(1, siteObs['coverage'] - 1)

            logLikelihood = self._logpdf(t, df).item()
            # logLikelihood = s.t.logpdf(t, df).item()
        else:
            logLikelihood = 0
        
        return logLikelihood


    def scorePosition(self, position, context):
        """ Compute the likelihood of the observed IPDs at position, given the context"""

        # Handle the prior for a modification at the current base here
        # unmodified bases get a prior of 0, modified bases get a prior less than 0.
        prior = 0.0
        if self.modPriors.has_key(context[self.pre]):
            prior = self.modPriors[context[self.pre]]

        # Handle positions where we don't have enough coverage
        if not self.rawKinetics.has_key( position ):
            return prior

        ll = self.singleScore( position, context )
        # return logLikelihood.item() + prior
        return ll + prior



    # Return expected IPDs for a portion [start, end] of the sequence.

    def getContextMeans(self, start, end, sequence):
        meanVector = []
        for pos in xrange(start, end+1):
            ctx = sequence[(pos-self.pre):(pos + self.post + 1)].tostring()
            if self.contextMeanTable.has_key(ctx):
                meanVector.append(self.contextMeanTable[ctx])
            else:
                meanVector.append(self.gbmModel.getPredictions([ctx]))
        return meanVector



