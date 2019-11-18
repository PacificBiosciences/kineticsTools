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
import sys

from numpy import log, pi, log10, e, log1p, exp
import numpy as np
import re

log10e = log10(e)

canonicalBaseMap = {'A': 'A', 'C': 'C', 'G': 'G',
                    'T': 'T', 'H': 'A', 'I': 'C', 'J': 'C', 'K': 'C'}
modNames = {'H': 'm6A', 'I': 'm5C', 'J': 'm4C', 'K': 'm5C'}

m5CCode = 'I'

iupacMap = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'K': 'GT',
    'M': 'AC',
    'R': 'AG',
    'Y': 'CT',
    'S': 'CG',
    'W': 'AT',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}


def findMotifPositions(seq, motifs):
    regexs = []

    # Generate a regex for each motif, honouring degenerate bases
    for m in motifs:
        regex = ''

        for c in m:
            regex = regex + "[" + iupacMap[c] + "]"

        regexs.append(regex)

    allMatches = []

    # Return a list of matching positions in the sequence
    for r in regexs:
        rr = re.compile(r)
        matches = [x.start() for x in rr.finditer(seq)]
        allMatches.extend(matches)

    allMatches.sort()

    return allMatches


class MultiSiteDetection(object):

    def __init__(self, gbmModel, sequence, rawKinetics, callBounds, methylMinCov, motifs=['CG']):
        """

        """

        self.methylMinCov = methylMinCov
        self.motifs = motifs

        self.gbmModel = gbmModel
        self.sequence = sequence

        self.callStart = callBounds[0]
        self.callEnd = callBounds[1]

        # Extents that we will attempt to call a modification
        self.callRange = range(self.callStart, self.callEnd)

        # These switch because we changing viewpoints
        self.pre = gbmModel.post
        self.post = gbmModel.pre

        self.lStart = self.pre
        self.lEnd = len(self.sequence) - self.post

        # Extents that we will use for likelihoods
        self.likelihoodRange = range(self.lStart, self.lEnd)

        self.alternateBases = dict(
            (x, list(sequence[x])) for x in range(len(sequence)))

        self.rawKinetics = rawKinetics

    def getConfigs(self, centerIdx):
        ''' Enumerate all the contexts centered at centerIdx with one
            modification added '''
        start = centerIdx - self.pre
        end = centerIdx + self.post
        return self._possibleConfigs(start, end)

    def _possibleConfigs(self, start, end):
        ''' Enumerate all the contexts coming from the substring self.sequence[start,end] with one
            modification added '''

        if start == end:
            return self.alternateBases[start]
        else:
            r = []
            allSuffixes = self._possibleConfigs(start + 1, end)

            # The first suffix is alway the one with no modifications
            # Only add the alternate to that one -- that way we only
            # get configurations with a single modification, not all combos

            noModsSuffix = allSuffixes[0]
            if len(allSuffixes) > 1:
                restSuffixes = allSuffixes[1:]
            else:
                restSuffixes = []

            # The noMods suffix get the alternates
            for c in self.alternateBases[start]:
                r.append(c + noModsSuffix)

            # the other suffixes already have mods -- they just get the unmodified base
            for suffix in restSuffixes:
                r.append(self.alternateBases[start][0] + suffix)

            return r

        # Compute something for all the windows in [start, end]
    def getContexts(self, start, end, sequence):
        contexts = []

        for pos in range(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring()
            contexts.append(ctx)

        return contexts

    def computeContextMeans(self):
        """Generate a hash of the mean ipd for all candidate contexts"""

        allContexts = []

        for pos in self.motifPositions:
            for offsetPos in range(pos - self.post, pos + self.pre + 1):
                cfgs = self.getConfigs(offsetPos)
                allContexts.extend(cfgs)

        predictions = self.gbmModel.getPredictions(allContexts)
        self.contextMeanTable = dict(zip(allContexts, predictions))

    def decode(self):
        """Use this method to do the full modification finding protocol"""

        # Find sites matching the desired motif
        self.findMotifs()

        # Compute all the required mean ipds under all possible composite hypotheses
        self.computeContextMeans()

        # Compute a confidence for each mod and return results
        return self.scorePositions()

    def findMotifs(self):
        """ Mark all the positions matching the requested motif """

        # Generate list of matching positions
        allMotifPositions = findMotifPositions(self.sequence, self.motifs)
        self.motifPositions = []

        for pos in allMotifPositions:
            # Only use bases that are inside the callBounds
            if self.callStart <= pos < self.callEnd:
                self.alternateBases[pos].append('I')
                self.motifPositions.append(pos)

    def multiSiteDetection(self, positions, nullPred, modPred, centerPosition):
        ''' kinetics, nullPred, and modifiedPred are parallel arrays 
            containing the observations and predictions surrounding a 
            single candidate motif site.  Estimate the p-value of
            modification and the modified fraction here'''

        # Apply the error model to the predictions
        nullErr = 0.01 + 0.03 * nullPred + 0.06 * nullPred ** (1.7)
        modErr = 0.01 + 0.03 * modPred + 0.06 * modPred ** (1.7)

        obsMean = np.zeros(nullPred.shape)
        obsErr = np.zeros(nullPred.shape)

        # Get the observations into the same array format
        for i in range(len(positions)):
            position = positions[i]

            if position in self.rawKinetics:
                siteObs = self.rawKinetics[position]
                obsMean[i] = siteObs['tMean']
                obsErr[i] = siteObs['tErr']
            else:
                # Crank up the variance -- we don't have an observation at this
                # position, so we should ignore it.
                obsMean[i] = 0.0
                obsErr[i] = 999999999

        # Subtract off the background model from the observations and the modified prediction
        dObs = obsMean - nullPred
        # Error of observation and prediction are uncorrelated
        obsSigma = obsErr ** 2 + nullErr ** 2
        invObsSigma = 1.0 / obsSigma

        # Error of null prediction and mod prediction are probably correlated -- need a better estimate of the error of the difference!!
        dPred = modPred - nullPred
        # Just stubbing in a factor of 2 here...
        dPredSigma = (obsErr ** 2 + nullErr ** 2) / 2

        weightsNumerator = invObsSigma * dPred
        weights = weightsNumerator / (dPred * weightsNumerator).sum()

        signalEstimate = (weights * dObs).sum()
        varianceEstimate = (np.abs(weights) * obsSigma).sum()

        maxSignal = (weights * dPred).sum()
        maxSignalVariance = (np.abs(weights) * dPredSigma).sum()

        # Now just run the standard erf on this Gaussian to quantify the probability that there is some signal
        # What we want now:
        #
        # 1. p-value that dObs * dPred (dot product) is greater than 0.
        # 2. Distribution of \alpha, where dObs = \alpha dPred, where \alpha \in [0,1], with appropriate error propagation
        # 2a. Is it possible to summarize 2 with a Beta distribution?

        pvalue = s.norm._cdf(-signalEstimate / varianceEstimate)
        pvalue = max(sys.float_info.min, pvalue)
        score = -10.0 * log10(pvalue)

        centerPosition['MSscore'] = score
        centerPosition['MSpvalue'] = pvalue

        centerPosition['signal'] = signalEstimate
        centerPosition['variance'] = varianceEstimate

        centerPosition['modelSignal'] = maxSignal
        centerPosition['modelVariance'] = maxSignalVariance

        centerPosition['Mask'] = []

        return centerPosition

    def scorePositions(self):
        """
        Score each motif site in the sequence.
        """

        qvModCalls = dict()

        dnaSeq = a.array('c')
        dnaSeq.fromstring(self.sequence)

        for pos in self.motifPositions:
            if pos in self.rawKinetics:

                # Fetch unmodified positions
                nullPred = self.getRegionPredictions(
                    pos - self.post, pos + self.pre, dnaSeq)

                # Fetch modified positions and reset sequence
                originalBase = dnaSeq[pos]
                dnaSeq[pos] = m5CCode
                modifiedPred = self.getRegionPredictions(
                    pos - self.post, pos + self.pre, dnaSeq)
                dnaSeq[pos] = originalBase

                # Position that contribute to this call
                positions = range(pos - self.post, pos + self.pre + 1)

                # Run the multi-site detection and save the results
                centerStats = self.rawKinetics[pos]
                centerStats = self.multiSiteDetection(
                    positions, nullPred, modifiedPred, centerStats)

                qvModCalls[pos] = centerStats

        return qvModCalls

    def getRegionPredictions(self, start, end, sequence):
        predictions = np.zeros(end - start + 1)

        for pos in range(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring()
            predictions[pos - start] = self.contextMeanTable[ctx]

        return predictions
