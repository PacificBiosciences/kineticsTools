from __future__ import absolute_import
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

from .MultiSiteCommon import MultiSiteCommon, canonicalBaseMap, modNames, ModificationPeakMask, FRAC, FRAClow, FRACup, log10e
from .MixtureEstimationMethods import MixtureEstimationMethods


class ModificationDecode(MultiSiteCommon):

    def __init__(self, gbmModel, sequence, rawKinetics, callBounds, methylMinCov, modsToCall=['H', 'J', 'K'], methylFractionFlag=False, useLDAFlag=False):

        MultiSiteCommon.__init__(self, gbmModel, sequence, rawKinetics)

        # Extents that we will attemp to call a modification
        self.callStart = callBounds[0]
        self.callEnd = callBounds[1]
        self.callRange = range(self.callStart, self.callEnd)

        self.methylMinCov = methylMinCov
        self.modsToCall = modsToCall
        self.methylFractionFlag = methylFractionFlag
        self.useLDA = useLDAFlag

    def decode(self):
        """Use this method to do the full modification finding protocol"""

        # Find potential modification sites
        self.findAlternates()

        # Compute all the required mean ipds under all possible composite hypotheses
        self.computeContextMeans()

        # Fill out the forward matrix
        self.fwdRecursion()

        # Trace back the fwd matrix and return modification calls
        modCalls = self.traceback()

        # Compute a confidence for each mod and return results
        return self.scoreMods(modCalls)

    def findAlternates(self):
        """ Use rules about where IPD peaks appear to generate list
            the set of possible modified bases that we will test during decoding."""

        scoreThresholdLow = 16
        scoreThresholdHigh = 19
        seq = self.sequence

        for (pos, peak) in self.rawKinetics.items():
            score = peak['score']

            # if self.useLDA:
            # Try using LDA model to identify putative Ca5C, regardless of scores
            #     if peak.has_key('Ca5C'):
            #         if peak['Ca5C'] < 0:
            #             self.alternateBases[pos].add('K')

            # Exclude points with low score
            if score < scoreThresholdLow:
                continue

            # note -- don't use the tpl in the actual array - use the dict key
            # we have reversed the indexing to deal with the reverse strand

            if self.callStart <= pos < self.callEnd:
                c = seq[pos]

                # On-target A peak
                if 'H' in self.modsToCall and c == 'A' and score > scoreThresholdHigh:
                    self.alternateBases[pos].add('H')

                # On-target C peak
                if 'J' in self.modsToCall and c == 'C' and score > scoreThresholdHigh:
                    self.alternateBases[pos].add('J')

                if 'K' in self.modsToCall:
                    if c == 'C':
                        self.alternateBases[pos].add('K')

                  # peak at -1 or  -2 or -6 of a C -- 5caC

                    if seq[pos - 2] == 'C' and pos - 2 >= self.callStart:
                        self.alternateBases[pos - 2].add('K')

                    if seq[pos + 1] == 'C' and pos + 1 < self.callEnd:
                        self.alternateBases[pos + 1].add('K')

                    if seq[pos + 2] == 'C' and pos + 2 < self.callEnd:
                        self.alternateBases[pos + 2].add('K')

                    if seq[pos + 5] == 'C' and pos + 5 < self.callEnd:
                        self.alternateBases[pos + 5].add('K')

                    if seq[pos + 6] == 'C' and pos + 6 < self.callEnd:
                        self.alternateBases[pos + 6].add('K')

    def fwdRecursion(self):
        start = self.lStart
        end = self.lEnd

        # Likelihood of each configuration at each position
        scores = dict()

        # fwd score matrix and fwd lookback matrix
        fwdScore = dict()
        fwdPrevState = dict()

        # Fill out first column of score & fwd matrix
        scores[start] = dict((cfg, self.scorePosition(start, cfg)) for cfg in self.getConfigs(start))

        # First column of fwd matrix is same a score matrix, with 'None' in the index matrix
        fwdScore[start] = scores[start]
        fwdPrevState[start] = dict((x, None) for x in scores[start].keys())

        for centerPos in range(start + 1, end):

            # Score and fwd column for current position
            scoreCol = dict()
            fwdScoreCol = dict()
            fwdPrevStateCol = dict()

            # Loop over current state options
            for cfg in self.getConfigs(centerPos):
                score = self.scorePosition(centerPos, cfg)
                scoreCol[cfg] = score

                bestPrevState = None
                bestScore = -1e20

                # Loop over previous state options
                for (prevCfg, prevScore) in fwdScore[centerPos - 1].items():
                    if self.compareStates(cfg, prevCfg) and prevScore + score > bestScore:
                        bestScore = prevScore + score
                        bestPrevState = prevCfg

                fwdScoreCol[cfg] = bestScore
                fwdPrevStateCol[cfg] = bestPrevState

            scores[centerPos] = scoreCol
            fwdScore[centerPos] = fwdScoreCol
            fwdPrevState[centerPos] = fwdPrevStateCol

        self.fwdScore = fwdScore
        self.fwdPrevState = fwdPrevState
        self.scores = scores

    def traceback(self):
        """
        Traceback the fwd matrix to get the bset scoring configuration of modifications
        """
        start = self.lStart
        end = self.lEnd

        modCalls = dict()

        def cogBase(cfg):
            return cfg[self.pre]

        pos = end - 1
        currentCol = self.fwdScore[end - 1]
        bestConfig = max(currentCol, key=lambda x: currentCol[x])

        while True:
            if cogBase(bestConfig) != self.sequence[pos]:
                # Found a modification - save it!
                modCalls[pos] = cogBase(bestConfig)

            bestConfig = self.fwdPrevState[pos][bestConfig]
            pos -= 1

            if bestConfig is None:
                break

        # if self.useLDA:
        # allow LDA-predicted sites through to GFF file
        #     for pos in range(start, end):
        #         if self.rawKinetics.has_key(pos):
        #             if self.rawKinetics[pos].has_key('Ca5C'):
        #                 if 'K' in self.modsToCall:
        # cutoff = min( 0, self.rawKinetics[pos]['coverage']/20.0 - 3.0 )
        #                     cutoff = 0
        #                     echoSites = [modCalls[pos + i] for i in [-6, -5, -2, -1, 2] if modCalls.has_key(pos + i)]
        #                 else:
        #                     cutoff = -2.25
        #                     echoSites = [modCalls[pos + i] for i in range(-10,11) if modCalls.has_key(pos + i)]
        #                 if self.rawKinetics[pos]['Ca5C'] < cutoff:
        # so long as those sites are not in the vicinity of a m6A/m4C call
        #                     if 'H' not in echoSites and 'J' not in echoSites:
        #                         modCalls[pos] = 'K'
        # else:
        # remove any non-LDA-predicted sites from modCalls dictionary?
        # if modCalls.has_key(pos):
        # del modCalls[pos]
        # correct adjacent calls:
        # if self.useLDA:
        #     for pos in range(start + 2, end - 2, 2 ):
        #         x = [pos + i for i in range(-2,3) if modCalls.has_key( pos + i) and self.rawKinetics.has_key( pos + i) ]
        #         y = [modCalls[j] for j in x]
        #         if y.count('K') > 1:
        #             tmp = [self.rawKinetics[j]['Ca5C'] for j in x if self.rawKinetics[j].has_key('Ca5C') ]
        #             if len(tmp) > 0:
        #                 lowest = min(tmp)
        #                 for j in x:
        #                     if self.rawKinetics[j].has_key('Ca5C'):
        #                         if self.rawKinetics[j]['Ca5C'] > lowest:
        #                             del modCalls[j]
        #
        # if adjacent m5C calls are made by the LDA, select the one that has the lower LDA score (Ca5C)
        # for pos in range(start, end):
        # if modCalls.has_key(pos) and self.rawKinetics.has_key(pos) and modCalls.has_key(pos+1) and self.rawKinetics.has_key(pos+1):
        # if self.rawKinetics[pos].has_key('Ca5C') and self.rawKinetics[pos+1].has_key('Ca5C'):
        # if self.rawKinetics[pos]['Ca5C'] < self.rawKinetics[pos+1]['Ca5C']:
        # del modCalls[pos+1]
        # else:
        # del modCalls[pos]
        return modCalls

    def scoreMods(self, modCalls):
        """
        For each modification in the best scoring configuration, score a config excluding the current mod against the winning config
        use this value as the Qmod for the deleted modification
        """

        qvModCalls = dict()

        modSeq = a.array('b')
        modSeq.fromstring(self.sequence)

        # Apply the found modifications to the raw sequence
        for (pos, call) in modCalls.items():
            modSeq[pos] = ord(call)

        for (pos, call) in modCalls.items():

            # Score the modified template at all positions affected by this mod
            modScore = self.scoreRegion(pos - self.post, pos + self.pre, modSeq)
            modScores = self.getRegionScores(pos - self.post, pos + self.pre, modSeq)

            if self.methylFractionFlag and pos in self.rawKinetics:
                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:
                    modifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq)

            # Switch back to the unmodified base and re-score
            modSeq[pos] = ord(canonicalBaseMap[call])
            noModScore = self.scoreRegion(pos - self.post, pos + self.pre, modSeq)
            noModScores = self.getRegionScores(pos - self.post, pos + self.pre, modSeq)

            if self.methylFractionFlag and pos in self.rawKinetics:
                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:
                    unModifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq)

            # Put back the modified base
            modSeq[pos] = ord(call)

            # Compute score difference
            llr = modScore - noModScore

            # Convert from LLR to phred-scaled probability of modification
            qModScore = 10 * llr * log10e + 10 * log1p(exp(-llr)) * log10e

            # Figure out which secondary peaks were likely generated by this modification
            # What is the posterior that the peak was generated by this mod?
            maskPos = self.findMaskPositions(pos, modScores, noModScores)

            # FIXME:  Without this, currently, the identificationQv score is too low for many Ca5C sites
            # if self.useLDA:
            #     if self.rawKinetics.has_key(pos):
            #         if self.rawKinetics[pos].has_key('Ca5C'):
            #             llr = -self.rawKinetics[pos]['Ca5C']
            #             qModScore = 100 * llr * log10e + 100*log1p(exp(-llr))*log10e
            if self.methylFractionFlag and pos in self.rawKinetics:

                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:

                    # Instantiate mixture estimation methods:
                    mixture = MixtureEstimationMethods(self.gbmModel.post, self.gbmModel.pre, self.rawKinetics, self.methylMinCov)

                    # Use modifiedMeanVectors and unmodifiedMeanVectors to calculate mixing proportion, and 95% CI limits.
                    methylFracEst, methylFracLow, methylFracUpp = mixture.estimateMethylatedFractions(pos, unModifiedMeanVectors, modifiedMeanVectors, ModificationPeakMask[modNames[call]])

                    qvModCalls[pos] = {'modification': modNames[call], 'QMod': qModScore, 'LLR': llr, 'Mask': maskPos,
                                       FRAC: methylFracEst, FRAClow: methylFracLow, FRACup: methylFracUpp}

                else:
                    qvModCalls[pos] = {'modification': modNames[call], 'QMod': qModScore, 'LLR': llr, 'Mask': maskPos}

            else:
                # Store the full results
                qvModCalls[pos] = {'modification': modNames[call], 'QMod': qModScore, 'LLR': llr, 'Mask': maskPos}

        return qvModCalls

    def scoreRegion(self, start, end, sequence):

        sc = 0.0
        for pos in range(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring().decode("ascii")
            if pos in self.scores:
                sc += self.scores[pos][ctx]

        return sc

    def getRegionScores(self, start, end, sequence):
        scores = np.zeros(end - start + 1)

        for pos in range(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring().decode("ascii")
            if pos in self.scores:
                scores[pos - start] = self.scores[pos][ctx]

        return scores

    def findMaskPositions(self, pos, modScores, noModScores):

        maskPos = []
        start = pos - self.post
        end = pos + self.pre

        for i in range(start, end + 1):
            # Add a neighboring peak to the mask if
            # a) it has a single-site qv > 20
            # b) the observed IPDs are somewhat more likely under the modified hypothesis than the unmodified hypothesis
            if i in self.rawKinetics and self.rawKinetics[i]["score"] > 20:
                if modScores[i - start] - noModScores[i - start] > 1.0:
                    maskPos.append(i - pos)

        return maskPos

    def compareStates(self, current, prev):
        return current[0:-1] == prev[1:]
