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


class ModificationDecode(object):

    def __init__(self, gbmModel, sequence, rawKinetics, callBounds, methylMinCov, modsToCall = ['H', 'J', 'K'], methylFractionFlag = False, useLDAFlag = False):
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

        self.methylMinCov = methylMinCov

        # Temporary:
        self.useLDA = useLDAFlag

        self.modsToCall = modsToCall

        self.methylFractionFlag = methylFractionFlag	
   
        log1p = math.log(0.05)
        self.modPriors = { 'H': log1p, 'I': log1p, 'J': log1p, 'K': log1p }

        self.gbmModel = gbmModel
        self.sequence = sequence

        self.callStart = callBounds[0]
        self.callEnd = callBounds[1]

        # Extents that we will attemp to call a modification
        self.callRange = xrange(self.callStart, self.callEnd)

        # These switch because we changing viewpoints
        self.pre = gbmModel.post
        self.post = gbmModel.pre

        self.lStart = self.pre
        self.lEnd = len(self.sequence) - self.post

        # Extents that we will use for likelihoods
        self.likelihoodRange = xrange(self.lStart, self.lEnd)

        self.alternateBases = dict((x, set(sequence[x])) for x in xrange(len(sequence)))


        self.rawKinetics = rawKinetics


    def getConfigs(self, centerIdx):
        start = centerIdx - self.pre
        end = centerIdx + self.post
        return self._possibleConfigs(start, end)

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

            if self.useLDA:
                # Try using LDA model to identify putative Ca5C, regardless of scores
                if peak.has_key('Ca5C'):
                    if peak['Ca5C'] < 0:
                        self.alternateBases[pos].add('K')

            # Exclude points with low score
            if score < scoreThresholdLow:
                continue

            # note -- don't use the tpl in the actual array - use the dict key
            # we have reversed the indexing to deal with the reverse strand

            if self.callStart <= pos < self.callEnd:
                c = seq[ pos ]

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

        # FIXME:  temporary change to try LDA for Ca5C detection?       
        # if self.useLDA:
        #     # sum scores over a window around current position:
        #     ll = 0
        #     for offset in xrange(-3, 3):
        #         ll += self.singleScore(position+offset, context) 
        # else:
        #     ll = self.singleScore(position, context)
        # Doesn't seem to work?

        ll = self.singleScore( position, context )
        # return logLikelihood.item() + prior
        return ll + prior


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

        for centerPos in xrange(start + 1, end):

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


        pos = end-1
        currentCol = self.fwdScore[end-1]
        bestConfig = max(currentCol, key=lambda x:currentCol[x])

        while True:
            if cogBase(bestConfig) != self.sequence[pos]:
                # Found a modification - save it!
                modCalls[pos] = cogBase(bestConfig)

            bestConfig = self.fwdPrevState[pos][bestConfig]
            pos -= 1

            if bestConfig is None:
                break


        if self.useLDA:
            # allow LDA-predicted sites through to GFF file
            for pos in range(start, end):
                if self.rawKinetics.has_key(pos):
                    if self.rawKinetics[pos].has_key('Ca5C'):
                        if 'K' in self.modsToCall:
                            # cutoff = min( 0, self.rawKinetics[pos]['coverage']/20.0 - 3.0 )
                            cutoff = 0
                            echoSites = [modCalls[pos + i] for i in [-6, -5, -2, -1, 2] if modCalls.has_key(pos + i)]
                        else:
                            cutoff = -2.25
                            echoSites = [modCalls[pos + i] for i in range(-10,11) if modCalls.has_key(pos + i)]
                        if self.rawKinetics[pos]['Ca5C'] < cutoff:
                            # so long as those sites are not in the vicinity of a m6A/m4C call
                            if 'H' not in echoSites and 'J' not in echoSites:
                                modCalls[pos] = 'K'
                        # else:
                        #     # remove any non-LDA-predicted sites from modCalls dictionary?
                        #     if modCalls.has_key(pos):
                        #         del modCalls[pos]


          
        # correct adjacent calls:
        if self.useLDA:
            for pos in range(start + 2, end - 2, 2 ):
                x = [pos + i for i in range(-2,3) if modCalls.has_key( pos + i) and self.rawKinetics.has_key( pos + i) ]
                y = [modCalls[j] for j in x]
                if y.count('K') > 1:
                    tmp = [self.rawKinetics[j]['Ca5C'] for j in x if self.rawKinetics[j].has_key('Ca5C') ]
                    if len(tmp) > 0:
                        lowest = min(tmp)
                        for j in x:
                            if self.rawKinetics[j].has_key('Ca5C'):
                                if self.rawKinetics[j]['Ca5C'] > lowest:
                                    del modCalls[j]
                    


            # if adjacent m5C calls are made by the LDA, select the one that has the lower LDA score (Ca5C)
            # for pos in range(start, end):
            #     if modCalls.has_key(pos) and self.rawKinetics.has_key(pos) and modCalls.has_key(pos+1) and self.rawKinetics.has_key(pos+1):
            #         if self.rawKinetics[pos].has_key('Ca5C') and self.rawKinetics[pos+1].has_key('Ca5C'):
            #             if self.rawKinetics[pos]['Ca5C'] < self.rawKinetics[pos+1]['Ca5C']:
            #                 del modCalls[pos+1] 
            #             else:
            #                 del modCalls[pos] 

            
        return modCalls


    def scoreMods(self, modCalls):
        """
        For each modification in the best scoring configuration, score a config excluding the current mod against the winning config
        use this value as the Qmod for the deleted modification
        """

        qvModCalls = dict()

        modSeq = a.array('c')
        modSeq.fromstring(self.sequence)


        # Apply the found modifications to the raw sequence
        for (pos, call) in modCalls.items():
            modSeq[pos] = call


        for (pos, call) in modCalls.items():

            # Score the modified template at all positions affected by this mod
            modScore = self.scoreRegion(pos - self.post, pos + self.pre, modSeq)
            modScores = self.getRegionScores(pos - self.post, pos + self.pre, modSeq)

            if self.methylFractionFlag and self.rawKinetics.has_key(pos):
                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:
                    modifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq)		

            # Switch back to the unmodified base and re-score
            modSeq[pos] = canonicalBaseMap[call]
            noModScore = self.scoreRegion(pos - self.post, pos + self.pre, modSeq)
            noModScores = self.getRegionScores(pos - self.post, pos + self.pre, modSeq)

            if self.methylFractionFlag and self.rawKinetics.has_key(pos):
                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:
                    unModifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq) 	

            # Put back the modified base
            modSeq[pos] = call

            # Compute score difference
            llr = modScore - noModScore

            # Convert from LLR to phred-scaled probability of modification
            qModScore = 10 * llr * log10e + 10*log1p(exp(-llr))*log10e


            # Figure out which secondary peaks were likely generated by this modification
            # What is the posterior that the peak was generated by this mod?
            maskPos = self.findMaskPositions(pos, modScores, noModScores)


            # FIXME:  Without this, currently, the identificationQv score is too low for many Ca5C sites
            if self.useLDA:
                if self.rawKinetics.has_key(pos):
                    if self.rawKinetics[pos].has_key('Ca5C'):
                        llr = -self.rawKinetics[pos]['Ca5C']
                        qModScore = 100 * llr * log10e + 100*log1p(exp(-llr))*log10e



            if self.methylFractionFlag and self.rawKinetics.has_key(pos):
                if self.rawKinetics[pos]["coverage"] > self.methylMinCov:
                    # Use modifiedMeanVectors and unmodifiedMeanVectors to calculate mixing proportion, and 95% CI limits.
                    methylFracEst,methylFracLow,methylFracUpp = self.estimateMethylatedFractions(pos, unModifiedMeanVectors, modifiedMeanVectors, ModificationPeakMask[modNames[call]] )
                    qvModCalls[pos] = { 'modification' : modNames[call], 'QMod' : qModScore, 'LLR' : llr, 'Mask': maskPos, \
                                    FRAC: methylFracEst, FRAClow: methylFracLow, FRACup: methylFracUpp }

                else:
                    qvModCalls[pos] = { 'modification' : modNames[call], 'QMod' : qModScore, 'LLR' : llr, 'Mask': maskPos }

            else:
                # Store the full results
                qvModCalls[pos] = { 'modification' : modNames[call], 'QMod' : qModScore, 'LLR' : llr, 'Mask': maskPos }

        return qvModCalls


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


    # Return value of mixture model log likelihood function

    def mixModelFn(self, p, a0, a1):
        tmp = (1-p)*a0 + p*a1
        return -np.log( tmp[ np.nonzero(tmp) ] ).sum()
        # return -np.ma.log( tmp ).sum()


    # Try to speed up calculation by avoiding a call to scipy.stats.norm.pdf() 

    def replaceScipyNormPdf( self, data, mu ):
        # return np.exp( -np.divide( data, mu) ) / mu
        tmp = np.divide( data, mu )
        return np.exp( np.subtract( tmp, np.power(tmp, 2) / 2.0 ) ) / mu 
        # pdf for normal distribution: res = res / sqrt( 2 * pi ) (can factor out sqrt(2 * pi))


    # Return optimum argument (mixing proportion) of mixture model log likelihood function.

    def estimateSingleFraction(self, mu1, data, mu0, L ):
        a0 = self.replaceScipyNormPdf( data, mu0 )
        a1 = self.replaceScipyNormPdf( data, mu1 )
        # if f'(0) < 0 (equ. a1/a0 < L), then f'(1) < 0 as well and solution p-hat <= 0
        if np.divide(a1, a0).sum() <= L:
            return 0.0
        # if f'(1) > 0 (equ. a0/a1 < L), then f'(0) > 0 as well and solution p-hat >= 1
        if np.divide(a0, a1).sum() <= L:
            return 1.0
        # unconstrained minimization of convex, single-variable function
        res = fminbound(self.mixModelFn, 0.01, 0.99, args=(a0, a1), xtol=1e-02)
        return res



    # Try bias-corrected, accelerated quantiles for bootstrap confidence intervals

    def bcaQuantile( self, estimate, bootDist, data, mu0, mu1, nSamples, n ):

        tmp = sum( y <= estimate for y in bootDist ) / float(nSamples + 1)
        if tmp > 0 and tmp < 1:

            # bias correction
            z0 = s.norm.ppf( tmp )

            # acceleration
            x = np.zeros(n)
            for i in range(n):
                x[i] = self.estimateSingleFraction(mu1, np.delete(data, i), mu0, n-1)
            xbar = np.mean(x)
            denom =  np.power( np.sum( np.power( x - xbar, 2) ), 1.5 )
            if abs(denom) < 1e-4:
                q1 = 2.5
                q2 = 97.5
            else:
                a = np.divide( np.sum( np.power( x - xbar, 3) ),  denom ) / 6.0

                # quantiles: (k1 and k2 are defined globally)
                q1 = 100*s.norm.cdf( z0 + (z0 + k1)/(1 - a*(z0 + k1)) )
                q2 = 100*s.norm.cdf( z0 + (z0 + k2)/(1 - a*(z0 + k2)) )

        elif tmp == 0.0:
            q1 = 0
            q2 = 0

        elif tmp == 1.0:
            q1 = 100
            q2 = 100

        return (q1, q2)


    # Bootstraps mix prop estimates to return estimate and simple bounds for 95% confidence interval

    def bootstrap(self, pos, mu0, mu1, nSamples = 500):

        if not self.rawKinetics.has_key( pos ):
            return np.array( [ float('nan'), float('nan'), float('nan') ] )

        res = np.zeros(3)
        sample = self.rawKinetics[pos]["rawData"]
        L = len(sample)
        X = np.zeros(nSamples+1)
        res[0] = self.estimateSingleFraction(mu1, sample, mu0, L)
        X[nSamples] = res[0]

        for i in range(nSamples):
            bootstrappedSamples = sample[s.randint.rvs(0, L-1, size=L)]
            X[i] = self.estimateSingleFraction(mu1, bootstrappedSamples, mu0, L)

        q1,q2 = self.bcaQuantile( res[0], X, sample, mu0, mu1, (nSamples+1), L )
        res[1] = np.percentile(X, q1)
        res[2] = np.percentile(X, q2)
        return res


    # Returns [estimate, 95% CI lower bnd, 95% CI upper bound] using a weighted sum
    # The hope is that this would work better for a multi-site signature, such as m5C_TET

    def estimateMethylatedFractions(self, pos, meanVector, modMeanVector, maskPos ):

        maskPos = np.array(maskPos)
        L = len(maskPos)
        if L == 0:
            res = self.bootstrap(pos, meanVector[self.post], modMeanVector[self.post] )
        else:
            est = np.zeros(L)
            low = np.zeros(L)
            upp = np.zeros(L)
            res = np.zeros(3)
            wts = np.zeros(L)

            # for offset in maskPos:
            for count in range(L):
                offset = maskPos[count]
                mu0 = meanVector[ self.post + offset ]
                mu1 = modMeanVector[ self.post + offset ]
                if mu1 > mu0:
                    k = self.bootstrap( (pos + offset), mu0, mu1 )
                    wts[count] = k[0] * (mu1 - mu0)
                    est[count] = k[0]
                    low[count] = k[1]
                    upp[count] = k[2]

            if sum(wts) > 1e-3:
                wts = wts/sum(wts)
                res[0] = np.multiply(est, wts).sum()
                res[1] = np.multiply(low, wts).sum()
                res[2] = np.multiply(upp, wts).sum()

        # print str(res)
        return res



    def scoreRegion(self, start, end, sequence):

        sc = 0.0
        for pos in xrange(start, end+1):
            ctx = sequence[(pos-self.pre):(pos + self.post + 1)].tostring()
            if self.scores.has_key(pos):
                sc += self.scores[pos][ctx]

        return sc


    def getRegionScores(self, start, end, sequence):
        scores = np.zeros(end-start+1)

        for pos in xrange(start, end+1):
            ctx = sequence[(pos-self.pre):(pos + self.post + 1)].tostring()
            if self.scores.has_key(pos):
                scores[pos - start] = self.scores[pos][ctx]

        return scores

    def findMaskPositions(self, pos, modScores, noModScores):

        maskPos = []
        start = pos - self.post
        end = pos + self.pre

        for i in xrange(start, end + 1):
            # Add a neighboring peak to the mask if
            # a) it has a single-site qv > 20
            # b) the observed IPDs are somewhat more likely under the modified hypothesis than the unmodified hypothesis
            if self.rawKinetics.has_key(i) and self.rawKinetics[i]["score"] > 20:
                if modScores[i - start] - noModScores[i - start] > 1.0:
                    maskPos.append(i - pos)

        return maskPos


    def compareStates(self, current, prev):
        return current[0:-1] == prev[1:]



    # Everything below here is unused for now:


    # Return second derivative of mixture model log likelihood function - unused for now
    def mixModelFnPrime2(self, p, a0, a1):
        tmp = np.square( (1-p)*a0 + p*a1 )
        nonzero_indices = np.nonzero(tmp)
        return np.divide( np.square(a1 - a0)[nonzero_indices], tmp[nonzero_indices] ).sum()


    # Return third derivative of mixture model log likelihood function - unused for now
    def mixModelFnPrime3(self, p, a0, a1):
        tmp = np.power( (1-p)*a0 + p*a1, 3 )
        nonzero_indices = np.nonzero(tmp)
        return -np.divide( np.power(a1 - a0, 3)[nonzero_indices], tmp[nonzero_indices] ).sum()


    # Try removing very large values before case resampling for bootstrap estimation - unused for now
    def processSample( self, sample ):
        q1 = np.percentile(sample, 25)
        q2 = np.percentile(sample, 75)
        iqr = 1.5*(q2 - q1)
        uif = q2 + iqr
        lif = q1 - iqr
        def removeBoxplotOutliers(x):
            if (x > lif) and (x < uif):
                return x
        return filter(removeBoxplotOutliers, sample)


    # Return derivative of mixture model log likelihood function -- unused for now
    def mixModelFnPrime(self, p, a0, a1):
        tmp = (1-p)*a0 + p*a1
        nonzero_indices = np.nonzero(tmp)
        return -np.divide( (a1 - a0)[nonzero_indices], tmp[nonzero_indices] ).sum()


    # unconstrained minimization of convex, single-variable function - unused for now
    # much slower than fminbound
    def homeMadeMinimization( self, a0, a1, low, up, xtol = 1e-02, maxIters = 500 ):
        nIters = 0
        while (up - low) > xtol and nIters < maxIters:
            p0 = (up - low)/2.0
            if self.mixModelFnPrime(p0, a0, a1) <= 0:
                low = p0
            else:
                up = p0
            nIters += 1
        return p0


