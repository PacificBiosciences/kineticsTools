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
from numpy import log, pi, log10, e, log1p, exp
import numpy as np


log10e = log10(e)

# canonicalBaseMap = { 'A': 'A', 'C':'C', 'G':'G', 'T':'T', 'H':'A', 'I':'C', 'J':'C', 'K':'C' }
# modNames = { 'H':'m6A', 'I':'m5C', 'J':'m4C', 'K':'m5C' }

# ModificationPeakMask = { 'm6A' : [0, -5], 'm4C': [0, -5], 'm5C': [2, 0, -1, -2, -4, -5, -6]  }

# Labels for modified fraction:

# FRAC = 'frac'
# FRAClow = 'fracLow'
# FRACup = 'fracUp'

# Try computing these only once

k1 = s.norm.ppf(0.025)
k2 = s.norm.ppf(0.975)


class MixtureEstimationMethods(object):

    def __init__(self, gbmModelPost, gbmModelPre, rawKinetics, methylMinCov):
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
        # self.useLDA = useLDAFlag
        # self.modsToCall = modsToCall
        # self.methylFractionFlag = methylFractionFlag

        # log1p = math.log(0.05)
        # self.modPriors = { 'H': log1p, 'I': log1p, 'J': log1p, 'K': log1p }

        # self.gbmModel = gbmModel
        # self.sequence = sequence

        # self.callStart = callBounds[0]
        # self.callEnd = callBounds[1]

        # Extents that we will attemp to call a modification
        # self.callRange = xrange(self.callStart, self.callEnd)

        # These switch because we changing viewpoints
        self.pre = gbmModelPost
        self.post = gbmModelPre

        # self.lStart = self.pre
        # self.lEnd = len(self.sequence) - self.post

        # Extents that we will use for likelihoods
        # self.likelihoodRange = xrange(self.lStart, self.lEnd)
        # self.alternateBases = dict((x, set(sequence[x])) for x in xrange(len(sequence)))

        self.rawKinetics = rawKinetics

    # Return value of mixture model log likelihood function
    def mixModelFn(self, p, a0, a1):
        tmp = (1 - p) * a0 + p * a1
        return -np.log(tmp[np.nonzero(tmp)]).sum()
        # return -np.ma.log( tmp ).sum()

    # Try to speed up calculation by avoiding a call to scipy.stats.norm.pdf()
    def replaceScipyNormPdf(self, data, mu):
        # return np.exp( -np.divide( data, mu) ) / mu
        tmp = np.divide(data, mu)
        return np.exp(np.subtract(tmp, np.power(tmp, 2) / 2.0)) / mu
        # pdf for normal distribution: res = res / sqrt( 2 * pi ) (can factor out sqrt(2 * pi))

    # Return optimum argument (mixing proportion) of mixture model log likelihood function.
    def estimateSingleFraction(self, mu1, data, mu0, L):
        a0 = self.replaceScipyNormPdf(data, mu0)
        a1 = self.replaceScipyNormPdf(data, mu1)
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
    def bcaQuantile(self, estimate, bootDist, data, mu0, mu1, nSamples, n):

        tmp = sum(y <= estimate for y in bootDist) / float(nSamples + 1)
        if tmp > 0 and tmp < 1:

            # bias correction
            z0 = s.norm.ppf(tmp)

            # acceleration
            x = np.zeros(n)
            for i in range(n):
                x[i] = self.estimateSingleFraction(mu1, np.delete(data, i), mu0, n - 1)
            xbar = np.mean(x)
            denom = np.power(np.sum(np.power(x - xbar, 2)), 1.5)
            if abs(denom) < 1e-4:
                q1 = 2.5
                q2 = 97.5
            else:
                a = np.divide(np.sum(np.power(x - xbar, 3)), denom) / 6.0

                # quantiles: (k1 and k2 are defined globally)
                q1 = 100 * s.norm.cdf(z0 + (z0 + k1) / (1 - a * (z0 + k1)))
                q2 = 100 * s.norm.cdf(z0 + (z0 + k2) / (1 - a * (z0 + k2)))

        elif tmp == 0.0:
            q1 = 0
            q2 = 0

        elif tmp == 1.0:
            q1 = 100
            q2 = 100

        return (q1, q2)

    # Bootstraps mix prop estimates to return estimate and simple bounds for 95% confidence interval
    def bootstrap(self, pos, mu0, mu1, nSamples=500):

        if not self.rawKinetics.has_key(pos):
            return np.array([float('nan'), float('nan'), float('nan')])

        res = np.zeros(3)
        sample = self.rawKinetics[pos]["rawData"]
        L = len(sample)
        X = np.zeros(nSamples + 1)
        res[0] = self.estimateSingleFraction(mu1, sample, mu0, L)
        X[nSamples] = res[0]

        for i in range(nSamples):
            bootstrappedSamples = sample[s.randint.rvs(0, L - 1, size=L)]
            X[i] = self.estimateSingleFraction(mu1, bootstrappedSamples, mu0, L)

        q1, q2 = self.bcaQuantile(res[0], X, sample, mu0, mu1, (nSamples + 1), L)
        res[1] = np.percentile(X, q1)
        res[2] = np.percentile(X, q2)
        return res

    # Returns [estimate, 95% CI lower bnd, 95% CI upper bound] using a weighted sum
    # The hope is that this would work better for a multi-site signature, such as m5C_TET
    def estimateMethylatedFractions(self, pos, meanVector, modMeanVector, maskPos):

        maskPos = np.array(maskPos)
        L = len(maskPos)
        if L == 0:
            res = self.bootstrap(pos, meanVector[self.post], modMeanVector[self.post])
        else:
            est = np.zeros(L)
            low = np.zeros(L)
            upp = np.zeros(L)
            res = np.zeros(3)
            wts = np.zeros(L)

            # for offset in maskPos:
            for count in range(L):
                offset = maskPos[count]
                mu0 = meanVector[self.post + offset]
                mu1 = modMeanVector[self.post + offset]
                if mu1 > mu0:
                    k = self.bootstrap((pos + offset), mu0, mu1)
                    wts[count] = k[0] * (mu1 - mu0)
                    est[count] = k[0]
                    low[count] = k[1]
                    upp[count] = k[2]

            if sum(wts) > 1e-3:
                wts = wts / sum(wts)
                res[0] = np.multiply(est, wts).sum()
                res[1] = np.multiply(low, wts).sum()
                res[2] = np.multiply(upp, wts).sum()

        # print str(res)
        return res

    # Return the optimal mixing proportion in the detection case: estimate both p and mu1
    def optimalMixProportion(self, data, mu0, L):
        mu1 = fminbound(self.estimateSingleFraction, mu0, 10.0 * mu0, args=(data, mu0, L), xtol=1e-01)
        res = self.estimateSingleFraction(mu1, data, mu0, L)
        return res

    # Bootstraps mix prop estimates to return estimate and simple bounds for 95% confidence interval
    def detectionMixModelBootstrap(self, modelPrediction, data, nSamples=100):

        # Case-resampled bootstrapped estimates:
        L = len(data)
        res = np.zeros(4)
        res[0] = self.optimalMixProportion(data, modelPrediction, L)
        X = np.zeros(nSamples + 1)
        X[nSamples] = res[0]
        for i in range(nSamples):
            resampledData = [data[j] for j in s.randint.rvs(0, L - 1, size=L)]
            X[i] = self.optimalMixProportion(resampledData, modelPrediction, L)

        # A very basic way to estimate the 95% confidence interval:
        res[1] = np.percentile(X, 2.5)
        res[2] = np.percentile(X, 97.5)

        # Estimate a weight:
        # weight = np.maximum( (x[1] - modelPrediction), 0 )
        res[3] = 1.0
        return res

    # Everything below here is unused for now:
    # Return second derivative of mixture model log likelihood function - unused for now
    def mixModelFnPrime2(self, p, a0, a1):
        tmp = np.square((1 - p) * a0 + p * a1)
        nonzero_indices = np.nonzero(tmp)
        return np.divide(np.square(a1 - a0)[nonzero_indices], tmp[nonzero_indices]).sum()

    # Return third derivative of mixture model log likelihood function - unused for now
    def mixModelFnPrime3(self, p, a0, a1):
        tmp = np.power((1 - p) * a0 + p * a1, 3)
        nonzero_indices = np.nonzero(tmp)
        return -np.divide(np.power(a1 - a0, 3)[nonzero_indices], tmp[nonzero_indices]).sum()

    # Try removing very large values before case resampling for bootstrap estimation - unused for now
    def processSample(self, sample):
        q1 = np.percentile(sample, 25)
        q2 = np.percentile(sample, 75)
        iqr = 1.5 * (q2 - q1)
        uif = q2 + iqr
        lif = q1 - iqr

        def removeBoxplotOutliers(x):
            if (x > lif) and (x < uif):
                return x
        return filter(removeBoxplotOutliers, sample)

    # Return derivative of mixture model log likelihood function -- unused for now
    def mixModelFnPrime(self, p, a0, a1):
        tmp = (1 - p) * a0 + p * a1
        nonzero_indices = np.nonzero(tmp)
        return -np.divide((a1 - a0)[nonzero_indices], tmp[nonzero_indices]).sum()

    # unconstrained minimization of convex, single-variable function - unused for now
    # much slower than fminbound
    def homeMadeMinimization(self, a0, a1, low, up, xtol=1e-02, maxIters=500):
        nIters = 0
        while (up - low) > xtol and nIters < maxIters:
            p0 = (up - low) / 2.0
            if self.mixModelFnPrime(p0, a0, a1) <= 0:
                low = p0
            else:
                up = p0
            nIters += 1
        return p0
