# pylint: skip-file
# FIXME this is currently non-functional and not used anywhere, should we
# remove it?

from math import sqrt
import math
import scipy.stats as s
import array as a

from scipy.optimize import fminbound
from scipy.special import gammaln as gamln
from numpy import *
from numpy import log, pi, log10, e, log1p, exp
import numpy as np

from kineticsTools.MultiSiteCommon import MultiSiteCommon
from kineticsTools.MixtureEstimationMethods import MixtureEstimationMethods


class PositiveControlEnricher(MultiSiteCommon):

    def __init__(self, gbmModel, sequence, rawKinetics):

        MultiSiteCommon.__init__(self, gbmModel, sequence, rawKinetics)
        self.fwd_model = np.genfromtxt(
            "/home/UNIXHOME/obanerjee/initial_lr_model_weights_fwd.csv", delimiter=',')
        self.rev_model = np.genfromtxt(
            "/home/UNIXHOME/obanerjee/initial_lr_model_weights_rev.csv", delimiter=',')
        self.fwd_model = np.squeeze(np.asarray(self.fwd_model))
        self.rev_model = np.squeeze(np.asarray(self.rev_model))

    def fn(self, l):
        if l == "A":
            return 1
        if l == "C":
            return 2
        if l == "G":
            return 3
        return 4

    def tStatisticDenominator(self, mu0, tErr):

        em = 0.01 + 0.03 * mu0 + 0.06 * mu0 ** 1.7
        den = sqrt(em ** 2 + tErr ** 2)
        return den

    def applyLRmodel(self, kinetics, pos, unmodIPDs, modifIPDs, model, up, down, context):
        """ Test out LDA model """

        res = np.zeros((up + down + 1, 7))
        ind = 0

        # range from -down to +up
        for offset in range(-down, (up + 1)):
            a = pos + offset
            tmp = np.squeeze(np.asarray([kinetics[a]["tMean"], kinetics[a]["tErr"], kinetics[a]
                                         ["coverage"], unmodIPDs[offset + down], modifIPDs[offset + down]]))

            # get t-statistics corresponding to mu0 and mu1:
            den = self.tStatisticDenominator(tmp[3], tmp[1])
            k1 = -(tmp[0] - tmp[3]) / den
            k2 = -(tmp[0] - tmp[4]) / den
            tmp = np.append(np.log(tmp + 0.01), k1)
            tmp = np.append(tmp, k2)

            res[ind, ] = tmp
            ind += 1

        # collected features for prediction:
        apply = np.hstack(res.transpose())
        apply = np.squeeze(np.asarray(apply))
        apply[isnan(apply)] = 0

        # include context:
        context = context[(pos - 10):(pos + 11)]
        del context[11]
        context = np.array(map(self.fn, context))
        apply = np.concatenate([apply, context])

        # calculate logistic regression score:
        z = sum(np.multiply(apply, model[1:])) + model[0]
        score = -z - np.log(1 + np.exp(-z))

        return score

    def callLRstrand(self, kinetics, strand, model, up, down):

        tmp = [d for d in kinetics if d["strand"] == strand]
        tmp.sort(key=lambda x: x["tpl"])

        modSeq = a.array('c')
        modSeq.fromstring(self.sequence)

        L = len(tmp)

        for pos in range(down, (L - up)):

            if tmp[pos]["base"] == 'C':

                # Get H1 means:
                modSeq[pos] = 'I'
                modifIPDs = self.getContextMeans(pos - 10, pos + 10, modSeq)

                # Get H0 means:
                modSeq[pos] = 'C'
                unmodIPDs = self.getContextMeans(pos - 10, pos + 10, modSeq)

                tmp[pos]["Ca5C"] = self.applyLRmodel(
                    tmp, pos, unmodIPDs, modifIPDs, model, up, down, modSeq)

        return tmp

    def callEnricherFunction(self, kinetics, up=10, down=10):

        # Compute all the required mean ipds under all possible composite hypotheses
        self.computeContextMeans()

        fwd = self.callLRstrand(kinetics, 0, self.fwd_model, up, down)
        rev = self.callLRstrand(kinetics, 1, self.rev_model, up, down)
        res = fwd + rev
        res.sort(key=lambda x: x["tpl"])
        return res

    def scoreMods(self, modCalls):
        """
        For each modification in the best scoring configuration, score a config excluding the current mod against the winning config
        use this value as the Qmod for the deleted modification
        """

        qvModCalls = dict()

        modSeq = a.array('c')
        modSeq.fromstring(self.sequence)

        modCalls = [i for i in range(len(modSeq)) if modSeq[i] == 'C']

        for pos in modCalls:

            # now apply the modification at this position:
            modSeq[pos] = 'K'
            modifIPDs = self.getContextMeans(
                pos - self.post, pos + self.pre, modSeq)

            # using canonical base map, try to get H0 means:
            modSeq[pos] = canonicalBaseMap[call]
            unmodIPDs = self.getContextMeans(
                pos - self.post, pos + self.pre, modSeq)

            # try to collect related statistics:  tMean, tErr, tStatistic
            tmp = self.applyLRmodel(
                kinetics, pos, unmodIPDs, modifIPDs, modSeq, up, down)

            # now try to use these vectors to make a basic decision:
            basicDecision[pos] = {'score': tmp}

        return basicDecision
