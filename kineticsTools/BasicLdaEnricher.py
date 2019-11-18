# Basic LDA Enricher class

from math import sqrt
import math
import scipy.stats as s
import array as a

from scipy.optimize import fminbound
from scipy.special import gammaln as gamln
from numpy import log, pi, log10, e, log1p, exp
import numpy as np

from .MultiSiteCommon import MultiSiteCommon
from .MixtureEstimationMethods import MixtureEstimationMethods


class BasicLdaEnricher(MultiSiteCommon):

    def __init__(self, gbmModel, sequence, rawKinetics, identifyFlag, modsToCall=['H', 'J', 'K']):

        MultiSiteCommon.__init__(self, gbmModel, sequence, rawKinetics)

        # FIXME: For debugging LDA, load in parameters for forward and reverse strands:

        self.fwd_model = np.genfromtxt(
            "/home/UNIXHOME/obanerjee/nat_fwd_model_expanded.csv", delimiter=',')
        self.rev_model = np.genfromtxt(
            "/home/UNIXHOME/obanerjee/nat_rev_model_expanded.csv", delimiter=',')

        if identifyFlag:
            if 'K' in modsToCall:
                self.fwd_model = np.genfromtxt(
                    "/home/UNIXHOME/obanerjee/tet_fwd_model_expanded.csv", delimiter=',')
                self.rev_model = np.genfromtxt(
                    "/home/UNIXHOME/obanerjee/tet_rev_model_expanded.csv", delimiter=',')

    # write a method to take perSiteResults dictionary in and add a column Ca5C
    def useLDAmodel(self, kinetics, pos, model, up, down):
        """ Test out LDA model """

        res = np.zeros((up + down + 1, 5))
        ind = 0

        # range from -down to +up
        for offset in range(-down, (up + 1)):
            a = pos + offset
            # res[ind,] = [kinetics[a]["tMean"], kinetics[a]["modelPrediction"], kinetics[a]["tErr"], kinetics[a]["coverage"]]
            res[ind, ] = [kinetics[a]["tMean"], kinetics[a]["modelPrediction"], kinetics[a]
                          ["tErr"], kinetics[a]["coverage"], np.exp(kinetics[a]["tStatistic"]) - 0.01]
            ind += 1

        apply = np.hstack(np.log(res + 0.01).transpose())
        tmp = sum(np.multiply(apply, model[1:])) + model[0]
        return tmp

    def callLDAstrand(self, kinetics, strand, model, up, down):
        tmp = [d for d in kinetics if d["strand"] == strand]
        tmp.sort(key=lambda x: x["tpl"])

        L = len(tmp)
        for pos in range(down, (L - up)):
            if tmp[pos]["base"] == 'C':
                tmp[pos]["Ca5C"] = self.useLDAmodel(tmp, pos, model, up, down)

        return tmp

    def callEnricherFunction(self, kinetics, up=10, down=10):

        fwd = self.callLDAstrand(kinetics, 0, self.fwd_model, up, down)
        rev = self.callLDAstrand(kinetics, 1, self.rev_model, up, down)
        res = fwd + rev
        res.sort(key=lambda x: x["tpl"])
        return res
