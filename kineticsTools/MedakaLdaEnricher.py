from __future__ import print_function
from __future__ import absolute_import
# Try to implement method used in Morishita et al.'s Medaka fish genome paper here

from collections import defaultdict, Counter

import os
from math import sqrt
import math
import scipy.stats as s
import array as a

from scipy.optimize import fminbound
from scipy.special import gammaln as gamln
from numpy import log, pi, log10, e, log1p, exp
import numpy as np

from .MultiSiteCommon import MultiSiteCommon


class MedakaLdaEnricher(MultiSiteCommon):

    def __init__(self, gbmModel, sequence, rawKinetics, m5Cclassifier):

        MultiSiteCommon.__init__(self, gbmModel, sequence, rawKinetics)

        models = np.genfromtxt(m5Cclassifier, delimiter=',' )
        self.fwd_model = models[:,0]
        self.rev_model = models[:,1]



    # write a method to take perSiteResults dictionary in and add a column Ca5C
    def useLDAmodel(self, kinetics, pos, model, up, down ):
        """ Test out LDA model """

        print("From use LDA model.\n")

        res = np.zeros((up + down + 1, 6))
        ind = 0

        # range from -down to +up
        for offset in range(-down, (up + 1)):
            a = pos + offset

            std = kinetics[a]["tErr"] * sqrt( kinetics[a]["coverage"] )
            mErr = 0.01 + 0.03 * kinetics[a]["modelPrediction"] + 0.06 * kinetics[a]["modelPrediction"] ** 2
            den = sqrt( mErr **2 + std **2 )
            t0 = ( kinetics[a]["tMean"] - kinetics[a]["modelPrediction"] ) / den

            res[ind, ] = [kinetics[a]["tMean"], kinetics[a]["modelPrediction"], std, np.exp( t0 ) - 0.01, kinetics[a]["ipdRatio"], den]
            ind += 1

        predictors = np.hstack(np.log(res + 0.01).transpose())
        tmp = sum( np.multiply( predictors, model[1:] ))  + model[0]

        return tmp


    def callLDAstrand(self, kinetics, strand, model, up, down):

        print("From callLDAstrand.\n")
      
        tmp = [d for d in kinetics if d["strand"] == strand]
        tmp.sort(key=lambda x: x["tpl"])

        L = len(tmp)
        for pos in range(down, (L - up)):
            if (strand == 0 and tmp[pos]["base"] == 'C' and tmp[pos+1]["base"] == 'G'):
                tmp[pos]["Ca5C"] = self.useLDAmodel(tmp, pos, model, up, down ) 

            if (strand == 1 and tmp[pos]["base"] == 'C' and tmp[pos-1]["base"] == 'G'):
                tmp[pos-1]["Ca5C"] = self.useLDAmodel(tmp, pos, model, up, down )

        return tmp


    def aggregate(self, dataset, group_by_key, sum_value_key):

        print("From aggregate.\n")
        emp = {}
        for item in dataset:
            if sum_value_key in item:
                if item[group_by_key] in emp:
                    emp[ item[group_by_key] ] += item[sum_value_key]
                else:
                    emp[ item[group_by_key] ] = item[sum_value_key]

        # Need to go back over the set again?
        for item in dataset:
            if sum_value_key in item:
                item[ sum_value_key ] = emp[ item[group_by_key] ]

        return dataset



    def callEnricherFunction(self, kinetics, up=10, down=10):

        print("From callEnricher function.\n")

        fwd = self.callLDAstrand(kinetics, 0, self.fwd_model, up, down) 
        rev = self.callLDAstrand(kinetics, 1, self.rev_model, up, down)
        res = fwd + rev
        res.sort( key = lambda x: x["tpl"] )

        # Would like to (1) find rows where tpl is the same, (2) add the Ca5C columns for these rows, 
        # (3) store the same result under Ca5C for each row.
        # In R, this would be:  a <- aggregate(Ca5C ~ tpl, df, sum )

        res = self.aggregate( res, 'tpl', 'Ca5C')
        return res
