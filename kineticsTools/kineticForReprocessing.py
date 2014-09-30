from itertools import ifilter
import json
from scipy.optimize import fminbound
import array as a
from math import sqrt
import math
from scipy.special import erfc
import logging

# from scipy.optimize import fmin_tnc	# Conjugate gradient
import scipy.stats as s
import numpy as np
import scipy.stats.mstats as mstats
import sys

from WorkerProcess import WorkerProcess, WorkerThread

from pbcore.io import GffReader

# from ModificationDecode?
canonicalBaseMap = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'H': 'A', 'I': 'C', 'J': 'C', 'K': 'C'}


toMod = {'m6A': 'H', 'm5C': 'I', 'm4C': 'J'}
modNames = {'H': 'm6A', 'I': 'm5C', 'J': 'm4C', 'K': 'm5C'}

ModificationPeakMask = {'m6A': [0, -5], 'm4C': [0, -5], 'm5C': [2, 0, -1, -2, -4, -5, -6]}

# Raw ipd record
ipdRec = [('tpl', '<u4'), ('strand', '<i8'), ('ipd', '<f4')]


# Labels for modified fraction:
FRAC = 'frac'
FRAClow = 'fracLow'
FRACup = 'fracUp'

k1 = s.norm.ppf(0.025)
k2 = s.norm.ppf(0.975)


def dictsToRecArray(dictList):
    """
    Convert a list of identical dictionaries to a recArray (not used?)
    """

    d1 = dictList[0]
    n = len(dictList)

    colNames = d1.keys()
    types = [np.array([d1[x]]).dtype for x in colNames]
    recArrayType = zip(colNames, types)

    res = np.zeros(n, dtype=recArrayType)

    for i in xrange(n):
        d = dictList[i]
        for c in colNames:
            res[c][i] = d[c]

    return res


class HitSelection:

    """
    Static methods for selecting sets of reads to use
    """

    @staticmethod
    def readFilter(alnInfo, minAcc=0.8, minLength=50):
        """Map the alnInfo recArray to bools indicating whether the hit passes a quality filter"""
        rl = np.abs(alnInfo.tEnd - alnInfo.tStart)
        acc = 1.0 - (alnInfo.nIns + alnInfo.nDel + alnInfo.nMM) / rl
        return np.logical_and(acc > minAcc, rl > minLength)


class KineticReprocessWorker(object):

    """
    Manages the summarization of pulse features over a single reference
    """

    def __init__(self, ipdModel):
        self.ipdModel = ipdModel
        self.debug = False

    def _prepForReferenceWindow(self, referenceWindow):
        """
        Helper function for testing
        """
        (reference, start, end) = referenceWindow
        self.refId = reference

        # Each chunk is from a single reference -- fire up meanIpd func on the current reference
        self.meanIpdFunc = self.ipdModel.predictIpdFunc(reference)

        # Get the cognate base at a given position
        self.cognateBaseFunc = self.ipdModel.cognateBaseFunc(reference)

    def getReferenceIndex(self, d, oldData):

        if oldData:
            # d["seqID"] = ref000001, for example
            s = int(d["seqID"].split("ref")[1]) - 1

        else:
            # match the 1st 64 characters (at most) of the ref name
            # tmp = d["seqID"].__len__()
            # if tmp > 64:
            #     s = self.selectRef.index( d["seqID"][:64] )
            # else:
            #     s = self.selectRef.index(d["seqID"])
            s = self.selectRef.index(d["seqID"])

        return s

    def fillOutEasyInformation(self, d, motifInfo):

        k = {}
        k['refName'] = d["seqID"]
        k['refId'] = self.refId
        k['tpl'] = d["pos"] - 1
        k['score'] = float(d["score"])
        k['modification'] = d["type"]

        if d["strand"] == '+':
            k['strand'] = 0
        else:
            k['strand'] = 1

        u = d["attributes"]
        k['coverage'] = u["coverage"]
        k['ipdRatio'] = u["IPDRatio"]

        if "identificationQv" in u:
            k['modificationScore'] = u["identificationQv"]

        if "motif" in u and u["motif"] in motifInfo.keys():
            self.motif = u["motif"]
            k['motif'] = self.motif
            k['id'] = self.motif

        # FIXME: Placeholders for these fields in output?

        k[FRAC] = np.nan
        k[FRAClow] = np.nan
        k[FRACup] = np.nan

        return k

    # Allows the use of data collected prior to 1.3.3 release (no modification type information)
    def oldDataModificationType(self, motifInfo, TET_TREATED=False):

        modifiedPosition = int(motifInfo[self.motif])

        modifiedBaseIdentity = list(self.motif)[modifiedPosition]

        if modifiedBaseIdentity == 'A':
            modificationType = 'm6A'

        elif modifiedBaseIdentity == 'C':

            if not TET_TREATED:
                modificationType = 'm4C'
            else:
                modificationType = 'm5C'

        else:
            modificationType = 'modified_base'

        return modificationType

    def onChunk(self, referenceWindow):

        (motifDicts, refInfo, motifInfo, modifFile, undetectedOnly, oldData, start, end) = referenceWindow

        print "Start = ", start, " End = ", end

        # Read in the modifications GFF if needed:

        if not undetectedOnly:
            modReader = GffReader(modifFile)
            modifDicts = [{"seqID": x.seqid, "type": x.type, "score": x.score, "pos": x.start, "strand": x.strand, "attributes": x.attributes}
                          for x in modReader]

        # To help find the correct reference for an entry in motifs.gff
        self.selectRef = [x.FullName for x in refInfo]

        # Loop through the rows of the motifs GFF file:
        self.pre = self.ipdModel.gbmModel.post
        self.post = self.ipdModel.gbmModel.pre
        collectResults = []
        modificationsLinecount = 0

        for d in motifDicts:

            if d["type"] != '.' and undetectedOnly:
                # Go on to the next row
                continue

            ref = refInfo[self.getReferenceIndex(d, oldData)]
            self.refId = ref.ID

            k = self.fillOutEasyInformation(d, motifInfo)
            if not "motif" in k.keys():
                # If no motif is listed in this row, go on to the next row
                continue

            if d["type"] != '.' and not undetectedOnly:

                # Search sorted list of template positions in modifications GFF for a match:
                for y in modifDicts[modificationsLinecount:]:
                    if y["pos"] == d["pos"] and y["strand"] == d["strand"] and y["seqID"] == d["seqID"]:
                        break
                    modificationsLinecount += 1

                # Once the match is found, copy in the modified fraction estimate
                if modificationsLinecount <= len(modifDicts):
                    u = y["attributes"]

                    if FRAC in u:
                        k[FRAC] = float(u[FRAC])
                        k[FRAClow] = float(u[FRAClow])
                        k[FRACup] = float(u[FRACup])

            if d["type"] == '.':

                # See update to ResultsWriter: 'nMd' is 'not modified'
                k['modification'] = 'nMd'

                # Figure out modification type:
                if oldData:
                    self.modificationType = self.oldDataModificationType(motifInfo)
                else:
                    self.modificationType = motifInfo[k['motif']]

                # Select a window around the current position to use for estimation
                stop = k["tpl"] + self.post
                start = min(max(1, (k["tpl"] - self.pre)), stop)

                # Trim end coordinate to length of current template
                # end = min(end, self.ipdModel.refLength(self.refId))

                # Try to estimate the modified fraction:
                if self.modificationType == 'modified_base':
                    # In this case, we'll need the mean Ipd function:
                    self.meanIpdFunc = self.ipdModel.predictIpdFunc(self.refId)

                self.strand = k['strand']
                perSiteResults = self._summarizeReferenceRegion((start, stop))

                if self.modificationType == 'modified_base':
                    k[FRAC] = perSiteResults[self.post - 1][FRAC]
                    k[FRAClow] = perSiteResults[self.post - 1][FRAClow]
                    k[FRACup] = perSiteResults[self.post - 1][FRACup]

                else:
                    mods = self._decodePositiveControl(perSiteResults, (start, stop))
                    k[FRAC] = mods[0]
                    k[FRAClow] = mods[1]
                    k[FRACup] = mods[2]

            collectResults.append(k)

        return collectResults

    def _summarizeReferenceRegion(self, targetBounds):
        """Compute the ipd stats for a chunk of the reference"""
        (start, end) = targetBounds
        logging.info('Making summary: %d to %d' % (start, end))

        caseReferenceGroupId = self.caseCmpH5.referenceInfo(self.refId).ID
        (caseChunks, capValue) = self._fetchChunks(caseReferenceGroupId, targetBounds, self.caseCmpH5)
        self.refName = self.caseCmpH5.referenceInfo(self.refId).FullName

        # FIXME: need to pass in the correct strand
        strand = 1 if self.strand == 0 else 0
        # return [self._computePositionSyntheticControl(x, capValue, start, end) for x in caseChunks if x['strand'] == strand]

        count = 1
        res = []
        for x in caseChunks:
            if x['strand'] == strand:
                res.append(self._computePositionSyntheticControl(x, capValue, count, start, end))
                count += 1

        return res

    def _fetchChunks(self, refGroupId, targetBounds, cmpH5File):
        """Get the IPDs for each position/strand on the given reference in the given window, from the given cmpH5 file"""
        (start, end) = targetBounds

        selV = (cmpH5File.RefGroupID == refGroupId) & \
               (np.logical_not((cmpH5File.tStart > end) | (cmpH5File.tEnd < start))) & \
               (cmpH5File.MapQV > self.options.mapQvThreshold) & \
               (HitSelection.readFilter(cmpH5File, minAcc=0.82))

        hits = cmpH5File[selV]

        # FIXME -- we are dealing with the IPD format change from seconds to frames here
        # Should be handled in pbcore
        ver = cmpH5File.version[0:3]

        if ver == '1.2':
            factor = 1.0
        elif ver == '1.3' or ver == '1.4':
            # NOTE -- assuming that all movies have the same frame rate!
            fr = cmpH5File.movieInfo(1).FrameRate
            factor = 1.0 / fr
        else:
            raise Exception('Unrecognized cmp.h5 version')

        rawIpds = self._loadRawIpds(hits, start, end, factor)
        ipdVect = rawIpds['ipd']

        if ipdVect.size < 10:
            # Default is there is no coverage
            capValue = 5.0
        else:
            # Compute IPD quantiles on the current block -- will be used for trimming extreme IPDs
            capValue = np.percentile(ipdVect, self.options.cap_percentile)

        chunks = self._chunkRawIpds(rawIpds)
        return chunks, capValue

    def _loadRawIpds(self, alnHitIter, targetStart=-1, targetEnd=3e12, factor=1.0):
        """
        Get a DataFrame of the raw ipds in the give alignment hits, indexed by template position and strand.
        Factor is a normalization factor to the get units into seconds.
        """

        # Put in an empty 'starter' array -- the np.concatenate call below will fail on an empty list
        array0 = np.zeros(0, dtype=ipdRec)

        # Maintain separate lists for each strand to speed up sorting
        s0list = [array0]
        s1list = [array0]

        for aln in alnHitIter:
            # Pull out error-free position
            matched = np.logical_and(np.array([x != '-' for x in aln.read()]), np.array([x != '-' for x in aln.reference()]))

            # Normalize kinetics of the entire subread
            rawIpd = aln.IPD() * factor
            np.logical_and(np.logical_not(np.isnan(rawIpd)), matched, out=matched)

            normalization = self._subreadNormalizationFactor(rawIpd[matched])
            rawIpd /= normalization

            # Trim down to just the position that cover our interval
            referencePositions = aln.referencePositions()
            np.logical_and(referencePositions < targetEnd, matched, matched)
            np.logical_and(referencePositions >= targetStart, matched, matched)
            nm = matched.sum()

            # Bail out if we don't have any samples
            if nm == 0:
                continue

            ipd = rawIpd[matched]
            tpl = referencePositions[matched]

            dfTemp = np.zeros(nm, dtype=ipdRec)
            dfTemp['ipd'] = ipd
            dfTemp['tpl'] = tpl
            dfTemp['strand'] = aln.RCRefStrand

            if aln.RCRefStrand == 0:
                s0list.append(dfTemp)
            else:
                s1list.append(dfTemp)

        # Sort the set of ipd observations
        s0Ipds = np.concatenate(s0list)
        sortOrder = np.argsort(s0Ipds['tpl'])
        s0Ipds = s0Ipds[sortOrder]

        s1Ipds = np.concatenate(s1list)
        sortOrder = np.argsort(s1Ipds['tpl'])
        s1Ipds = s1Ipds[sortOrder]

        return np.concatenate([s0Ipds, s1Ipds])

    def _chunkRawIpds(self, rawIpds):
        """
        Return a list of view recarrays into the rawIpds recarray, one for each unique (tpl, stand) level
        """
        views = []

        # Bail out if we have no data
        if rawIpds.size == 0:
            return views

        start = 0
        tpl = rawIpds['tpl']
        strand = rawIpds['strand']

        # Start off at the first chunk
        curIdx = (tpl[0], strand[0])
        for i in xrange(1, rawIpds.shape[0]):
            newIdx = (tpl[i], strand[i])

            # In this case we are still int he same chunk -- continue
            if curIdx == newIdx:
                continue

            # In this case we have completed the chunk -- emit the chunk
            else:
                obj = {'tpl': curIdx[0], 'strand': curIdx[1], 'data': rawIpds[start:i]}
                views.append(obj)
                start = i
                curIdx = newIdx

        # Make sure to return final chunk
        obj = {'tpl': curIdx[0], 'strand': curIdx[1], 'data': rawIpds[start:]}
        views.append(obj)

        # If the user has specified a maximum coverage level to use, enforce it here -- just take the first n reads
        if self.options.maxCoverage is not None:
            maxCov = self.options.maxCoverage
            for x in views:
                d = x['data']
                d = d[0:maxCov]
                x['data'] = d

        return views

    def _subreadNormalizationFactor(self, rawIpds):
        """
        Normalize subread ipds
        """

        # Default normalization factor -- this value should very rarely get used
        if rawIpds.size < 2:
            return 0.1

        if np.isnan(rawIpds).any():
            print "got nan: %s" % str(rawIpds)

        if rawIpds.mean() < 0.0001:
            print "small"
            print "got small: %s" % str(rawIpds)

        capValue = min(10, np.percentile(rawIpds, 99))
        capIpds = np.minimum(rawIpds, capValue)
        return capIpds.mean()

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
            for suffix in self._possibleConfigs(start + 1, end):
                for c in currentChars:
                    r.append(c + suffix)

            return r

    # The following methods are used to estimate the mixture proportion when the type is known:
    def computeContextMeans(self):
        """Generate a hash of the mean ipd for all candidate contexts"""
        # print "from computeContext: ", self.likelihoodRange
        allContexts = list(set([cfg for pos in self.likelihoodRange for cfg in self.getConfigs(pos)]))
        predictions = self.ipdModel.gbmModel.getPredictions(allContexts)
        self.contextMeanTable = dict(zip(allContexts, predictions))

    # Return expected IPDs for a portion [start, end] of the sequence.
    def getContextMeans(self, start, end, sequence):
        meanVector = []
        for pos in xrange(start, end + 1):
            ctx = sequence[(pos - self.pre):(pos + self.post + 1)].tostring()
            if self.contextMeanTable.has_key(ctx):

                meanVector.append(self.contextMeanTable[ctx])
            else:
                meanVector.append(self.ipdModel.gbmModel.getPredictions([ctx]))
        return meanVector

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
    # def estimateSingleFraction(self, mu1, data, mu0, L ):
    #     a0 = self.replaceScipyNormPdf( data, mu0 )
    #     a1 = self.replaceScipyNormPdf( data, mu1 )
    #     tmp = np.multiply( 1.0/mu0 - 1.0/mu1, data )
    #     if np.exp( tmp ).sum() <= L * mu1 / mu0:
    #         return 0.0
    #     if np.exp( -tmp ).sum() <= L * mu0 / mu1:
    #         return 1.0
    #     res = fminbound(self.mixModelFn, 0.01, 0.99, args=(a0, a1), xtol=1e-02)
    #     return res
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

        print str(res)
        return res

    # End of mixture model methods for the case where the modification type is known
    # The following methods are used to estimate modified fraction for an unknown modification type:
    # Return the optimal mixing proportion in the detection case: estimate p and mu1 together.
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

    # End of unknown-type modified fraction estimation methods
    # This function replaces the ModificationDecode class used by KineticWorker:
    def replaceModificationDecode(self, sequence, pos, rawKinetics):

        self.lStart = self.pre
        self.lEnd = len(sequence) - self.post

        # Extents that we will use for likelihoods
        self.likelihoodRange = xrange(self.lStart, self.lEnd)
        self.alternateBases = dict((x, set(sequence[x])) for x in xrange(len(sequence)))

        # Compute all the required mean ipds under all possible composite hypotheses
        self.computeContextMeans()

        qvModCalls = dict()

        modSeq = a.array('c')
        modSeq.fromstring(sequence)

        # Apply the found modifications to the raw sequence
        modSeq[self.post - 1] = toMod[self.modificationType]
        modifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq)

        # Switch back to the unmodified base and re-score
        modSeq[self.post - 1] = canonicalBaseMap[toMod[self.modificationType]]
        unModifiedMeanVectors = self.getContextMeans(pos - self.post, pos + self.pre, modSeq)

        self.rawKinetics = rawKinetics
        methylFraction = self.estimateMethylatedFractions(pos, unModifiedMeanVectors, modifiedMeanVectors, ModificationPeakMask[self.modificationType])

        return methylFraction

    def _decodePositiveControl(self, kinetics, bounds):
        (kinStart, kinEnd) = bounds
        callBounds = (8, kinEnd - kinStart + 8)

        tpl = kinStart + self.post

        chunkFwd = dict((x['tpl'], x) for x in kinetics if x['strand'] == self.strand and x['coverage'] > 3)

        modCalls = []

        # Fwd sequence window
        canonicalSequence = self.ipdModel.getReferenceWindow(self.refId, 0, kinStart, kinEnd)

        # Map the raw kinetics into the frame-of reference of our sequence snippets
        def toRef(p):
            return p - (kinStart)

        def fromRef(r):
            return r + (kinStart)

        mappedChunk = dict((toRef(pos), k) for (pos, k) in chunkFwd.items())

        # Decode the modifications
        decoder = self.replaceModificationDecode(canonicalSequence, self.post - 1, mappedChunk)
        return decoder

    def _computePositionSyntheticControl(self, caseObservations, capValue, count, start, end):
        """Summarize the observed ipds at one template position/strand, using the synthetic ipd model"""

        # Compute stats on the observed ipds
        d = caseObservations['data']['ipd']
        res = dict()

        # NOTE -- this is where the strand flipping occurs -- make sure to reproduce this in the all calling methods
        strand = res['strand'] = 1 - caseObservations['strand']
        tpl = res['tpl'] = caseObservations['tpl']
        res['coverage'] = d.size

        # Store in case of methylated fraction estimtion:
        res['rawData'] = d

        # For modificationType listed as 'modified_base', use detection-mode mixture model:
        if self.modificationType == 'modified_base' and count == self.post:
            if res['coverage'] > 3:
                modelPrediction = self.meanIpdFunc(tpl, strand).item()
                x = self.detectionMixModelBootstrap(modelPrediction, d)
                res[FRAC] = x[0]
                res[FRAClow] = x[1]
                res[FRACup] = x[2]
                res['fracWts'] = x[3]
            else:
                res[FRAC] = np.nan
                res[FRACup] = np.nan
                res[FRAClow] = np.nan
                res['fracWts'] = np.nan

        return res


class KineticWorkerProcess(KineticReprocessWorker, WorkerProcess):

    """Worker that executes as a process."""

    def __init__(self, options, workQueue, resultsQueue, ipdModel):
        WorkerProcess.__init__(self, options, workQueue, resultsQueue)
        KineticReprocessWorker.__init__(self, ipdModel)


class KineticWorkerThread(KineticReprocessWorker, WorkerThread):

    """Worker that executes as a thread (for debugging purposes only)."""

    def __init__(self, options, workQueue, resultsQueue, ipdModel):
        WorkerThread.__init__(self, options, workQueue, resultsQueue)
        KineticReprocessWorker.__init__(self, ipdModel)
