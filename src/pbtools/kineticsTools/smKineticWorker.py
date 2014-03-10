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
from scipy.special import erfc
import logging

import scipy.stats as s
import numpy as np
import scipy.stats.mstats as mstats
import sys

from MixtureEstimationMethods import MixtureEstimationMethods
from MultiSiteCommon import MultiSiteCommon, canonicalBaseMap, modNames, ModificationPeakMask, FRAC, FRAClow, FRACup, log10e

from BasicLdaEnricher import BasicLdaEnricher
from PositiveControlEnricher import PositiveControlEnricher

from pbtools.kineticsTools.ModificationDecode import ModificationDecode, ModificationPeakMask

from WorkerProcess import WorkerProcess, WorkerThread


# Raw ipd record
ipdRec = [('tpl', '<u4'), ('strand', '<i8'), ('ipd', '<f4')]


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


class KineticWorker(object):

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


    def onChunk(self, referenceWindow):


        # start and end are the windows of the reference that we are responsible for reporting data from.
        # We may elect to pull data from a wider window for use with positive control

        if self.options.smBaseMod:
		(reference, smId, start, end) = referenceWindow
	else:
		(reference, start, end) = referenceWindow

	targetBounds = (start,end)
        # Trim end coordinate to length of current template
        end = min(end,self.ipdModel.refLength(reference))

        # Each chunk is from a single reference -- fire up meanIpd func on the current reference
        self.meanIpdFunc = self.ipdModel.predictIpdFunc(reference)

        # Get the cognate base at a given position
        self.cognateBaseFunc = self.ipdModel.cognateBaseFunc(reference)

        self.refId = reference

        self.sequence = self.ipdModel.getReferenceWindow(self.refId, 0, start, end)

        # Compute the data for this chunk

        if self.options.identify:
            # If we are attempting to identify modifications, get the raw data for a slightly expanded window
            # then do the decoding, then weave the modification results back into the main results

            padStart = start - 8
            padEnd = end + 8
            perSiteResults = self._summarizeReferenceRegion((padStart, padEnd), self.options.methylFraction, self.options.identify)

            if self.options.useLDA:

                # FIXME: add on a column "Ca5C" containing LDA score for each C-residue site
                # Below is an example of how to use an alternative, the BasicLdaEnricher, which does not use the positive control model
                # PositiveControlEnricher currently uses a logistic regression model trained using SMRTportal job 65203 (native E. coli)

                # lda = BasicLdaEnricher( self.ipdModel.gbmModel, self.sequence, perSiteResults, self.options.identify, self.options.modsToCall )
                lda = PositiveControlEnricher( self.ipdModel.gbmModel, self.sequence, perSiteResults )
                perSiteResults = lda.callEnricherFunction( perSiteResults )

            mods = self._decodePositiveControl(perSiteResults, (start, end))

            finalCalls = []

            # Weave together results
            for strand in [0, 1]:
                strandSign = 1 if strand == 0 else -1

                siteDict = dict((x['tpl'], x) for x in perSiteResults if start <= x['tpl'] < end and x['strand'] == strand)
                modDict = dict((x['tpl'], x) for x in mods if start <= x['tpl'] < end and x['strand'] == strand)

                # Go through the modifications - add tags for identified mods to per-site stats
                # add a 'offTarget' tag to the off target peaks.
                for (pos, mod) in modDict.items():

                    # Only convert to positive control call if we actually have enough
                    # coverage on the cognate base!
                    if siteDict.has_key(mod['tpl']):

                        # Copy mod identification data 
                        siteDict[mod['tpl']]['modificationScore'] = mod['QMod']
                        siteDict[mod['tpl']]['modification'] = mod['modification']
		
                        if self.options.methylFraction and mod.has_key(FRAC):
                            siteDict[mod['tpl']][FRAC] = mod[FRAC]
                            siteDict[mod['tpl']][FRAClow] = mod[FRAClow]
                            siteDict[mod['tpl']][FRACup] = mod[FRACup]


                    if mod.has_key('Mask'):
                        # The decoder should supply the off-target peak mask
                        mask = mod['Mask']
                        mask.append(0) # make sure we always mask the cognate position
                    else:
                        # If the decoder doesn't supply a mask - use a hard-coded version
                        # FIXME - this branch is deprecated
                        mask = ModificationPeakMask[mod['modification']]

                    # Mask out neighbor peaks that may have been caused by this mod
                    for offset in mask:
                        shadowPos = mod['tpl'] + strandSign * offset
                        if siteDict.has_key(shadowPos):
                            siteDict[shadowPos]['offTargetPeak'] = True

                finalCalls.extend(siteDict.values())

            # Sort by template position
            finalCalls.sort(key = lambda x: x['tpl'])
            return finalCalls

        else:
            if self.options.smBaseMod: 
		result = self._summarizeMolecule(smId, targetBounds, self.options.methylFraction, self.options.identify)
	    else:
		result = self._summarizeReferenceRegion(targetBounds, self.options.methylFraction, self.options.identify)		

            if self.options.useLDA and self.controlCmpH5 is None:

                # FIXME: add on a column "Ca5C" containing LDA score for each C-residue site
                # lda = BasicLdaEnricher(self.ipdModel.gbmModel, self.sequence, result, self.options.identify)
                lda = PositiveControlEnricher( self.ipdModel.gbmModel, self.sequence, result )
                results = lda.callEnricherFunction( result )

            if self.options.smBaseMod:
	 	pass
	    else:
		result.sort(key = lambda x: x['tpl'])
            
	    return result


    def _summarizeReferenceRegion(self, targetBounds, methylFractionFlag, identifyFlag):
        """Compute the ipd stats for a chunk of the reference"""
        (start, end) = targetBounds
        logging.info('Making summary: %d to %d' % (start, end))

        caseReferenceGroupId = self.caseCmpH5.referenceInfo(self.refId).ID
        (caseChunks, capValue) = self._fetchChunks(caseReferenceGroupId, targetBounds, self.caseCmpH5)
        self.refName = self.caseCmpH5.referenceInfo(self.refId).FullName
        
	if self.controlCmpH5 is None:
            # in silico control workflow -- only get data from the main 'case' cmp.h5
	    return [self._computePositionSyntheticControl(x, capValue, methylFractionFlag, identifyFlag) for x in caseChunks if x['data']['ipd'].size >= self.options.minCoverage]

        else:
            # case/control workflow -- get data from the case and control files and compare
            result = []

            contigName = self.caseCmpH5.referenceInfo(self.refId).FullName
            controlRefTable = self.controlCmpH5.referenceInfoTable

            # Make sure this RefId contains a refGroup in the control cmp.h5 file
            # if self.refId in self.controlCmpH5.referenceInfoTable.Name:
            # if self.refId in [ int( str.split('ref')[1] ) for str in self.controlCmpH5.referenceInfoTable.Name ]:
            if contigName in controlRefTable.FullName:

                controlRefRow = controlRefTable[ controlRefTable['FullName'] == contigName ][0]
                (controlChunks, _) = self._fetchChunks(controlRefRow.ID, targetBounds, self.controlCmpH5)
                controlSites = { (x['strand'], x['tpl']) : x for x in controlChunks }

                for caseChunk in caseChunks:
                    #try:
                        ## FIXME: catch None or the exception.
                        caseKey = (caseChunk['strand'], caseChunk['tpl'])
                        controlChunk = controlSites.get(caseKey) # , default = None)

                        if controlChunk and \
                            caseChunk['data']['ipd'].size > 2 and \
                            controlChunk['data']['ipd'].size > 2:
                            result.append(self._computePositionTraditionalControl(caseChunk, controlChunk, methylFractionFlag, identifyFlag))
                    #except:
                    #    pass

            return result
    
    def _summarizeMolecule(self, smId, targetBounds,methylFractionFlag, identifyFlag):
        """Compute the ipd stats for a chunk of the reference"""
        logging.info('Making summary: %d' % smId)
	
	(start,end) = targetBounds
        caseReferenceGroupId = self.caseCmpH5.referenceInfo(self.refId).ID
        (caseChunks, capValue) = self._fetchMolecule(smId,targetBounds, self.caseCmpH5)
        self.refName = self.caseCmpH5.referenceInfo(self.refId).FullName
        self.moleculeID = smId

        if self.controlCmpH5 is None:
            # in silico control workflow -- only get data from the main 'case' cmp.h5
            return [self._computePositionSyntheticControl(x, capValue, methylFractionFlag, identifyFlag) for x in caseChunks if x['data']['ipd'].size >= self.options.minCoverage]

        else:
            # case/control workflow -- get data from the case and control files and compare. may need to fix for smAnalysis. 
            result = []

            contigName = self.caseCmpH5.referenceInfo(self.refId).FullName
            controlRefTable = self.controlCmpH5.referenceInfoTable

            # Make sure this RefId contains a refGroup in the control cmp.h5 file
            # if self.refId in self.controlCmpH5.referenceInfoTable.Name:
            # if self.refId in [ int( str.split('ref')[1] ) for str in self.controlCmpH5.referenceInfoTable.Name ]:
            if contigName in controlRefTable.FullName:
                controlRefRow = controlRefTable[ controlRefTable['FullName'] == contigName ][0]
                (controlChunks, _) = self._fetchChunks(controlRefRow.ID, targetBounds, self.controlCmpH5)
                controlSites = { (x['strand'], x['tpl']) : x for x in controlChunks }

                for caseChunk in caseChunks:
                    #try:
                        ## FIXME: catch None or the exception.
                        caseKey = (caseChunk['strand'], caseChunk['tpl'])
                        controlChunk = controlSites.get(caseKey) # , default = None)

                        if controlChunk and \
                            caseChunk['data']['ipd'].size >= self.options.minCoverage and \
                            controlChunk['data']['ipd'].size >= self.options.minCoverage:
                            result.append(self._computePositionTraditionalControl(caseChunk, controlChunk, methylFractionFlag, identifyFlag))
                    #except:
                    #    pass

            return result 

    def _decodePositiveControl(self, kinetics, bounds):
        """Compute the ipd stats for a chunk of the reference"""

        (kinStart, kinEnd) = bounds
        callBounds = (8, kinEnd-kinStart + 8)

        chunkFwd = dict((x['tpl'], x) for x in kinetics if x['strand'] == 0 and x['coverage'] > self.options.identifyMinCov )
        chunkRev = dict((x['tpl'], x) for x in kinetics if x['strand'] == 1 and x['coverage'] > self.options.identifyMinCov )

        modCalls = []

        # Fwd sequence window
        canonicalSequence = self.ipdModel.getReferenceWindow(self.refId, 0, kinStart-8, kinEnd + 8)

        # Map the raw kinetics into the frame-of reference of our sequence snippets
        def toRef(p):
            return p - (kinStart - 8)
        def fromRef(r):
            return r + (kinStart - 8)

        mappedChunk = dict((toRef(pos), k) for (pos, k) in chunkFwd.items())

        # Decode the modifications
        decoder = ModificationDecode(self.ipdModel.gbmModel, canonicalSequence, mappedChunk, callBounds, self.options.methylMinCov, self.options.modsToCall, self.options.methylFraction, self.options.useLDA)

        # Map the modification positions back to normal template indices
        for (r, mod) in decoder.decode().items():
            mod["strand"] = 0
            mod['tpl'] = fromRef(r)
            modCalls.append(mod)


        # Repeat decoding on reverse sequence
        # Reverse sequence
        canonicalSequence = self.ipdModel.getReferenceWindow(self.refId, 1, kinStart-8, kinEnd + 8)

        # Map the raw kinetics into the frame-of reference of our sequence snippets
        def toRef(p):
            return len(canonicalSequence) - p + (kinStart - 8)
        def fromRef(r):
            return len(canonicalSequence) - r + (kinStart - 8)

        mappedChunk = dict((toRef(pos), k) for (pos, k) in chunkRev.items())
        decoder = ModificationDecode(self.ipdModel.gbmModel, canonicalSequence, mappedChunk, callBounds, self.options.methylMinCov, self.options.modsToCall, self.options.methylFraction, self.options.useLDA)

        for (r, mod) in decoder.decode().items():
            mod["strand"] = 1
            mod['tpl'] = fromRef(r)
            modCalls.append(mod)

        return modCalls



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
            fr = cmpH5File.movieInfoTable[0].FrameRate
            factor = 1.0/fr
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

    def _fetchMolecule(self, smId,targetBounds, cmpH5File):
        """Get the IPDs for each position/strand on the given reference in the given window, from the given cmpH5 file"""
        selV = (cmpH5File.MoleculeID == smId) 

        hits = cmpH5File[selV]
	(start,end) = targetBounds
        # FIXME -- we are dealing with the IPD format change from seconds to frames here
        # Should be handled in pbcore
        ver = cmpH5File.version[0:3]

        if ver == '1.2':
                factor = 1.0
        elif ver == '1.3' or ver == '1.4':
            # NOTE -- assuming that all movies have the same frame rate!
                fr = cmpH5File.movieInfoTable[0].FrameRate
                factor = 1.0/fr
        else:
                raise Exception('Unrecognized cmp.h5 version')

        rawIpds = self._loadRawIpds(hits, start,end,factor)
	ipdVect = rawIpds['ipd']
        if ipdVect.size < 10:
            # Default is there is no coverage
                capValue = 10
		rawIpds['ipd'] = np.minimum(rawIpds['ipd'],capValue)
        else:
           # Compute IPD quantiles on the current block -- will be used for trimming extreme IPDs
                capValue = min(10,np.percentile(ipdVect, self.options.cap_percentile))
		rawIpds['ipd'] = np.minimum(rawIpds['ipd'],capValue)

        chunks = self._chunkRawIpds(rawIpds)
        return chunks, capValue

    def _loadRawIpds(self, alnHitIter, targetStart = -1, targetEnd = 3e12, factor = 1.0):
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
            matched = np.logical_and( np.array([x != '-' for x in aln.read()]), np.array([x != '-' for x in aln.reference()]))
	    if (self.options.errWin):
		mIdx = []
		mIdx.extend(np.array([i-1 for i in range(0,len(matched)) if matched[i]==False]))
		mIdx.extend(np.array([i+1 for i in range(0,len(matched)) if matched[i]==False]))
		mIdx = np.maximum(mIdx,0)
		mIdx = np.minimum(mIdx,len(matched))
		mIdx = list(set(mIdx))
		matched[mIdx] = False

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
        for i in xrange(1,rawIpds.shape[0]):
            newIdx = (tpl[i], strand[i])

            # In this case we are still int he same chunk -- continue
            if curIdx == newIdx:
                continue

            # In this case we have completed the chunk -- emit the chunk
            else:
                obj = { 'tpl': curIdx[0], 'strand': curIdx[1], 'data': rawIpds[start:i] }
                views.append(obj)
                start = i
                curIdx = newIdx

        # Make sure to return final chunk
        obj = { 'tpl': curIdx[0], 'strand': curIdx[1], 'data': rawIpds[start:] }
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

        if self.options.smBaseMod:
		capValue = min(10, np.percentile(rawIpds, 99))
		capIpds = rawIpds[-(rawIpds>capValue)]
	else:
		capValue = min(10, np.percentile(rawIpds, 99))
        	capIpds = np.minimum(rawIpds, capValue)
        
	return capIpds.mean()

    def computeObservationPValue(self, siteObs):
        """
        Compute a p-value on the observation of a kinetic event
        """

        # p-value of detection -- FIXME needs much more thought here!
        # p-value computation (slightly robustified Gaussian model)
        #  emf - rms fractional error of background model
        #  em - rms error of background model = um * emf
        #  um - predicted mean of unmodified ipd from model
        #  uo - (trimmed) observed mean ipd
        #  eo - (trimmed) standard error of observed mean (std / sqrt(coverage))
        #  Null model is ~N(um, em^2 + eo^2)
        #  Then compute standard gaussian p-value = erfc((uo-um) / sqrt(2 * (em^2 + eo^2))) / 2
        # FIXME? -- right now we only detect the case where the ipd gets longer.

        um = siteObs['modelPrediction']

        # FIXME -- pipe through model error
        em = 0.1 * um
        # em = model.fractionalModelError * em

        uo = siteObs['tMean']
        eo = siteObs['tErr']

        pvalue = erfc((uo - um) / sqrt(2 * (em**2 + eo**2))) / 2
        return pvalue.item()

    
    def computeObservationTstatistic(self, siteObs):
        """
        Compute a p-value on the observation of a kinetic event
        """

        # p-value of detection -- FIXME needs much more thought here!
        # p-value computation (slightly robustified Gaussian model)
        #  emf - rms fractional error of background model
        #  em - rms error of background model = um * emf
        #  um - predicted mean of unmodified ipd from model
        #  uo - (trimmed) observed mean ipd
        #  eo - (trimmed) standard error of observed mean (std / sqrt(coverage))
        #  Null model is ~N(um, em^2 + eo^2)
        #  Then compute standard gaussian p-value = erfc((uo-um) / sqrt(2 * (em^2 + eo^2))) / 2
        # FIXME? -- right now we only detect the case where the ipd gets longer.

        um = siteObs['modelPrediction']

        # FIXME -- pipe through model error
        #em = 0.06 * um + 0.12 * um**2.0
        em = 0.01 + 0.03*um + 0.06*um**(1.7)
        # em = model.fractionalModelError * em

        uo = siteObs['tMean']
        eo = siteObs['tErr']

        import scipy.stats as s

        t = -(uo - um) / sqrt(em**2 + eo**2)
        return t


    def computeObservationPValueTTest(self, siteObs):
        t = siteObs['tStatistic']
        df = max(1, siteObs['coverage'] - 1)

        pvalue = s.t.cdf(t, df)
        return pvalue.item()




    def _computePositionSyntheticControl(self, caseObservations, capValue, methylFractionFlag, identifyFlag):
        """Summarize the observed ipds at one template position/strand, using the synthetic ipd model"""

        # Compute stats on the observed ipds
        d = caseObservations['data']['ipd']
        res = dict()

        # ref00000x name
        res['refId'] = self.refId
	
        # FASTA header name
        res['refName'] = self.refName

	# single molecule id if in smBaseMod model
	if self.options.smBaseMod:
		res['moleculeID'] = self.moleculeID

        # NOTE -- this is where the strand flipping occurs -- make sure to reproduce this in the all calling methods
        strand = res['strand'] = 1 - caseObservations['strand']
        tpl = res['tpl'] = caseObservations['tpl']
        res['coverage'] = d.size


        # Don't compute these stats - they just take time and confuse things
        # res['mean'] = d.mean().item()
        # res['median'] = np.median(d).item()
        # res['std'] = np.std(d).item()

        # Compute the predicted IPD from the model
        # # NOTE! The ipd model is in the observed read strand
        modelPrediction = self.meanIpdFunc(tpl, strand).item()
        res['modelPrediction'] = modelPrediction

        res['base'] = self.cognateBaseFunc(tpl, strand)

	# Store in case of methylated fraction estimtion:
	res['rawData'] = d

        # Try a hybrid capping approach -- cap at the higher of
        #  - 5x the model prediction
        #  - 90th percentile of the local data (at low coverage we pick a lower percentile to ensure we trim the highest datapoint
        #  - global cap value

        percentile = min(90, (1.0 - 1.0/(d.size - 1))*100)
        localPercentile = np.percentile(d, percentile)
        capValue = max(capValue, 4.0 * modelPrediction, localPercentile)

        # np.minimum(d, capValue, out=d)  # this version will send capped IPDs to modified fraction estimator
        d = np.minimum(d, capValue)

        # Trimmed stats
        res['tMean'] = d.mean().item()
        res['tErr'] = np.std(d).item() / sqrt(d.size)

	res['ipdRatio'] = res['tMean'] / res['modelPrediction']

        # Don't know the modification yet
        res["modification"] = "."

        # use ttest-based pvalue
        # res['pvalue'] = self.computeObservationPValue(res)
        res['tStatistic'] = self.computeObservationTstatistic(res)
        res['pvalue'] = self.computeObservationPValueTTest(res)


        pvalue = max(sys.float_info.min, res['pvalue'])
        score = round(-10.0 * math.log10(pvalue))
        res['score'] = score

        # If the methylFractionFlag is set, then estimate fraction using just modelPrediction in the detection case.
        if methylFractionFlag and pvalue < self.options.pvalue and not identifyFlag:
            if res['coverage'] > self.options.methylMinCov:
                modelPrediction = self.meanIpdFunc(tpl, strand).item()

                # Instantiate mixture estimation methods:
                mixture = MixtureEstimationMethods(self.ipdModel.gbmModel.post, self.ipdModel.gbmModel.pre, res, self.options.methylMinCov)
                x = mixture.detectionMixModelBootstrap(modelPrediction, d)
                # x = self.detectionMixModelBootstrap(modelPrediction, d)

                res[FRAC] = x[0]
                res[FRAClow] = x[1]
                res[FRACup] = x[2]
            else:
                res[FRAC] = np.nan
                res[FRACup] = np.nan
                res[FRAClow] = np.nan

        # print res
        return res


##
## straight port from R's p.adjust. 
##
    def _cummin(v):
        r = array([v[0]]*len(v))
        r[0] = v[0]
        for i in xrange(1, len(v)):
            r[i] = r[i-1] if r[i-1] < v[i] else v[i]
        return r
        
    def _BH_FDR(pvals):
        s = array(range(len(pvals), 0, -1))
        o = array(argsort(pvals)[::-1])
        r = array(argsort(o))
        return [ 1 if x > 1 else x for x in (_cummin(float(len(pvals))/s * pvals[o]))[r] ]
    
##
## Null simulation. the test below assumes that IPDs are normal after
## capping and logging. FIXME: permutation based
##
# def sim(N=100):
#     return [ _tTest(np.exp(x1),np.exp(x2), 100)['pvalue'] for x1,x2 in
#              zip([ np.random.normal(size=100) for g in range(0, N)],
#                  [ np.random.normal(size=100) for g in range(0, N) ] )]
#
##  _tTest(np.exp(np.random.normal(1.5, size = 100)), np.exp(np.random.normal(1., size = 100)))
##
    def _tTest(x, y, exclude = 95):
        """Compute a one-sided Welsh t-statistic."""
        def cappedSlog(v):
            q = np.percentile(v, exclude)
            v2 = v.copy()
            v2 = v2[~np.isnan(v2)]
            v2[v2 <= 0] = 1./(75 + 1)
            v2[v2 > q] = q
            return np.log(v2)
        x1 = cappedSlog(x)
        x2 = cappedSlog(y)
        sx1 = np.var(x1)/len(x1)
        sx2 = np.var(x2)/len(x2)
        stat = (np.mean(x1)-np.mean(x2))/np.sqrt(sx1 + sx2)
        #df   = (sx1 + sx2)**2 / (sx1**2/(len(x1)-1) + sx2**2/(len(x2) - 1))
        #pval = 1 - scidist.t.cdf(stat, df)

        # Scipy's t distribution CDF implementaton has inadequate
        # precision.  We have switched to the normal distribution for
        # better behaved p values.
        pval = 0.5*erfc(stat/sqrt(2))

        return {'testStatistic':stat, 'pvalue':pval}

    def _computePositionTraditionalControl(self, caseObservations, controlObservations, methylFractionFlag, identifyFlag, testProcedure = _tTest):
        """Summarize the observed ipds at one template position/strand, using a case-control analysis"""
        # Compute stats on the observed ipds
        caseData = caseObservations['data']['ipd']
        controlData = controlObservations['data']['ipd']

        res = dict()
        res['refId'] = self.refId

        # FASTA header name
        res['refName'] = self.refName
	
	# single molecule id if in smBaseMod model
        if self.options.smBaseMod:
                res['moleculeID'] = self.moleculeID

        strand = res['strand'] = 1 - caseObservations['strand']
        tpl = res['tpl'] = caseObservations['tpl']
        res['base'] = self.cognateBaseFunc(tpl, strand)

        res['coverage'] = int(round((caseData.size + controlData.size)/2.0)) ## need a coverage annotation

        res['caseCoverage'] = caseData.size
        res['controlCoverage'] = controlData.size

        res['caseMean'] = caseData.mean().item()
        res['caseMedian'] = np.median(caseData).item()
        res['caseStd'] = np.std(caseData).item()

        res['controlMean'] = controlData.mean().item()
        res['controlMedian'] = np.median(controlData).item()
        res['controlStd'] = np.std(controlData).item()

        trim = (0.001, 0.03)
        ctrlMean = mstats.trimmed_mean(controlData, trim).item()
        if abs(ctrlMean) > 1e-3:
            res['ipdRatio'] = ( mstats.trimmed_mean(caseData, trim).item() / ctrlMean )
        else:
            res['ipdRatio'] = 1.0

        testResults = testProcedure(caseData, controlData)
        res['testStatistic'] = testResults['testStatistic']
        res['pvalue'] = testResults['pvalue']

        pvalue = max(sys.float_info.min, res['pvalue'])
        res['score'] = round(-10.0 * math.log10(pvalue))

        # If the methylFractionFlag is set, then estimate fraction using just modelPrediction in the detection case.
        if methylFractionFlag and pvalue < self.options.pvalue and not identifyFlag:
            if res['controlCoverage'] > self.options.methylMinCov and res['caseCoverage'] > self.options.methylMinCov:

                # Instantiate mixture estimation methods:
                mixture = MixtureEstimationMethods(self.ipdModel.gbmModel.post, self.ipdModel.gbmModel.pre, res, self.options.methylMinCov)
                x = mixture.detectionMixModelBootstrap( res['controlMean'], caseData )

                res[FRAC] = x[0]
                res[FRAClow] = x[1]
                res[FRACup] = x[2]
            else:
                res[FRAC] = np.nan
                res[FRACup] = np.nan
                res[FRAClow] = np.nan

        return res



class KineticWorkerProcess(KineticWorker, WorkerProcess):
    """Worker that executes as a process."""
    def __init__(self, options, workQueue, resultsQueue, ipdModel):
        WorkerProcess.__init__(self, options, workQueue, resultsQueue)
        KineticWorker.__init__(self, ipdModel)


class KineticWorkerThread(KineticWorker, WorkerThread):
    """Worker that executes as a thread (for debugging purposes only)."""
    def __init__(self, options, workQueue, resultsQueue, ipdModel):
        WorkerThread.__init__(self, options, workQueue, resultsQueue)
        KineticWorker.__init__(self, ipdModel)




