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

import logging
import os
import re
import h5py
import numpy as np
import ctypes as C
from pbtools.kineticsTools.sharedArray import SharedArray

byte = np.dtype('byte')
float32 = np.dtype('float32')
uint8 = np.dtype('uint8')

# Map for ascii encoded bases to integers 0-3 -- will be used to define a 24-bit lookup code
# for fetching predicted IPDs from the kinetic LUT.

# We start everything at 0, so anything will map to 'A' unless it appears in this table
lutCodeMap = np.zeros(256, dtype=uint8)
maps = { 'a': 0, 'A':0, 'c':1, 'C':1, 'g':2, 'G':2, 't':3, 'T':3 }
for k in maps:
    lutCodeMap[ord(k)] = maps[k]
lutReverseMap = { 0: 'A', 1:'C', 2:'G', 3:'T'}

seqCodeMap = np.ones(256, dtype=uint8) * 4
for k in maps:
    seqCodeMap[ord(k)] = maps[k]
seqMap = { 0: 'A', 1:'C', 2:'G', 3:'T', 4:'N'}
seqMapComplement = { 0: 'T', 1:'G', 2:'C', 3:'A', 4:'N'}

# Base letters for modification calling
# 'H' : m6A, 'I' : m5C, 'J' : m4C, 'K' : m5C/TET
baseToCode = { 'N': 0, 'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'H' : 4, 'I':5, 'J':6, 'K' :7 }
baseToCanonicalCode =  { 'N':0, 'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3, 'H' : 0, 'I':1, 'J':1, 'K' :1 }

codeToBase = dict([(y,x) for (x,y) in baseToCode.items()])



class GbmContextModel(object):
    """
    Class for computing ipd predictions on contexts. Evaluate the GBM tree model for a list of contexts
    Contexts may contain arbitrary combinations of modified bases
    """

    def __init__(self, modelH5Group, modelIterations=-1):

        # This will hold the ctypes function pointer
        # It will be lazily initialized
        self.nativeInnerPredict = None
        self.nativeInnerPredictCtx = None

        def ds(name):
            return modelH5Group[name][:]

        self.varNames = ds("VarNames")
        self.modFeatureIdx = dict( (int(self.varNames[x][1:]), x) for x in range(len(self.varNames)) if self.varNames[x][0] == 'M' )
        self.canonicalFeatureIdx = dict( (int(self.varNames[x][1:]), x) for x in range(len(self.varNames)) if self.varNames[x][0] == 'R' )

        self.pre = 10
        self.post = 4
        self.ctxSize = self.pre + self.post + 1

        self.splitVar = ds("Variables")
        self.leftNodes = ds("LeftNodes")
        self.rightNodes = ds("RightNodes")
        self.missingNodes = ds("MissingNodes")

        self.splitVar16 = self.splitVar.astype(np.int16)


        self.splitCodes = ds("SplitCodes").astype(np.float32)

        self.cSplits = ds("CSplits")
        self.maxCSplits = self.cSplits.shape[1]

        self.initialValue = ds("InitialValue").astype(np.float32)[0]

        exp = 2**np.arange(self.cSplits.shape[1] - 1, -1, -1)
        self.bSplits = ((self.cSplits > 0) * exp).sum(1)

        # total number of trees in model
        self.nTrees = self.splitVar.shape[0]
        self.treeSize = self.splitVar.shape[1]


        offsets = np.floor(np.arange(0, self.leftNodes.size)/self.treeSize)*self.treeSize
        offsets = offsets.astype(np.int32)

        self.leftNodesOffset = self.leftNodes.flatten().astype(np.int32) + offsets
        self.rightNodesOffset = self.rightNodes.flatten().astype(np.int32) + offsets
        self.missingNodesOffset = self.missingNodes.flatten().astype(np.int32) + offsets


        self.splitCodesCtx = self.splitCodes.copy().flatten()

        splitCodesCtxView = self.splitCodesCtx.view()
        splitCodesCtxView.dtype = np.uint32

        # Pack the cSplits as a bit array directly into the splitCode array
        # using an uin32 view of the splitCode array
        flatSplitVar = self.splitVar.flatten()

        powOfTwo = 2**np.arange(self.maxCSplits)

        for i in xrange(self.splitCodesCtx.shape[0]):

            if flatSplitVar[i] != -1:
                # This is a pointer to a cSplit row -- pack the csplit into a unit32, then overwirte
                # this slot of the ctxSplitCodes
               cs = self.cSplits[int(self.splitCodesCtx[i]),:]
               v = (powOfTwo * (cs > 0)).sum()

               splitCodesCtxView[i] = v

        # If the user has requested fewer iterations, update nTrees
        if modelIterations > 0:
            self.nTrees = modelIterations



    def _initNativeTreePredict(self):
        """
        Initialization routine the C tree-predict method
        Needs to be invoked lazily because the native function pointer cannot be pickled
        """

        import platform

        if platform.system() == "Windows":

            libfn = "tree_predict.dll"
            path = os.path.dirname(os.path.abspath(__file__))
            windowsLib = path + os.path.sep + libfn

            if os.path.exists(windowsLib):
                self._lib = np.ctypeslib.load_library(libfn, path)
            else:
                raise Exception("can't find tree_predict.dll")
        else:
            curPath = os.path.dirname(os.path.abspath(__file__))
            srcPath = os.path.dirname(os.path.dirname(curPath))

            DLL_PATH1 = curPath + os.path.sep + "tree_predict.so"
            DLL_PATH2 = srcPath + os.path.sep + "tree_predict.so"

            if os.path.exists(DLL_PATH1):
                self._lib = np.ctypeslib.load_library("tree_predict.so", curPath)
            elif os.path.exists(DLL_PATH2):
                self._lib = np.ctypeslib.load_library("tree_predict.so", srcPath)
            else:
                raise Exception("can't find tree_predict.so")

        lpb = self._lib

        lpb.init_native.argtypes = [C.c_int]

        fp = C.POINTER(C.c_float)
        fpp = C.POINTER(fp)
        ip = C.POINTER(C.c_int)
        sp = C.POINTER(C.c_int16)
        ui64p = C.POINTER(C.c_uint64)

        args = [fp, fpp, C.c_int, ip, ip, ip, fp, ip, ip, ip, C.c_float, C.c_int, C.c_int, C.c_int ]
        lpb.innerPredict.argtypes = args
        self.nativeInnerPredict = lpb.innerPredict

        # Fast version

        #void innerPredictCtx(
        #    int ctxSize, float radPredF[], uint64_t contextPack[], int cRows, 
        #    int16 left[], int16 right[], int16 missing[], float splitCode[], int16 splitVar[], 
        #    int varTypes[], float initialValue, int treeSize, int numTrees, int maxCSplitSize)

        args = [C.c_int, fp, ui64p, C.c_int, ip, ip, ip, fp, sp, ip, C.c_float, C.c_int, C.c_int, C.c_int ]
        lpb.innerPredictCtx.argtypes = args
        self.nativeInnerPredictCtx = lpb.innerPredictCtx



    def getPredictionsSlow(self, ctxStrings, nTrees = None):
        """Compute IPD predictions for arbitrary methylation-containing contexts."""
        # C prototype that we call:
        # void innerPredict(
        #   float[] radPredF,
        #   IntPtr[] dataMatrix,
        #   int cRows, int[] left, int[] right, int[] missing,
        #   float[] splitCode, int[] splitVar, int[] cSplits,
        #   int[] varTypes, float initialValue,
        #   int treeSize, int numTrees, int maxCSplitSize);

        # Make sure native library is initialized
        if self.nativeInnerPredict is None:
            self._initNativeTreePredict()


        def fp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_float))

        def ip(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int))

        if nTrees is None:
            nTrees = self.nTrees

        n = len(ctxStrings)

        mCols = [ np.zeros(n, dtype=np.float32) for x in xrange(self.ctxSize) ]
        rCols = [ np.zeros(n, dtype=np.float32) for x in xrange(self.ctxSize) ]

        for stringIdx in xrange(len(ctxStrings)):
            s = ctxStrings[stringIdx]

            for i in xrange(len(s)):
                mCols[i][stringIdx] = baseToCode[s[i]]
                rCols[i][stringIdx] = baseToCanonicalCode[s[i]]


        dataPtrs = (C.POINTER(C.c_float)*(2*self.ctxSize))()

        varTypes = np.zeros(2*self.ctxSize, dtype=np.int32)

        for i in xrange(self.ctxSize):
            dataPtrs[self.modFeatureIdx[i]] = mCols[i].ctypes.data_as(C.POINTER(C.c_float))
            dataPtrs[self.canonicalFeatureIdx[i]] = rCols[i].ctypes.data_as(C.POINTER(C.c_float))

            varTypes[self.modFeatureIdx[i]] = 8
            varTypes[self.canonicalFeatureIdx[i]] = 4


        self.predictions = np.zeros(len(ctxStrings), dtype=np.float32)

        self.nativeInnerPredict(
            fp(self.predictions), dataPtrs,
            n, ip(self.leftNodes), ip(self.rightNodes), ip(self.missingNodes),
            fp(self.splitCodes), ip(self.splitVar), ip(self.cSplits),
            ip(varTypes), self.initialValue, self.treeSize, nTrees, self.maxCSplits)

        return np.exp(self.predictions)


    def getPredictions(self, ctxStrings, nTrees = None):
        """Compute IPD predictions for arbitrary methylation-containing contexts."""
        # C prototype that we call:
        # void innerPredictCtx(
        #   int ctxSize, float[] radPredF,
        #   int[] contextPack,
        #   int cRows, int[] left, int[] right, int[] missing,
        #   float[] splitCode, int[] splitVar,
        #   int[] varTypes, float initialValue,
        #   int treeSize, int numTrees, int maxCSplitSize);

        # Make sure native library is initialized
        if self.nativeInnerPredictCtx is None:
            self._initNativeTreePredict()


        def fp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_float))

        def ip(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int))

        def ulp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_uint64))

        def sp(arr):
            return arr.ctypes.data_as(C.POINTER(C.c_int16))


        n = len(ctxStrings)

        if nTrees is None:
            nTrees = self.nTrees

        packCol = np.zeros(n, dtype=np.uint64)

        for stringIdx in xrange(len(ctxStrings)):
            s = ctxStrings[stringIdx]
            code = 0

            for i in xrange(len(s)):
                modBits = baseToCode[s[i]]

                slotForPosition = self.modFeatureIdx[i]

                code = code | (modBits << (4*slotForPosition))

            packCol[stringIdx] = code

        # print packed base codes
        #for v in packCol.flatten():
        #    print v
        #    for i in np.arange(12):
        #        print "%d: %o" % (i,  (v.item() >> (5*i)) & 0x1f)

        varTypes = np.zeros(2*self.ctxSize, dtype=np.int32)

        for i in xrange(self.ctxSize):
            varTypes[self.modFeatureIdx[i]] = 8
            varTypes[self.canonicalFeatureIdx[i]] = 4


        self.predictions = np.zeros(len(ctxStrings), dtype=np.float32)

        self.nativeInnerPredictCtx(
            self.ctxSize, fp(self.predictions), ulp(packCol),
            n, ip(self.leftNodesOffset), ip(self.rightNodesOffset), ip(self.missingNodesOffset),
            fp(self.splitCodesCtx), sp(self.splitVar16),
            ip(varTypes), self.initialValue, self.treeSize, nTrees, self.maxCSplits)

        return np.exp(self.predictions)


class IpdModel:
    """
    Predicts the IPD of an any context, possibly containing multiple modifications.
    We use a 4^12 entry LUT to get the predictions for contexts without modifications,
    then we use the GbmModel to get predictions in the presence of arbitrary mods.
	Note on the coding scheme.  For each contig we store a byte-array that has size = contig.length + 2*self.pad
	The upper 4 bits contain a lookup into seqReverseMap, which can contains N's. This is used for giving
	template snippets that may contains N's if the reference sequence does, or if the snippet
	The lowe 4 bits contain a lookup into lutReverseMap, which
    """

    def __init__(self, fastaRecords, modelFile=None, modelIterations=-1):
        """
        Load the reference sequences and the ipd lut into shared arrays that can be
        used as numpy arrays in worker processes.
        fastaRecords is a list of FastaRecords, in the cmp.h5 file order
        """

        self.pre = 10
        self.post = 4

        self.pad = 30
        self.base4 = 4 ** np.array(range(self.pre + self.post + 1))


        self.refDict = {}
        self.refLengthDict = {}

        for contig in fastaRecords:
            if contig.id is None:
                # This contig has no mapped reads -- skip it
                continue

            rawSeq = contig.sequence
            refSeq = np.fromstring(rawSeq, dtype=byte)

            # Store the reference length
            self.refLengthDict[contig.id] = len(rawSeq)

            # Make a shared array
            sa = SharedArray(dtype='B', shape=len(rawSeq) + self.pad * 2)
            saWrap = sa.getNumpyWrapper()

            # Lut Codes convert Ns to As so that we don't put Ns into the Gbm Model
            # Seq Codes leaves Ns as Ns for getting reference snippets out
            innerLutCodes = lutCodeMap[refSeq]
            innerSeqCodes = seqCodeMap[refSeq]
            innerCodes = np.bitwise_or(innerLutCodes, np.left_shift(innerSeqCodes, 4))

            saWrap[self.pad:(len(rawSeq) + self.pad)] = innerCodes

            # Padding codes -- the lut array is padded with 0s the sequence array is padded with N's (4)
            outerCodes = np.left_shift(np.ones(self.pad, dtype=uint8) * 4, 4)
            saWrap[0:self.pad] = outerCodes
            saWrap[(len(rawSeq) + self.pad):(len(rawSeq) + 2*self.pad)] = outerCodes

            self.refDict[contig.id] = sa

        # No correction factor for IPDs everything is normalized to 1
        self.meanIpd = 1

        # Find and open the ipd model file

        if modelFile:
            self.lutPath = modelFile
        else:
            self.lutPath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + "kineticLut.h5"

        if os.path.exists(self.lutPath):
            h5File = h5py.File(self.lutPath, mode='r')

            gbmModelGroup = h5File["/AllMods_GbmModel"]
            self.gbmModel = GbmContextModel(gbmModelGroup, modelIterations)

            # We will use the LUT for the null model if it's available
            if "NullModel_Lut" in h5File.keys():
                nullModelGroup = h5File["/NullModel_Lut"]
                self._loadIpdTable(nullModelGroup)
                self.haveLut = True
                self.predictIpdFunc = self.predictIpdFuncLut
            else:
                self.haveLut = False
                self.predictIpdFunc = self.predictIpdFuncModel
        else:
            logging.info("Couldn't find model file: %s" % self.lutPath)


    def _loadIpdTable(self, nullModelGroup):
        """
        Read the null kinetic model into a shared numpy array dataset
        """
        nullModelDataset = nullModelGroup["KineticValues"]

        # assert that the dataset is a uint8
        assert(nullModelDataset.dtype == uint8)

        # Construct a 'shared array' (a numpy wrapper around some shared memory
        # Read the LUT into this table
        self.sharedArray = SharedArray('B', nullModelDataset.shape[0])
        lutArray = self.sharedArray.getNumpyWrapper()
        nullModelDataset.read_direct(lutArray)

        # Load the second-level LUT
        self.floatLut = nullModelGroup["Lut"][:]

    def refLength(self, refId):
        return self.refLengthDict[refId]

    def cognateBaseFunc(self, refId):
        """
        Return a function that returns a snippet of the reference sequence around a given position
        """

        # FIXME -- what is the correct strand to return?!
        # FIXME -- what to do about padding when the snippet runs off the end of the reference
        # how do we account for / indicate what is happening
        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, tplStrand):

            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand==0:
                slc = refArray[tplPos]
                slc = np.right_shift(slc, 4)
                return seqMap[slc]

            # Reverse strand
            else:
                slc = refArray[tplPos]
                slc = np.right_shift(slc, 4)
                return seqMapComplement[slc]

        return f

    def snippetFunc(self, refId, pre, post):
        """
        Return a function that returns a snippet of the reference sequence around a given position
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, tplStrand):
            """Closure for returning a reference snippet. The reference is padded with N's for bases falling outside the extents of the reference"""
            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand==0:
                slc = refArray[(tplPos-pre):(tplPos+1+post)]
                slc = np.right_shift(slc, 4)
                return "".join(seqMap[x] for x in slc)

            # Reverse strand
            else:
                slc = refArray[(tplPos+pre):(tplPos-post-1):-1]
                slc = np.right_shift(slc, 4)
                return "".join(seqMapComplement[x] for x in slc)

        return f


    def getReferenceWindow(self, refId, tplStrand, start, end):
        """
        Return  a snippet of the reference sequence
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        # adjust position for reference padding
        start += self.pad
        end += self.pad

        # Forward strand
        if tplStrand==0:
            slc = refArray[start:end]
            slc = np.right_shift(slc, 4)
            return "".join(seqMap[x] for x in slc)

        # Reverse strand
        else:
            slc = refArray[end:start:-1]
            slc = np.right_shift(slc, 4)
            return "".join(seqMapComplement[x] for x in slc)



    def predictIpdFuncLut(self, refId):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        # Materialized the numpy wrapper around the shared data
        refArray = self.refDict[refId].getNumpyWrapper()
        lutArray = self.sharedArray.getNumpyWrapper()
        floatLut = self.floatLut

        def f(tplPos, tplStrand):

            # skip over the padding
            tplPos += self.pad

            # Forward strand
            if tplStrand==0:
                slc = np.bitwise_and(refArray[(tplPos+self.pre):(tplPos-self.post-1):-1], 0xf)

            # Reverse strand
            else:
                slc = 3 - np.bitwise_and(refArray[(tplPos-self.pre):(tplPos+1+self.post)], 0xf)

            code = (self.base4 * slc).sum()
            return floatLut[max(1, lutArray[code])]

        return f

    def predictIpdFuncModel(self, refId):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        # Materialized the numpy wrapper around the shared data
        snipFunction = self.snippetFunc(refId, self.post, self.pre)

        def f(tplPos, tplStrand):
            # Get context string
            context = snipFunction(tplPos, tplStrand)

            # Get prediction
            return self.gbmModel.getPredictions([context])[0]

        return f

    def modPredictIpdFunc(self, refId, mod):
        """
        Each (pre+post+1) base context gets mapped to an integer
        by converting each nucleotide to a base-4 number A=0, C=1, etc,
        and treating the 'pre' end of the context of the least significant
        digit.  This code is used to lookup the expected IPD in a
        pre-computed table.  Contexts near the ends of the reference
        are coded by padding the context with 0
        """

        refArray = self.refDict[refId].getNumpyWrapper()

        def f(tplPos, relativeModPos, readStrand):

            # skip over the padding
            tplPos += self.pad

            # Read sequence matches forward strand
            if readStrand==0:
                slc = 3 - np.bitwise_and(refArray[(tplPos-self.pre):(tplPos+1+self.post)], 0xf)

            # Reverse strand
            else:
                slc = np.bitwise_and(refArray[(tplPos+self.pre):(tplPos-self.post-1):-1], 0xf)

            # Modify the indicated position
            slc[relativeModPos + self.pre] = baseToCode[mod]

            slcString = "".join([codeToBase[x] for x in slc])

            # Get the prediction for this context
            #return self.gbmModel.getPredictions([slcString])[0]
            return 0.0

        return f


