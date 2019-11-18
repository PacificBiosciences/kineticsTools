from __future__ import print_function
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

from collections import namedtuple, defaultdict
import cProfile
import logging
import os.path
import os
from multiprocessing import Process
import pickle
import sys
import math
import re

import numpy as np

from pbcore.io import GffWriter, Gff3Record

from kineticsTools.pipelineTools import consumer

DEFAULT_NCHUNKS = 256
# Labels for modified fraction:
FRAC = 'frac'
FRAClow = 'fracLow'
FRACup = 'fracUp'

log = logging.getLogger(__name__)


class ResultCollectorProcess(Process):

    """
    Gathers results and writes to a file.
    """

    def __init__(self, options, resultsQueue):
        Process.__init__(self)
        self.daemon = True
        self.options = options
        self._resultsQueue = resultsQueue

    def _run(self):
        log.info("Process %s (PID=%d) started running" % (self.name, self.pid))

        self.onStart()

        nextChunkId = 0
        chunkCache = {}

        sentinelsReceived = 0
        while sentinelsReceived < self.options.numWorkers:
            result = self._resultsQueue.get()
            self._resultsQueue.task_done()

            if result is None:
                sentinelsReceived += 1
            else:
                # Write out chunks in chunkId order.
                # Buffer received chunks until they can be written in order
                (chunkId, datum) = result
                chunkCache[chunkId] = datum

                # The rawData field is large and unused. Delete it to mitigate
                # risk of OOM problems
                for column in datum:
                    if 'rawData' in column:
                        del column['rawData']

                # Write out all the chunks that we can
                while nextChunkId in chunkCache:
                    nextChunk = chunkCache.pop(nextChunkId)
                    self.onResult(nextChunk)

                    nextChunkId += 1

        log.info("Result thread shutting down...")
        self.onFinish()

    def run(self):

        if self.options.doProfiling:
            cProfile.runctx("self._run()",
                            globals=globals(),
                            locals=locals(),
                            filename="profile-%s.out" % self.name)
        else:
            self._run()

    # ==================================
    # Overridable interface begins here.
    #
    def onStart(self):
        pass

    def onResult(self, result):
        pass

    def onFinish(self):
        pass


class KineticsWriter(ResultCollectorProcess):

    def __init__(self, options, resultQueue, refInfo, ipdModel):
        ResultCollectorProcess.__init__(self, options, resultQueue)

        self.refInfo = refInfo
        self.ipdModel = ipdModel

    @consumer
    def msCsvConsumer(self, filename):
        """
        Consume IPD summary rows and write them to csv
        """

        # Open the csv file
        f = self.openWriteHandle(filename)
        delim = ","

        cols = ["refName", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction", "ipdRatio", "coverage", "signal", "variance", "MSscore"]

        # Special cases for formatting columns of the csv
        handlers = dict()
        threeF = lambda x: "%.3f" % x

        handlers["refName"] = lambda x: "\"%s\"" % x

        handlers["tpl"] = lambda x: str(x.item() + 1)
        handlers["score"] = lambda x: "%d" % x

        handlers["tMean"] = threeF
        handlers["modelPrediction"] = threeF
        handlers["caseMean"] = threeF
        handlers["controlMean"] = threeF
        handlers["ipdRatio"] = threeF
        handlers["pvalue"] = lambda x: "%.3e" % x

        handlers["controlStd"] = threeF
        handlers["controlStd"] = threeF
        handlers["tErr"] = threeF

        fourF = lambda x: "%.4f" % x
        handlers["signal"] = fourF
        handlers["variance"] = fourF
        handlers["MSscore"] = lambda x: "%d" % x

        print(delim.join(cols), file=f)

        def fmt(rowData, colName):
            if colName not in rowData:
                return ""

            if colName in handlers:
                return handlers[colName](rowData[colName])
            else:
                return str(rowData[colName])

        try:
            while True:
                # Pull a list of record in from the producer
                itemList = (yield)

                for item in itemList:
                    if "signal" in item:
                        values = [fmt(item, col) for col in cols]
                        print(delim.join(values), file=f)

        except GeneratorExit:
            f.close()
            return
        except Exception as e:
            print(e)


    @consumer
    def csvConsumer(self, filename):
        """
        Consume IPD summary rows and write them to csv
        """

        # Open the csv file
        f = self.openWriteHandle(filename)
        delim = ","

        if self.options.control is None:

            # Columns for in-silico control
            if self.options.methylFraction:
                cols = ["refName", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction", "ipdRatio", "coverage", FRAC, FRAClow, FRACup]
            else:
                if self.options.useLDA:
                    # FIXME: For testing LDA model, to look at LDA scores in csv output (run without --methylFraction or --control):
                    cols = ["refName", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction", "ipdRatio", "coverage", "Ca5C"]
                else:
                    cols = ["refName", "tpl", "strand", "base", "score", "tMean", "tErr", "modelPrediction", "ipdRatio", "coverage"]

        else:
            # Columns for case-control
            if self.options.methylFraction:
                cols = ["refName", "tpl", "strand", "base", "score", "pvalue", "caseMean", "controlMean", "caseStd", "controlStd", "ipdRatio", "testStatistic", "coverage", "controlCoverage", "caseCoverage", FRAC, FRAClow, FRACup]
            else:
                cols = ["refName", "tpl", "strand", "base", "score", "pvalue", "caseMean", "controlMean", "caseStd", "controlStd", "ipdRatio", "testStatistic", "coverage", "controlCoverage", "caseCoverage"]

        # Special cases for formatting columns of the csv
        handlers = dict()
        threeF = lambda x: "%.3f" % x

        handlers["refName"] = lambda x: "\"%s\"" % x

        handlers["tpl"] = lambda x: str(x.item() + 1)
        handlers["score"] = lambda x: "%d" % x

        handlers["tMean"] = threeF
        handlers["modelPrediction"] = threeF
        handlers["caseMean"] = threeF
        handlers["controlMean"] = threeF
        handlers["ipdRatio"] = threeF
        handlers["pvalue"] = lambda x: "%.3e" % x

        handlers["controlStd"] = threeF
        handlers["controlStd"] = threeF
        handlers["tErr"] = threeF

        # FIXME: remove this line later:
        handlers["Ca5C"] = threeF

        handlers[FRAC] = threeF
        handlers[FRAClow] = threeF
        handlers[FRACup] = threeF

        print(delim.join(cols), file=f)

        def fmt(rowData, colName):
            if colName not in rowData:
                return ""

            if colName in handlers:
                return handlers[colName](rowData[colName])
            else:
                return str(rowData[colName])

        try:
            while True:
                # Pull a list of record in from the producer
                itemList = (yield)

                for item in itemList:
                    values = [fmt(item, col) for col in cols]
                    print(delim.join(values), file=f)

        except GeneratorExit:
            f.close()
            return
        except Exception as e:
            print(e)

    @consumer
    def bigWigConsumer(self, filename):
        import pyBigWig
        records = []
        records_by_pos = defaultdict(list)
        ranges = {}
        BaseInfo = namedtuple("BaseInfo", ("seqid", "pos", "sense", "ipd"))
        try:
            while True:
                chunk = (yield)
                if len(chunk) == 0:
                    continue
                # Fill out the ipd observations into the dataset
                for x in chunk:
                    pos = int(x['tpl']) + 1
                    seqid = x['refName']
                    ranges.setdefault(seqid, (sys.maxsize, 0))
                    ranges[seqid] = (min(ranges[seqid][0], pos),
                                     max(ranges[seqid][1], pos+1))
                    rec = BaseInfo(
                        seqid=seqid,
                        pos=pos,
                        sense=int(x['strand']),
                        ipd=float(x['ipdRatio']))
                    records.append(rec)
                    records_by_pos[(rec.seqid, rec.pos)].append(rec)
        except GeneratorExit:
            records.sort(key=lambda x: x.pos)
            records.sort(key=lambda x: x.seqid)
            regions = [(s, ranges[s][1]-1) for s in sorted(ranges.keys())]
            if len(regions) == 0:
                with open(filename, "wb") as _:
                    return
            bw = pyBigWig.open(filename, "w")
            bw.addHeader(regions)
            k = 0
            seqids = []
            starts = []
            ends = []
            ipd_enc = []
            # records are not necessarily consecutive or two per base!
            have_pos = set()
            def encode_ipds(plus, minus):
                enc = lambda x: min(65535, int(round(100*x)))
                return float(enc(minus) + 65536*enc(plus))
            for rec in records:
                if (rec.seqid, rec.pos) in have_pos:
                    continue
                have_pos.add((rec.seqid, rec.pos))
                strand_records = records_by_pos[(rec.seqid, rec.pos)]
                if len(strand_records) == 2:
                    rec_minus = strand_records[k] if strand_records[k].sense else strand_records[k + 1]
                    rec_plus = strand_records[k + 1] if strand_records[k].sense else strand_records[k]
                    assert rec_plus.pos == rec_minus.pos, (rec_plus, rec_minus)
                    seqids.append(rec_plus.seqid)
                    starts.append(rec_plus.pos-1)
                    ends.append(rec_plus.pos)
                    ipd_enc.append(encode_ipds(rec_plus.ipd, rec_minus.ipd))
                else:
                    seqids.append(rec.seqid)
                    starts.append(rec.pos-1)
                    ends.append(rec.pos)
                    if rec.sense == 0:
                        ipd_enc.append(encode_ipds(rec.ipd, 0))
                    else:
                        ipd_enc.append(encode_ipds(0, rec.ipd))
            log.info("Writing records for {n} bases".format(n=len(seqids)))
            bw.addEntries(seqids, starts, ends=ends, values=ipd_enc)
            bw.close()
            return

    def openWriteHandle(self, filename):
        if filename[-2:] == 'gz':
            import gzip
            fileobj = gzip.GzipFile(filename, mode="w", compresslevel=3)
        else:
            fileobj = open(filename, "w", 2 << 15)

        return fileobj

    @consumer
    def pickleConsumer(self, fileName):
        """
        Consume IPD summary rows and pickle to a 'None' terminated stream
        """

        f = open(fileName, "w")
        pickleStream = pickle.Pickler(f)

        try:
            while True:
                # Pickle a record
                n = (yield)
                pickleStream.dump(n)
                pickleStream.clear_memo()

        except GeneratorExit:
            # Write an end sentinel to the pickle stream
            pickleStream.dump(None)
            f.close()
            return



    def makeGffRecord(self, siteObs):
        """
        Convert the internal site observation object into a GFF entry
        """
        # Some useful attributes about the observation
        # - cognate base
        # - context snippet
        # - ipd ratio
        # - coverage
        snippet = self.snippetFunc(siteObs['tpl'], siteObs['strand'])
        attributes = [('coverage', siteObs['coverage']),
                      ('context', snippet),
                      ('IPDRatio', siteObs['ipdRatio'])]

        # Base of detected mod -- single position, closed,open
        # interval.
        # Note -- internally the tool uses 0-based reference
        # coordinates, however in gff the template indices are
        # 1-based.  Make that adjustment here.
        # On start vs. end: My reading of the gff spec
        # (http://www.sequenceontology.org/resources/gff3.html) says
        # to me that 1-base long feature (e.g. a modified base) should
        # have start + 1 == end, and 0-base long features
        # (e.g. insertions) should have start == end. This is not the
        # convention that Marco has apdopted in SMRTView, or the
        # convention that EviCons originally used.  We will adopt
        # their convention here, for now.
        start = siteObs['tpl'] + 1
        end = siteObs['tpl'] + 1

        if 'motif' in siteObs:
            attributes.append(('motif', "%s" % siteObs['motif']))

        if 'id' in siteObs:
            attributes.append(('id', "%s" % siteObs['id']))

        if self.options.methylFraction and FRAC in siteObs:
            attributes.append(('frac', "%.3f" % siteObs[FRAC]))
            attributes.append(('fracLow', "%.3f" % siteObs[FRAClow]))
            attributes.append(('fracUp', "%.3f" % siteObs[FRACup]))

        if 'modificationScore' in siteObs:
            # Report the QV from the modification identification module as a special tag
            attributes.append(('identificationQv', "%d" % int(round(siteObs['modificationScore']))))

        if 'modification' in siteObs:

            if siteObs['modification'] == '.':
                recordType = 'modified_base'

            elif siteObs['modification'] == 'nMd':
                recordType = '.'

            else:
                # if we have an identified mod, use it; otherwise use the old generic term
                recordType = siteObs['modification']

        else:
            recordType = 'modified_base'

        refName = siteObs['refName']
        score = int(round(siteObs['score']))
        strand = '+' if siteObs['strand'] == 0 else '-'

        return Gff3Record(refName, start, end,
                          type=recordType,
                          score=score,
                          strand=strand,
                          source='kinModCall',
                          attributes=attributes)

    @consumer
    def gffConsumer(self, filename):
        """
        Consume IPD summary rows, filter them and write to GFF
        """

        #f = file(filename, 'w', 2<<15)
        f = self.openWriteHandle(filename)
        gff = GffWriter(f)

        # write headers describing the program that generated the data
        gff.writeHeader('##source ipdSummary v2.0')
        gff.writeHeader('##source-commandline %s' % self.options.cmdLine)

        # Write the reference renaming info into the gff headers ala evicons
        for entry in self.refInfo:
            gff.writeHeader("##sequence-region %s 1 %d"
                            % (entry.Name, entry.Length))

        minScore = -10 * math.log10(self.options.pvalue)
        snippetRef = -1
        try:
            while True:
                # Pull a record in from the
                siteObsList = (yield)

                for siteObs in siteObsList:
                    # self.snippetFunc is a function that return a reference snippet given a template position and a strand
                    if snippetRef != siteObs['refId']:
                        self.snippetFunc = self.ipdModel.snippetFunc(siteObs['refId'], 20, 20)
                        snippetRef = siteObs['refId']

                    # Two cases for gff entries:
                    # 1. 'Identified modification' - will have a 'modification' key
                    #     - use the modification name as the gff event type
                    #     - use 'modificationScore' for the gff score
                    # 2. Detected - no 'modification' key
                    #     - use 'modified_base' as the event type
                    #     - use the single site 'score' property as the gff score
                    #     - do not put this kind into the gff if it contains the a 'offTargetPeak' tag

                    if siteObs['coverage'] > self.options.minCoverage:
                        # Case 1
                        if 'modification' in siteObs and siteObs['modification'] != '.':
                            gff.writeRecord(self.makeGffRecord(siteObs))

                        # Case 2
                        elif siteObs['score'] > minScore and 'offTargetPeak' not in siteObs:
                            gff.writeRecord(self.makeGffRecord(siteObs))

                    # FIXME: Try not filtering:
                    # gff.writeRecord(self.makeGffRecord(siteObs))

        except GeneratorExit:
            f.close()
            return


    def makeM5CgffRecord(self, siteObs):


        start = siteObs['tpl'] + 1
        end = siteObs['tpl'] + 1

        attributes = [('coverage', siteObs['coverage']),
                      ('IPDRatio', siteObs['ipdRatio'])]

        recordType = 'CG'
        refName = siteObs['refId']
        score = "%.3f" % siteObs['Ca5C']
        strand = '+' if siteObs['strand'] == 0 else '-'

        return Gff3Record(refName, start, end,
                          type=recordType,
                          score=score,
                          strand=strand,
                          source='kinModCall',attributes=attributes)


    @consumer
    def m5CgffConsumer( self, filename ):
    
        f = self.openWriteHandle( filename )        
        gff = GffWriter( f )

              
        # write headers describing the program that generated the data        
        gff.writeHeader('##source ipdSummary v2.0')        
        gff.writeHeader('##source-commandline %s' % self.options.cmdLine)
        
        # Write the reference renaming info into the gff headers ala evicons        
        # for entry in self.refInfo:        
        #     gff.writeHeader("##sequence-region %s 1 %d"
        #                     % (entry.Name, entry.Length))
                
        try:
            while True:
                # Pull in a single record?
                siteObsList = (yield)

                for siteObs in siteObsList:
                    if 'Ca5C' in siteObs and siteObs['strand'] == 0:
                        gff.writeRecord( self.makeM5CgffRecord( siteObs ) )

        except GeneratorExit:
            f.close()
            return




    def onStart(self):

        # Spec for what kinds of output files we can generate.
        # Entry format is (<option field name>, <extension>, <writer consumer function>)
        fileSpec = [
            ('m5Cgff', 'm5C.gff', self.m5CgffConsumer),
            ('gff', 'gff', self.gffConsumer),
            ('csv', 'csv', self.csvConsumer),
            ('bigwig', 'bw', self.bigWigConsumer),
            ('ms_csv', 'ms.csv', self.msCsvConsumer),
            ('pickle', 'pickle', self.csvConsumer),
        ]

        sinkList = []

        # Go through the possible output file types and
        # determine if they should be output
        for (fileType, ext, func) in fileSpec:
            name = None

            # The 'outfile argument causes all outputs to be generated
            if self.options.outfile:
                if ext == "bw":
                    try:
                        import pyBigWig
                    except ImportError:
                        pass
                    else:
                        name = self.options.outfile + '.' + ext
                else:
                    name = self.options.outfile + '.' + ext

            # Individual outputs can specified - these filename override the default
            if self.options.__getattribute__(fileType):
                name = self.options.__getattribute__(fileType)

            if name:
                sinkList.append(func(name))

        self.sinkList = sinkList

    def onResult(self, resultChunk):
        for sink in self.sinkList:
            sink.send(resultChunk)

    def onFinish(self):
        for sink in self.sinkList:
            sink.close()
