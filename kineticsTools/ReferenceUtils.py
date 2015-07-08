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

# FIXME most of this belongs somewhere else

import os, itertools, re, math
from collections import namedtuple

from pbcore.io import AlignmentSet, ReferenceSet

ReferenceWindow = namedtuple("ReferenceWindow", ["refId", "refName", "start", "end"])

class ReferenceUtils():

    @staticmethod
    def loadReferenceContigs(referencePath, alignmentPath):
        """Load the reference contigs, and tag each one with the ref.cmpH5ID it was assigned in the cmp.h5 file.  Return a list of contigs, which are used to set up IpdModel"""

        # Read contigs from FASTA file (or XML dataset)
        refReader = ReferenceSet(referencePath)
        contigs = [x for x in refReader]
        contigDict = dict([(x.id, x) for x in contigs])

        # Read reference info table from cmp.h5
        (refInfoTable, movieInfoTable) = ReferenceUtils.loadAlignmentTables(alignmentPath)

        # initially each contig has an id of None -- this will be overwritten with the id from the cmp.h5, if there are any
        # reads mapped to it.
        for x in contigs:
            x.cmph5ID = None

        # Mark each contig with it's ID from the cmp.h5 - match them up using MD5s
        for x in refInfoTable:
            contigDict[x.FullName].cmph5ID = x.ID

        return contigs

    @staticmethod
    def parseReferenceWindow(s, refInfoLookup):
        if s is None:
            return None
        m = re.match("(.*):(.*)-(.*)", s)
        if m:
            refContigInfo = refInfoLookup(m.group(1))
            refId    = refContigInfo.ID
            refName  = refContigInfo.Name
            refStart = int(m.group(2))
            refEnd   = min(int(m.group(3)), refContigInfo.Length)
        else:
            refContigInfo = refInfoLookup(s)
            refId    = refContigInfo.ID
            refName  = refContigInfo.Name
            refStart = 0
            refEnd   = refContigInfo.Length
        return ReferenceWindow(refId=refId, refName=refName, start=refStart,
            end=refEnd)

    @staticmethod
    def createReferenceWindows(refInfo):
        return [ ReferenceWindow(refId=r.ID,
                    refName=r.Name,
                    start=0,
                    end=r.Length) for r in refInfo ]


    @staticmethod
    def enumerateChunks(referenceStride, referenceWindow):
        """
        Enumerate all work chunks on this reference contig (restricted to
        the windows, if provided).
        """
        def intersection(int1, int2):
            s1, e1 = int1
            s2, e2 = int2
            si, ei = max(s1, s2), min(e1, e2)
            if si < ei:
                return (si, ei)
            else:
                return None

        def enumerateIntervals(bounds, stride):
            """
            Enumerate windows of size "stride", attempting to align window
            boundaries on multiple of stride.
            """
            def alignDown(chunk, x):
                return (x/chunk)*chunk
            def alignUp(chunk, x):
                return int(math.ceil(float(x)/chunk)*chunk)

            start, end = bounds
            roundStart = alignDown(stride, start)
            roundEnd   = alignUp  (stride, end)

            for s in xrange(roundStart, roundEnd, stride):
                roundWin = (s, s + stride)
                yield intersection(bounds, roundWin)

        for (s, e) in enumerateIntervals((referenceWindow.start,
                referenceWindow.end), referenceStride):
            yield ReferenceWindow(refId=referenceWindow.refId,
                refName=referenceWindow.refName,
                start=s, end=e)


    @staticmethod
    def loadAlignmentTables(alignmentFile):
        """
        Load the alignments and get the ReferenceInfo table, in order to
        correctly number the contigs.
        """
        with AlignmentSet(alignmentFile) as ds:
            return ds.referenceInfoTable, ds.readGroupTable

    @staticmethod
    def loadAlignmentChemistry(alignmentFile):
        with AlignmentSet(alignmentFile) as ds:
            chems = ds.sequencingChemistry
            chemCounts = {k: len(list(v)) for k, v in itertools.groupby(chems)}
            majorityChem = max(chemCounts, key=chemCounts.get)
            return majorityChem
