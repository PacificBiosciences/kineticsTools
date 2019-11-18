# FIXME all of this belongs somewhere else (probably either pbcore.io.dataset
# or a future base module for resequencing apps)

from collections import namedtuple
import itertools
import math
import re
import os

from pbcore.io import ReferenceSet

# FIXME pbcore keys contigs by name, but kineticsTools usually keys by ID
ReferenceWindow = namedtuple(
    "ReferenceWindow", ["refId", "refName", "start", "end"])


class ReferenceUtils():

    @staticmethod
    def loadReferenceContigs(referencePath, alignmentSet, windows=None):
        # FIXME we should get rid of this entirely, but I think it requires
        # fixing the inconsistency in how contigs are referenced here versus in
        # pbcore.io
        """
        Load the reference contigs, and tag each one with the ref.alignmentID it
        was assigned in the alignment file(s).  Return a list of contigs,
        which are used to set up IpdModel.
        """

        # Read contigs from FASTA file (or XML dataset)
        refReader = ReferenceSet(referencePath)
        contigs = []
        if windows is not None:
            refNames = set([rw.refName for rw in windows])
            for contig in refReader:
                if contig.id in refNames:
                    contigs.append(contig)
        else:
            contigs.extend([x for x in refReader])
        contigDict = dict([(x.id, x) for x in contigs])

        # initially each contig has an id of None -- this will be overwritten with the id from the alignments if there are any
        # reads mapped to it.
        for x in contigs:
            x.alignmentID = None

        # Mark each contig with it's ID from the alignments - match them up using MD5s
        for x in alignmentSet.referenceInfoTable:
            if x.FullName in contigDict:
                contigDict[x.FullName].alignmentID = x.ID

        return contigs

    @staticmethod
    def referenceWindowsFromAlignment(ds, refInfoLookup):
        return [ReferenceWindow(refId=refInfoLookup(w[0]).ID,
                                refName=w[0],
                                start=w[1],
                                end=w[2]) for w in ds.refWindows]

    @staticmethod
    def parseReferenceWindow(s, refInfoLookup):
        if s is None:
            return None
        m = re.match("(.*):(.*)-(.*)", s)
        if m:
            refContigInfo = refInfoLookup(m.group(1))
            refId = refContigInfo.ID
            refName = refContigInfo.Name
            refStart = int(m.group(2))
            refEnd = min(int(m.group(3)), refContigInfo.Length)
        else:
            refContigInfo = refInfoLookup(s)
            refId = refContigInfo.ID
            refName = refContigInfo.Name
            refStart = 0
            refEnd = refContigInfo.Length
        return ReferenceWindow(refId=refId, refName=refName, start=refStart,
                               end=refEnd)

    @staticmethod
    def createReferenceWindows(refInfo):
        return [ReferenceWindow(refId=r.ID,
                                refName=r.Name,
                                start=0,
                                end=r.Length) for r in refInfo]

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
                return (x//chunk)*chunk

            def alignUp(chunk, x):
                return int(math.ceil(float(x)/chunk)*chunk)

            start, end = bounds
            roundStart = alignDown(stride, start)
            roundEnd = alignUp(stride, end)

            for s in range(roundStart, roundEnd, stride):
                roundWin = (s, s + stride)
                yield intersection(bounds, roundWin)

        for (s, e) in enumerateIntervals((referenceWindow.start,
                                          referenceWindow.end), referenceStride):
            yield ReferenceWindow(refId=referenceWindow.refId,
                                  refName=referenceWindow.refName,
                                  start=s, end=e)

    @staticmethod
    def loadAlignmentChemistry(alignmentSet):
        chems = alignmentSet.sequencingChemistry
        chemCounts = {k: len(list(v)) for k, v in itertools.groupby(chems)}
        majorityChem = max(chemCounts, key=chemCounts.get)
        return majorityChem
