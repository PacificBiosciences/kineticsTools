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

import cProfile
import logging
import os.path
import copy
from multiprocessing import Process
from multiprocessing.process import current_process
import warnings

from future.moves.urllib.parse import urlparse
import numpy as np

from pbcore.io import AlignmentSet
from pbcore.io.opener import (openAlignmentFile, openIndexedAlignmentFile)


# FIXME this should ultimately go somewhere else.  actually, so should the
# rest of this module.
def _openFiles(self, refFile=None, sharedIndices=None):
    """
    Hack to enable sharing of indices (but not filehandles!) between dataset
    instances.
    """
    log = logging.getLogger()
    log.debug("Opening resources")
    for k, extRes in enumerate(self.externalResources):
        location = urlparse(extRes.resourceId).path
        sharedIndex = None
        if sharedIndices is not None:
            sharedIndex = sharedIndices[k]
        try:
            resource = openIndexedAlignmentFile(
                location,
                referenceFastaFname=refFile,
                sharedIndex=sharedIndex)
        except (IOError, ValueError):
            log.info("pbi file missing for {f}, operating with "
                     "reduced speed and functionality".format(
                         f=location))
            resource = openAlignmentFile(location,
                                         referenceFastaFname=refFile)
        if len(resource) == 0:
            log.warn("{f} has no mapped reads".format(f=location))
        else:
            self._openReaders.append(resource)
    if len(self._openReaders) == 0:
        raise IOError("No mapped reads found")
    log.debug("Done opening resources")


def _reopen(self):
    """
    Force re-opening of underlying alignment files, preserving the
    reference and indices if present, and return a copy of the
    AlignmentSet.  This is a workaround to allow us to share the index
    file(s) already loaded in memory while avoiding multiprocessing
    problems related to .bam files.
    """
    refFile = self._referenceFile
    newSet = copy.deepcopy(self)
    newSet._referenceFastaFname = refFile
    indices = [f.index for f in self.resourceReaders()]
    self.close()
    _openFiles(newSet, refFile=refFile, sharedIndices=indices)
    return newSet


class WorkerProcess(Process):

    """
    Base class for worker processes that read reference coordinates
    from the task queue, perform variant calling, then push results
    back to another queue, to be written to a GFF file by another
    process.

    All tasks that are O(genome length * coverage depth) should be
    distributed to Worker processes, leaving the ResultCollector
    process only O(genome length) work to do.
    """

    def __init__(self, options, workQueue, resultsQueue,
                 sharedAlignmentSet=None):
        Process.__init__(self)
        self.options = options
        self.daemon = True
        self._workQueue = workQueue
        self._resultsQueue = resultsQueue
        self._sharedAlignmentSet = sharedAlignmentSet

    def _run(self):
        logging.info("Worker %s (PID=%d) started running" %
                     (self.name, self.pid))
        if self._sharedAlignmentSet is not None:
            # XXX this will create an entirely new AlignmentSet object, but
            # keeping any indices already loaded into memory
            self.caseAlignments = _reopen(self._sharedAlignmentSet)
            # `self._sharedAlignmentSet.close()
            self._sharedAlignmentSet = None
        else:
            warnings.warn("Shared AlignmentSet not used")
            self.caseAlignments = AlignmentSet(self.options.infile,
                                               referenceFastaFname=self.options.reference)

        self.controlAlignments = None
        if not self.options.control is None:
            self.controlAlignments = AlignmentSet(self.options.control,
                                                  referenceFastaFname=self.options.reference)

        if self.options.randomSeed is None:
            np.random.seed(42)
        self.onStart()

        while True:
            if self.isTerminated():
                break

            chunkDesc = self._workQueue.get()
            if chunkDesc is None:
                # Sentinel indicating end of input.  Place a sentinel
                # on the results queue and end this worker process.
                self._resultsQueue.put(None)
                self._workQueue.task_done()
                break
            else:
                (chunkId, datum) = chunkDesc
                logging.info("Got chunk: (%s, %s) -- Process: %s" %
                             (chunkId, str(datum), current_process()))
                result = self.onChunk(
                    datum)  # pylint: disable=assignment-from-none

                logging.debug("Process %s: putting result." %
                              current_process())
                self._resultsQueue.put((chunkId, result))
                self._workQueue.task_done()

        self.onFinish()

        logging.info("Process %s (PID=%d) done; exiting." %
                     (self.name, self.pid))

    def run(self):
        # Make the workers run with lower priority -- hopefully the results writer will win
        # It is single threaded so it could become the bottleneck
        self._lowPriority()

        if self.options.doProfiling:
            cProfile.runctx("self._run()",
                            globals=globals(),
                            locals=locals(),
                            filename="profile-%s.out" % self.name)
        else:
            self._run()

    # ==
    # Begin overridable interface
    # ==
    def onStart(self):
        pass

    def onChunk(self, target):
        """
        This function is the heart of the matter.

        referenceWindow, alnHits -> result
        """
        return None

    def onFinish(self):
        pass

    def isTerminated(self):
        return False

    def _lowPriority(self):
        """
        Set the priority of the process to below-normal.
        """
        os.nice(10)
