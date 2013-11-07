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

import cProfile, logging, os.path
from multiprocessing import Process
from multiprocessing.process import current_process
from threading import Thread, Event

from pbcore.io import CmpH5Reader

class Worker(object):
    """
    Base class for worker processes that read reference coordinates
    from the task queue, perform variant calling, then push results
    back to another queue, to be written to a GFF file by another
    process.

    All tasks that are O(genome length * coverage depth) should be
    distributed to Worker processes, leaving the ResultCollector
    process only O(genome length) work to do.
    """

    def __init__(self, options, workQueue, resultsQueue):
        self.options = options
        self.daemon = True
        self._workQueue = workQueue
        self._resultsQueue = resultsQueue

    def _run(self):
        logging.info("Worker %s (PID=%d) started running" % (self.name, self.pid))

        self.caseCmpH5 = CmpH5Reader(self.options.infile)

        if not self.options.control is None:
            # We have a cmp.h5 with control vales -- load that cmp.h5
            self.controlCmpH5 = CmpH5Reader(self.options.control)
        else:
            self.controlCmpH5 = None

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
                logging.info("Got chunk: (%s, %s) -- Process: %s" % (chunkId, str(datum), current_process()))
                result = self.onChunk(datum)

                logging.debug("Process %s: putting result." % current_process())
                self._resultsQueue.put((chunkId,result))
                self._workQueue.task_done()

        self.onFinish()

        logging.info("Process %s (PID=%d) done; exiting." % (self.name, self.pid))

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


    #==
    # Begin overridable interface
    #==

    def onStart(self):
        pass

    def onChunk(self, target):
        """
        This function is the heart of the matter.

        referenceWindow, alnHits -> result
        """
        pass

    def onFinish(self):
        pass


class WorkerProcess(Worker, Process):
    """Worker that executes as a process."""
    def __init__(self, *args):
        Process.__init__(self)
        super(WorkerProcess,self).__init__(*args)
        self.daemon = True

    def _lowPriority(self):
        """
        Set the priority of the process to below-normal.
        """
        import sys
        try:
            sys.getwindowsversion()
        except:
            isWindows = False
        else:
            isWindows = True

        if isWindows:
            # Based on:
            #   "Recipe 496767: Set Process Priority In Windows" on ActiveState
            #   http://code.activestate.com/recipes/496767/
            import win32api,win32process,win32con

            pid = win32api.GetCurrentProcessId()
            handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
            win32process.SetPriorityClass(handle, win32process.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            os.nice(10)

    def isTerminated(self):
        return False


class WorkerThread(Worker, Thread):
    """Worker that executes as a thread (for debugging purposes only)."""
    def __init__(self, *args):
        Thread.__init__(self)
        super(WorkerThread,self).__init__(*args)
        self._stop = Event()
        self.daemon = True
        self.exitcode = 0

    def terminate(self):
        self._stop.set()

    def isTerminated(self):
        return self._stop.isSet()

    @property
    def pid(self):
        return -1

    def _lowPriority(self):
        pass
