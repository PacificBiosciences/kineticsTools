#!/usr/bin/env python
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
import functools
import gc
import itertools
import argparse

import os
import logging
import sys
import multiprocessing
import time
import threading
import numpy as np
import Queue
import traceback

from pbcore.io import CmpH5Reader

from pbtools.kineticsTools.KineticWorker import KineticWorkerThread, KineticWorkerProcess
from pbtools.kineticsTools.ResultWriter import KineticsWriter
from pbtools.kineticsTools.ipdModel import IpdModel
from pbtools.kineticsTools.ReferenceUtils import ReferenceUtils

# Version info
__p4revision__ = "$Revision: #1 $"
__p4change__ = "$Change: 100972 $"
revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum = int(__p4change__.strip("$").split(":")[-1])
__version__ = "2.2"


def _validateResource(func, p):
    """Basic func for validating files, dirs, etc..."""
    if func(p):
        return os.path.abspath(p)
    else:
        raise IOError("Unable to find {p}".format(p=p))


def _validateNoneOrResource(func, p):
    """
    Handle optional values. If a file or dir is explicitly provided, then
    it will validated.
    """
    if p is None:
        return p
    else:
        return _validateResource(func, p)


validateFile = functools.partial(_validateResource, os.path.isfile)
validateDir = functools.partial(_validateResource, os.path.isdir)

validateNoneOrFile = functools.partial(_validateNoneOrResource, os.path.isfile)
validateNoneOrDir = functools.partial(_validateNoneOrResource, os.path.isdir)


class KineticsToolsRunner(object):

    def __init__(self):
        desc = ['Tool for detecting DNA base-modifications from kinetic signatures',
                'Notes: For all command-line arguments, default values are listed in [].']
        description = '\n'.join(desc)

        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                              description=description,
                                              version=__version__)


        # Positional arguments:

        self.parser.add_argument('reference', type=validateFile,
                                 help='Path to reference FASTA file')

        self.parser.add_argument('infile',
                                 metavar='input.cmp.h5',
                                 type=validateFile,
                                 help='Input cmp.h5 filename')

        # Optional arguments:

        # Output options:

        self.parser.add_argument('--outfile',
                                 dest='outfile',
                                 default=None,
                                 help='Use this option to generate all possible output files. Argument here is the root filename of the output files.')

        self.parser.add_argument('--gff',
                                 dest='gff',
                                 default=None,
                                 help='Name of output GFF file')

        self.parser.add_argument('--csv',
                                 dest='csv',
                                 default=None,
                                 help='Name of output CSV file.')


        self.parser.add_argument('--csv_h5',
                                 dest='csv_h5',
                                 default=None,
                                 help='Name of csv output to be written in hdf5 format.')

        self.parser.add_argument('--pickle',
                                 dest='pickle',
                                 default=None,
                                 help='Name of output pickle file.')

        self.parser.add_argument('--summary_h5',
                                 dest='summary_h5',
                                 default=None,
                                 help='Name of output summary h5 file.')


        self.parser.add_argument('--ms_csv',
                                 dest='ms_csv',
                                 default=None,
                                 help='Multisite detection CSV file.')


        # Calculation options:


        self.parser.add_argument('--control',
                                 dest='control',
                                 default=None,
                                 type=validateNoneOrFile,
                                 help='cmph.h5 file containing a control sample. Tool will perform a case-control analysis')

        self.parser.add_argument('--identify',
                                 dest='identify',
                                 default=False,
                                 help='Identify modification types. Comma-separated list of know modification types. Current options are: m6A, m4C, m5C_TET. Cannot be used with --control')


        self.parser.add_argument("--methylFraction",
                                 action="store_true",
                                 dest="methylFraction",
                                 default=False,
                                 help="In the --identify mode, add --methylFraction to command line to estimate the methylated fraction, along with 95%% confidence interval bounds.")

        # Temporary addition to test LDA for Ca5C detection:
        self.parser.add_argument('--useLDA',
                                 action="store_true",
                                 dest='useLDA',
                                 default=False,
                                 help='Set this flag to debug LDA for m5C/Ca5C detection')



        # Parameter options:

        self.parser.add_argument('--paramsPath',
                                 dest='paramsPath',
                                 default=None,
                                 type=validateNoneOrDir,
                                 help='Directory containing in-silico trained model for each chemistry')

        self.parser.add_argument("--maxLength",
                                 default=3e12,
                                 type=int,
                                 help="Maximum number of bases to process per contig")

        self.parser.add_argument('--minCoverage',
                                 dest='minCoverage',
                                 default=3,
                                 type=int,
                                 help='Minimum coverage required to call a modified base')

        self.parser.add_argument('--maxQueueSize',
                                 dest='maxQueueSize',
                                 default=20,
                                 type=int,
                                 help='Max Queue Size')

        self.parser.add_argument('--maxCoverage',
                                 dest='maxCoverage',
                                 type=int, default=-1,
                                 help='Maximum coverage to use at each site')

        self.parser.add_argument('--mapQvThreshold',
                                 dest='mapQvThreshold',
                                 type=float,
                                 default=-1.0)

        self.parser.add_argument('--pvalue',
                                 dest='pvalue',
                                 default=0.01,
                                 type=float,
                                 help='p-value required to call a modified base')

        self.parser.add_argument('--ipdModel',
                                 dest='ipdModel',
                                 default=None,
                                 help='Alternate synthetic IPD model HDF5 file')

        self.parser.add_argument('--modelIters',
                                 dest='modelIters',
                                 type=int,
                                 default=-1,
                                 help='[Internal] Number of GBM model iteration to use')

        self.parser.add_argument('--cap_percentile',
                                 dest='cap_percentile',
                                 type=float,
                                 default=99.0,
                                 help='Global IPD percentile to cap IPDs at')


        self.parser.add_argument("--methylMinCov",
                                 type=int,
                                 dest='methylMinCov',
                                 default=10,
                                 help="Do not try to estimate methylFraction unless coverage is at least this.")

        self.parser.add_argument("--identifyMinCov",
                                 type=int,
                                 dest='identifyMinCov',
                                 default=5,
                                 help="Do not try to identify the modification type unless coverage is at least this.")


        # Computation management options:

        self.parser.add_argument("--refId",
                                 type=int,
                                 dest='refId',
                                 default=-1,
                                 help="Specify a single reference index (beginning with 0) rather than looping through all")


        self.parser.add_argument('--numWorkers',
                                 dest='numWorkers',
                                 default=-1,  # Defaults to using all logical CPUs
                                 type=int,
                                 help='Number of thread to use (-1 uses all logical cpus)')

        # Debugging help options:

        self.parser.add_argument("--threaded", "-T",
                                 action="store_true",
                                 dest="threaded",
                                 default=False,
                                 help="Run threads instead of processes (for debugging purposes only)")

        self.parser.add_argument("--profile",
                                 action="store_true",
                                 dest="doProfiling",
                                 default=False,
                                 help="Enable Python-level profiling (using cProfile).")

        self.parser.add_argument('--pdb',
                                 action='store_true',
                                 dest="usePdb",
                                 default=False,
                                 help="Enable dropping down into pdb debugger if an Exception is raised.")




    def parseArgs(self):
        self.args = self.parser.parse_args()

    def start(self):
        self.parseArgs()
        self.validateArgs()
        return self.run()

    def getVersion(self):
        return __version__

    def validateArgs(self):
        if not os.path.exists(self.args.infile):
            self.parser.error('input.cmp.h5 file provided does not exist')

        if self.args.identify and self.args.control:
            self.parser.error('--control and --identify are mutally exclusive. Please choose one or the other')

        # if self.args.methylFraction and not self.args.identify:
        #    self.parser.error('Currently, --methylFraction only works when the --identify option is specified.')

    def run(self):

        # Figure out what modifications to identify
        mods = self.args.identify
        modsToCall = []
        if mods:
            items = mods.split(",")

            if 'm6A' in items:
                modsToCall.append('H')

            if 'm4C' in items:
                modsToCall.append('J')

            if 'm5C_TET' in items:
                modsToCall.append('K')

            self.args.identify = True
            self.args.modsToCall = modsToCall

        self.options = self.args
        self.options.cmdLine = " ".join(sys.argv)
        self._workers = []

        # Log generously
        stdOutHandler = logging.StreamHandler(sys.stdout)
        logFormat = '%(asctime)s [%(levelname)s] %(message)s'
        logging.basicConfig(level=logging.INFO, format=logFormat)

        if self.args.doProfiling:
            cProfile.runctx("self._mainLoop()",
                            globals=globals(),
                            locals=locals(),
                            filename="profile.out")

        else:
            try:
                ret = self._mainLoop()
            finally:
                # Be sure to shutdown child processes if we get an exception on the main thread
                if not self.args.threaded:
                    for w in self._workers:
                        if w.is_alive():
                            w.terminate()

            return ret

    def _initQueues(self):
        if self.options.threaded:
            # Work chunks are created by the main thread and put on this queue
            # They will be consumed by KineticWorker threads, stored in self._workers
            self._workQueue = Queue.Queue(self.options.maxQueueSize)

            # Completed chunks are put on this queue by KineticWorker threads
            # They are consumed by the KineticsWriter process
            self._resultsQueue = multiprocessing.JoinableQueue(self.options.maxQueueSize)
        else:
            # Work chunks are created by the main thread and put on this queue
            # They will be consumed by KineticWorker threads, stored in self._workers
            self._workQueue = multiprocessing.JoinableQueue(self.options.maxQueueSize)

            # Completed chunks are put on this queue by KineticWorker threads
            # They are consumed by the KineticsWriter process
            self._resultsQueue = multiprocessing.JoinableQueue(self.options.maxQueueSize)

    def _launchSlaveProcesses(self):
        """
        Launch a group of worker processes (self._workers), the queue
        (self._workQueue) that will be used to send them chunks of
        work, and the queue that will be used to receive back the
        results (self._resultsQueue).

        Additionally, launch the result collector process.
        """
        availableCpus = multiprocessing.cpu_count()
        logging.info("Available CPUs: %d" % (availableCpus,))
        logging.info("Requested worker processes: %d" % (self.options.numWorkers,))

        # Use all CPUs if numWorkers < 1
        if self.options.numWorkers < 1:
            self.options.numWorkers = availableCpus

        # Warn if we make a bad numWorker argument is used
        if self.options.numWorkers > availableCpus:
            logging.warn("More worker processes requested (%d) than CPUs available (%d);"
                         " may result in suboptimal performance."
                         % (self.options.numWorkers, availableCpus))

        self._initQueues()

        if self.options.threaded:
            self.options.numWorkers = 1
            WorkerType = KineticWorkerThread
        else:
            WorkerType = KineticWorkerProcess

        # Launch the worker processes
        self._workers = []
        for i in xrange(self.options.numWorkers):
            p = WorkerType(self.options, self._workQueue, self._resultsQueue, self.ipdModel)
            self._workers.append(p)
            p.start()
        logging.info("Launched worker processes.")

        # Launch result collector
        self._resultCollectorProcess = KineticsWriter(self.options, self._resultsQueue, self.refInfo, self.ipdModel)
        self._resultCollectorProcess.start()
        logging.info("Launched result collector process.")

        # Spawn a thread that monitors worker threads for crashes
        self.monitoringThread = threading.Thread(target=monitorChildProcesses, args=(self._workers + [self._resultCollectorProcess],))
        self.monitoringThread.start()

    def _queueChunksForReference(self, refInfo):
        """
        Compute the chunk extents and queue up the work for a single reference
        """

        # Number of hits on current reference
        refGroupId = refInfo.ID
        numHits = (self.cmph5.RefGroupID == refGroupId).sum()

        # Don't process reference groups with 0 hits.  They may not exist?
        if numHits == 0:
            return

        # Maximum chunk size (set no larger than 1Mb for now)
        MAX_BLOCK_SIZE = 25000

        # Maximum number of hits per chunk
        MAX_HITS = 5000
        nBases = min(refInfo.Length, self.args.maxLength)

        # Adjust numHits if we are only doing part of the contig
        numHits = (numHits * nBases) / refInfo.Length

        nBlocks = max([numHits / MAX_HITS, nBases / (MAX_BLOCK_SIZE - 1) + 1])

        # Including nBases / (MAX_BLOCK_SIZE - 1) + 1 in nBlocks calculation:
        # E. coli genome: this should be ~ 10.
        # Human genome: ought to be largest & is meant to ensure that blockSize < MAX_BLOCK_SIZE.

        # Block layout
        blockSize = min(nBases, max(nBases / nBlocks + 1, 1000))
        blockStarts = np.arange(0, nBases, step=blockSize)
        blockEnds = blockStarts + blockSize
        blocks = zip(blockStarts, blockEnds)

        logging.info("Queueing chunks for ref: %d.  NumReads: %d, Block Size: %d " % (refGroupId, numHits, blockSize))

        # Queue up work blocks
        for block in blocks:
            # NOTE! The format of a work chunk is (refId <int>, refStartBase <int>, refEndBase <int>)
            chunk = (refInfo.ID, block[0], block[1])
            self._workQueue.put((self.workChunkCounter, chunk))
            self.workChunkCounter += 1

            if self.workChunkCounter % 10 == 0:
                logging.info("Queued chunk: %d.  Chunks in queue: %d" % (self.workChunkCounter, self._workQueue.qsize()))

    def loadReferenceAndModel(self, referencePath, cmpH5Path):

        # Load the reference contigs - annotated with their refID from the cmp.h5
        contigs = ReferenceUtils.loadReferenceContigs(referencePath, cmpH5Path)

        # Read reference info table from cmp.h5
        (refInfoTable, movieInfoTable) = ReferenceUtils.loadCmpH5Tables(cmpH5Path)
        self.refInfo = refInfoTable

        # There are three different ways the ipdModel can be loaded.
        # In order of precedence they are:
        # 1. Explicit path passed to --ipdModel
        # 2. Path to parameter bundle, model selected using the /MovieInfo/SequencingChemistry tags
        # 3. Fall back to built-in model.

        # By default, use built-in model
        ipdModel = None

        if self.args.ipdModel:
            ipdModel = self.args.ipdModel
            logging.info("Using passed in ipd model: %s" % self.args.ipdModel)
            if not os.path.exists(self.args.ipdModel):
                logging.error("Couldn't find model file: %s" % self.args.ipdModel)
                sys.exit(1)
        elif self.args.paramsPath:
            if not os.path.exists(self.args.paramsPath):
                logging.error("Params path doesn't exist: %s" % self.args.paramsPath)
                sys.exit(1)

            # Use the SequencingChemistry data to select an ipd model
            if 'SequencingChemistry' in movieInfoTable.dtype.fields.keys():
                # Pick majority chemistry
                chemistries = movieInfoTable.SequencingChemistry.tolist()
                chemCounts = dict([(k, len(list(v))) for (k, v) in itertools.groupby(chemistries)])
                majorityChem = max(chemCounts, key=chemCounts.get)

                # Find the appropriate model file:
                ipdModel = os.path.join(self.args.paramsPath, majorityChem + ".h5")

                if majorityChem == 'unknown':
                    logging.warning("Chemistry is unknown. Falling back to built-in model")
                    ipdModel = None
                elif not os.path.exists(ipdModel):
                    logging.warning("Model not found: %s" % ipdModel)
                    logging.warning("Falling back to built-in model")
                    ipdModel = None
                else:
                    logging.info("Using Chemistry matched IPD model: %s" % ipdModel)

        self.ipdModel = IpdModel(contigs, ipdModel, self.args.modelIters)

    def _mainLoop(self):
        """
        Main loop
        First launch the worker and writer processes
        Then we loop over ReferenceGroups in the cmp.h5.  For each contig we will:
        1. Load the sequence into the main memory of the parent process
        3. Chunk up the contig and submit the chunk descriptions to the work queue
        Finally, wait for the writer process to finish.
        """

        # This looks scary but it's not.  Python uses reference
        # counting and has a secondary, optional garbage collector for
        # collecting garbage cycles.  Unfortunately when a cyclic GC
        # happens when a thread is calling cPickle.dumps, the
        # interpreter crashes sometimes.  See Bug 19704.  Since we
        # don't leak garbage cycles, disabling the cyclic GC is
        # essentially harmless.
        gc.disable()

        # Load reference and IpdModel
        self.loadReferenceAndModel(self.args.reference, self.args.infile)

        # Spawn workers
        self._launchSlaveProcesses()

        # WARNING -- cmp.h5 file must be opened AFTER worker processes have been spawned
        # cmp.h5 we're using -- use this to orchestrate the work
        self.cmph5 = CmpH5Reader(self.args.infile)
        logging.info('Generating kinetics summary for [%s]' % self.args.infile)

        #self.referenceMap = self.cmph5['/RefGroup'].asDict('RefInfoID', 'ID')
        #self.alnInfo = self.cmph5['/AlnInfo'].asRecArray()

        # Main loop -- we loop over ReferenceGroups in the cmp.h5.  For each contig we will:
        # 1. Load the sequence into the main memory of the parent process
        # 2. Fork the workers
        # 3. chunk up the contig and

        self.workChunkCounter = 0

        if self.options.refId > -1:
            # Under the --refId option, rather than iterating over references, process
            # just the one specified reference.

            # ref = x[ self.options.refId ]
            ref = self.refInfo[self.options.refId]
            logging.info('Processing reference entry: [%s]' % ref.Name)
            self._queueChunksForReference(ref)

        else:
            # Iterate over references

            for ref in self.refInfo:
                logging.info('Processing reference entry: [%s]' % ref.Name)
                self._queueChunksForReference(ref)

        # Shutdown worker threads with None sentinels
        for i in xrange(self.args.numWorkers):
            self._workQueue.put(None)

        for w in self._workers:
            w.join()

        # Join on the result queue and the resultsCollector process.
        # This ensures all the results are written before shutdown.
        self.monitoringThread.join()
        self._resultsQueue.join()
        self._resultCollectorProcess.join()
        logging.info("ipdSummary.py finished. Exiting.")
        del self.cmph5
        return 0


def monitorChildProcesses(children):
    """
    Monitors child processes: promptly exits if a child is found to
    have exited with a nonzero exit code received; otherwise returns
    zero when all processes exit cleanly (0).

    This approach is portable--catching SIGCHLD doesn't work on
    Windows.
    """
    while True:
        all_exited = all(not p.is_alive() for p in children)
        nonzero_exits = [p.exitcode for p in children if p.exitcode]
        if nonzero_exits:
            exitcode = nonzero_exits[0]
            logging.error("Child process exited with exitcode=%d.  Aborting." % exitcode)

            # Kill all the child processes
            for p in children:
                if p.is_alive():
                    p.terminate()

            os._exit(exitcode)
        elif all_exited:
            return 0
        time.sleep(1)


def main():
    try:
        kt = KineticsToolsRunner()
        kt.parseArgs()
        rcode = kt.start()
        return rcode
    except Exception as e:
        type, value, tb = sys.exc_info()
        traceback.print_exc(file=sys.stderr)
        # Note: if kt.args.usePdb
        # This won't work. If an exception is raised in parseArgs,
        # then kt.args is not defined yet.
        if '--pdb' in sys.argv:
            try:
                # this has better integration with ipython and is nicer
                # pip install ipdb
                import ipdb
                ipdb.post_mortem(tb)
            except ImportError:
                import pdb
                pdb.post_mortem(tb)
        else:
            # exit non-zero
            raise


if __name__ == "__main__":
    sys.exit(main())
