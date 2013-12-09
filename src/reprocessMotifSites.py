#!/usr/bin/env python

#
# Purpose:
#
#    Estimates modified fraction at all undetected sites within each identified motif.
#
#    This script calls ./pbtools/kineticsTools/kineticForReprocessing.py instead of KineticWorker.py
#    Otherwise similar to ipdSummary.py, and uses WorkerProcess, ResultWriter
#
# Required Inputs:
#
#    motifs GFF file
#    motif_summary GFF file
#    reference file
#    input cmpH5 file
#    name of output GFF file
#    either --undetected flag or modifications GFF file
#       NOTE:  assume that modifications GFF file has methylFracEst fields filled in.
#
# Outputs:
#
#    a GFF file containing modified fraction estimates:
#      - includes all undetected sites for each identified motif in motifs GFF
#      - if modifications GFF was provided, then also includes detected sites in each motif
#
#
# Sample Commands:
#
# cd /mnt/secondary/Smrtanalysis/opt/smrtanalysis-1.4.0-ubuntu-117447/common/
#
# reprocessMotifSites.py  --reference ./references/SP48_HGAP_v7_Assembly_3contigs/ --motifs ./jobs/055/055789/data/motifs.gff.gz \
#       --motif_summary ./jobs/055/055789/data/motif_summary.csv --undetected --gff ~/retrain/SP48_undetected.gff ./jobs/055/055789/data/aligned_reads.cmp.h5
#
# reprocessMotifSites.py  --reference ./references/SP48_HGAP_v7_Assembly_3contigs/ --motifs ./jobs/055/055789/data/motifs.gff.gz \
#       --motif_summary ./jobs/055/055789/data/motif_summary.csv --modifications ./jobs/055/055789/data/modifications.gff.gz \
#       --gff ~/retrain/SP48_undetected.gff ./jobs/055/055789/data/aligned_reads.cmp.h5
#


import cProfile
import gc

from pbcore.deprecated import ReferenceEntry

import os
import logging
import sys
import multiprocessing
import time
import threading
import numpy as np
import Queue

# from pbcore.io.ReferenceEntry import ReferenceEntry
from pbcore.io import CmpH5Reader
from pbcore.util.ToolRunner import PBToolRunner

# Replace:
# from pbtools.kineticsTools.KineticWorker import KineticWorker, KineticWorkerThread, KineticWorkerProcess
from pbtools.kineticsTools.kineticForReprocessing import KineticReprocessWorker, KineticWorkerThread, KineticWorkerProcess

from pbtools.kineticsTools.ResultWriter import KineticsWriter
from pbtools.kineticsTools.ipdModel import IpdModel

# New:
from pbcore.io import GffReader
import csv


# Addition to ipdSummary.py since reprocessMotifSites.py was first written:
from pbtools.kineticsTools.ReferenceUtils import ReferenceUtils


# Version info
__p4revision__ = "$Revision: #1 $"
__p4change__ = "$Change: 100972 $"
revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum = int(__p4change__.strip("$").split(":")[-1])
__version__ = "1.4"


class ReprocessMotifSites(PBToolRunner):

    def __init__(self):
        desc = ['For all sites in motifs.gff, reports estimated methylated fraction and 95% confidence interval',
                'Notes: For all command-line arguments, default values are listed in [].']
        super(ReprocessMotifSites, self).__init__('\n'.join(desc))

        self.parser.add_argument('--numWorkers',
                                 dest='numWorkers',
                                 default=-1,  # Defaults to using all logical CPUs
                                 type=int,
                                 help='Number of thread to use (-1 uses all logical cpus)')

        self.parser.add_argument('infile',
                                 metavar='input.cmp.h5',
                                 help='Input cmp.h5 filename')

        # self.parser.add_argument('--control',
        #     dest='control',
        #     default=None,
        #     help='cmph.h5 file ')

        # self.parser.add_argument('--outfile',
        #     dest='outfile',
        #     default=None,
        #    help='Use this option to generate all possible output files. Argument here is the root filename of the output files.')

        self.parser.add_argument('--gff',
                                 dest='gff',
                                 default=None,
                                 help='Name of output GFF file [%(default)s]')

        # self.parser.add_argument('--identify',
        #     dest='identify',
        #     default=False,
        #     help='Identify modification types. Comma-separated list of know modification types. Current options are: m6A, m4C, m5C_TET. Cannot be used with --control')

        # self.parser.add_argument('--csv',
        #     dest='csv',
        #     default=None,
        #     help='Name of output CSV file [%(default)s]')

        # self.parser.add_argument('--pickle',
        #     dest='pickle',
        #     default=None,
        #     help='Name of output pickle file [%(default)s]')

        # self.parser.add_argument('--summary_h5',
        #     dest='summary_h5',
        #     default=None,
        #     help='Name of output summary h5 file [%(default)s]')

        self.parser.add_argument('--reference',
                                 dest='reference',
                                 required=True,
                                 help='Path to reference FASTA file')

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
                                 default=1000,
                                 type=int,
                                 help='Max Queue Size')

        self.parser.add_argument('--maxCoverage',
                                 dest='maxCoverage',
                                 type=int, default=None,
                                 help='Maximum coverage to use at each site')

        self.parser.add_argument('--mapQvThreshold',
                                 dest='mapQvThreshold',
                                 type=float,
                                 default=-1.0)

        # self.parser.add_argument('--pvalue',
        #     dest='pvalue',
        #     default=0.01,
        #     type=float,
        #     help='p-value required to call a modified base')

        self.parser.add_argument('--subread_norm',
                                 dest='subread_norm',
                                 default=True,
                                 type=lambda x: x != 'False',
                                 help='Normalized subread ipds')

        self.parser.add_argument('--ipdModel',
                                 dest='ipdModel',
                                 default=None,
                                 help='Alternate synthetic IPD model HDF5 file')

        self.parser.add_argument('--cap_percentile',
                                 dest='cap_percentile',
                                 type=float,
                                 default=99.0,
                                 help='Global IPD percentile to cap IPDs at')

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

        # self.parser.add_argument("--methylFraction",
        #     action="store_true",
        #     dest="methylFraction",
        #     default=False,
        #     help="In the --identify mode, add --methylFraction to command line to estimate the methylated fraction, along with 95% confidence interval bounds.")

        # The following are in addition to ipdSummary.py's inputs:
        self.parser.add_argument('--motifs',
                                 dest="motifs",
                                 required=True,
                                 help='Name of motifs GFF file [%(default)s]')

        self.parser.add_argument('--motif_summary',
                                 dest="motif_summary",
                                 required=True,
                                 help='Name of motif summary CSV file')

        self.parser.add_argument('--undetected',
                                 action="store_true",
                                 dest="undetected",
                                 default=False,
                                 help="Setting this flag yields output with only undetected motif sites.")

        self.parser.add_argument('--modifications',
                                 dest="modifications",
                                 default=None,
                                 help='Name of modifications GFF file [%(default)s]')

        self.parser.add_argument('--oldData',
                                 action="store_true",
                                 dest="oldData",
                                 default=False,
                                 help="For datasets prior to 1.3.3 (use this option to increase testing possibilities)")

        # A new addition to ipdSummary.py

        self.parser.add_argument('--paramsPath',
                                 dest='paramsPath',
                                 default=None,
                                 help='Directory containing in-silico trained model for each chemistry')

        self.parser.add_argument('--modelIters',
                                 dest='modelIters',
                                 type=int,
                                 default=-1,
                                 help='[Internal] Number of GBM model iteration to use')

    def getVersion(self):
        return __version__

    def validateArgs(self):
        if not os.path.exists(self.args.infile):
            self.parser.error('input.cmp.h5 file provided does not exist')

        # Add checks corresponding to new required inputs:
        if not os.path.exists(self.args.motifs):
            self.parser.error('input motifs gff file provided does not exist')

        if not os.path.exists(self.args.motif_summary):
            self.parser.error('input motif_summary csv file provided does not exist')

        if not self.args.undetected and not os.path.exists(self.args.modifications):
            self.parser.error('either the --undetected flag must be set, or a valid modifications.gff must be provided')

    def run(self):

        # The following arguments are set in order to use ResultWriter.py as is:
        self.args.methylFraction = True
        self.args.outfile = None
        self.args.csv = None
        self.args.control = None
        self.args.summary_h5 = None
        self.args.pickle = None
        self.args.identify = False
        self.args.pvalue = 1.0

        self.options = self.args
        self.options.cmdLine = " ".join(sys.argv)
        self._workers = []

        # Log generously
        logFormat = '%(asctime)s [%(levelname)s] %(message)s'
        logging.basicConfig(level=logging.INFO, format=logFormat)
        stdOutHandler = logging.StreamHandler(sys.stdout)
        # logging.Logger.root.addHandler(stdOutHandler)
        # logging.info("t1")

        if self.args.doProfiling:
            cProfile.runctx("self._mainLoop()",
                            globals=globals(),
                            locals=locals(),
                            filename="profile-main4.out")

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

    def _queueChunksForReference(self):

        # Read in motif_summary.csv
        motifInfo = {}
        reader = csv.reader(open(self.args.motif_summary, 'r'), delimiter=',')
        reader.next()
        if self.options.oldData:
            col = 1
        else:
            col = 2
        for row in reader:
            motifInfo[row[0]] = row[col]

        # Figure out the length of the motifs file:
        motReader = GffReader(self.args.motifs)
        if self.options.undetected:
            motifDicts = [{"seqID": x.seqid, "type": x.type, "score": x.score, "pos": x.start, "strand": x.strand, "attributes": x.attributes}
                          for x in motReader if x.type == '.']
        else:
            motifDicts = [{"seqID": x.seqid, "type": x.type, "score": x.score, "pos": x.start, "strand": x.strand, "attributes": x.attributes}
                          for x in motReader]

        refLength = len(motifDicts)

        # Maximum number of hits per chunk
        MAX_HITS = 500
        nBases = min(refLength, self.args.maxLength)
        nBlocks = max(self.options.numWorkers * 4, nBases / MAX_HITS)

        # Block layout
        blockSize = min(nBases, max(nBases / nBlocks + 1, 100))
        blockStarts = np.arange(0, nBases, step=blockSize)
        blockEnds = blockStarts + blockSize
        blocks = zip(blockStarts, blockEnds)

        if self.options.undetected:
            self.options.modifications = None

        # Queue up work blocks
        for block in blocks:
            # NOTE! The format of a work chunk is (refId <int>, refStartBase <int>, refEndBase <int>)
            # chunk = (refInfoId, block[0], block[1])
            # chunk = (self.options.motifs, self.refInfo, motifInfo, self.options.modifications, self.options.undetected, self.options.oldData, block[0], block[1])
            chunk = (motifDicts[block[0]:block[1]], self.refInfo, motifInfo, self.options.modifications, self.options.undetected, self.options.oldData, block[0], block[1])
            self._workQueue.put((self.workChunkCounter, chunk))
            self.workChunkCounter += 1

    def loadReference(self):
        # FIXME - support a bare fasta file as well?
        self.referenceEntry = ReferenceEntry(self.args.reference)
        self.refInfo = self.referenceEntry.contigs

        if self.args.ipdModel:
            self.lutPath = self.args.ipdModel
            if not os.path.exists(self.lutPath):
                logging.info("Couldn't find model file: %s" % self.lutPath)
                raise Exception("Couldn't find model file: %s" % self.lutPath)
        else:
            self.lutPath = None

        self.ipdModel = IpdModel(self.referenceEntry, self.lutPath)

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

        # See comments in ipdSummary.py
        gc.disable()

        # Load reference and IpdModel
        # self.loadReference()

        # Load reference and IpdModel
        self.loadReferenceAndModel(self.args.reference, self.args.infile)

        # Spawn workers
        self._launchSlaveProcesses()

        # cmp.h5 we're using -- use this to orchestrate the work
        self.cmph5 = CmpH5Reader(self.args.infile)
        logging.info('Generating kinetics summary for [%s]' % self.args.infile)

        self.workChunkCounter = 0
        self._queueChunksForReference()

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
        logging.info("reprocessMotifSites.py finished. Exiting.")
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


if __name__ == "__main__":
    kt = ReprocessMotifSites()
    kt.start()
