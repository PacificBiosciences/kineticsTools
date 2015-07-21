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
"""
Tool for detecting DNA base-modifications from kinetic signatures
Notes: For all command-line arguments, default values are listed in []
"""

import cProfile
import functools
import gc
import itertools
import argparse
import json

import os
import logging
import sys
import multiprocessing
import time
import threading
import numpy as np
import Queue
import traceback
from pkg_resources import Requirement, resource_filename

from pbcommand.cli import (pacbio_args_or_contract_runner,
                           get_default_argparser)
from pbcommand.models import TaskTypes, FileTypes, get_default_contract_parser
from pbcommand.utils import setup_log
from pbcommand.common_options import add_resolved_tool_contract_option

from pbcore.io import AlignmentSet
from kineticsTools.KineticWorker import KineticWorkerThread, KineticWorkerProcess
from kineticsTools.ResultWriter import KineticsWriter
from kineticsTools.ipdModel import IpdModel
from kineticsTools.ReferenceUtils import ReferenceUtils

# Version info
__p4revision__ = "$Revision: #1 $"
__p4change__ = "$Change: 100972 $"
revNum = int(__p4revision__.strip("$").split(" ")[1].strip("#"))
changeNum = int(__p4change__.strip("$").split(":")[-1])
__version__ = "2.2"

def _getResourcePath():
    return resource_filename(Requirement.parse('kineticsTools'),'kineticsTools/resources')

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

# XXX FUTURE: use this!  then we can generate the tool contract dynamically
# instead of editing a static file.
def get_contract_parser():
    nproc = 1
    resources = ()
    driver_exe = "ipdSummary.py --resolved-tool-contract "
    p = get_default_contract_parser(
        "kineticsTools.ipdSummary",
        __version__,
        __doc__,
        driver_exe,
        TaskTypes.DISTRIBUTED,
        nproc,
        resources)
    p.add_input_file_type(FileTypes.DS_BAM, "infile",
        "Alignment DataSet", "BAM or Alignment DataSet")
    p.add_input_file_type(FileTypes.DS_REF, "reference",
        "Reference DataSet", "Fasta or Reference DataSet")
    p.add_output_file_type(FileTypes.GFF, "gff",
        name="GFF file",
        description="GFF file of modified bases",
        default_name="basemods.gff")
    p.add_output_file_type(FileTypes.CSV, "csv",
        name="CSV file",
        description="CSV file of per-nucleotide information",
        default_name="basemods.csv")
    p.add_int("numWorkers", "numWorkers", nproc, "Number of processors",
        "Number of processors")
    p.add_float("pvalue", "pvalue", 0.01, "P-value", "P-value cutoff")
    p.add_int("maxLength", "maxLength", int(3e12), "Max sequence length",
        "Maximum number of bases to process per contig")
    p.add_str(option_id="identify",
        option_str="identify",
        default=None,
        name="Identify basemods",
        description="Specific modifications to identify (comma-separated list")
    _get_more_options(p.arg_parser.parser)
    return p

def get_argument_parser():
    parser = get_default_argparser(__version__, __doc__)

    # Positional arguments:
    parser.add_argument('infile',
                        metavar='aligned.subreads.xml',
                        type=validateFile,
                        help='Input AlignmentSet filename')
    # Note: reference is actually not optional:
    parser.add_argument('--reference', '-r',
                        type=validateFile,
                        help='Path to reference FASTA file')
    parser.add_argument('--gff',
                        dest='gff',
                        default=None,
                        help='Name of output GFF file')
    parser.add_argument('--csv',
                        dest='csv',
                        default=None,
                        help='Name of output CSV file.')
    parser.add_argument('--identify',
        dest='identify',
        default=False,
        help='Identify modification types. Comma-separated list of know modification types. Current options are: m6A, m4C, m5C_TET. Cannot be used with --control')

    parser.add_argument("--methylFraction",
        action="store_true",
        dest="methylFraction",
        default=False,
        help="In the --identify mode, add --methylFraction to command line to estimate the methylated fraction, along with 95%% confidence interval bounds.")
    parser.add_argument("--maxLength",
                        default=3e12,
                        type=int,
                        help="Maximum number of bases to process per contig")
    parser.add_argument('--pvalue',
                        dest='pvalue',
                        default=0.01,
                        type=float,
                        help='p-value required to call a modified base')
    parser.add_argument('--numWorkers', '-j',
        dest='numWorkers',
        default=1,
        type=int,
        help='Number of thread to use (-1 uses all logical cpus)')
    _get_more_options(parser)
    add_resolved_tool_contract_option(parser)
    # FIXME temporary workaround for parser chaos
    class EmitToolContractAction(argparse.Action):
        def __call__(self, parser_, namespace, values, option_string=None):
            parser2 = get_contract_parser()
            sys.stdout.write(json.dumps(parser2.to_contract(), indent=4)+'\n')
            sys.exit(0)
    parser.add_argument("--emit-tool-contract",
                        nargs=0,
                        action=EmitToolContractAction)
    return parser

def _get_more_options(parser):
    """
    Advanced options that won't be exposed via tool contract interface.
    """
    parser.add_argument('--outfile',
        dest='outfile',
        default=None,
        help='Use this option to generate all possible output files. Argument here is the root filename of the output files.')

    # FIXME: Need to add an extra check for this; it can only be used if --useLDA flag is set.
    parser.add_argument('--m5Cgff',
        dest='m5Cgff',
        default=None,
        help='Name of output GFF file containing m5C scores')

    # FIXME: Make sure that this is specified if --useLDA flag is set.
    parser.add_argument('--m5Cclassifer',
                             dest='m5Cclassifier',
                             default=None,
                             help='Specify csv file containing a 127 x 2 matrix')


    parser.add_argument('--csv_h5',
                             dest='csv_h5',
                             default=None,
                             help='Name of csv output to be written in hdf5 format.')

    parser.add_argument('--pickle',
                             dest='pickle',
                             default=None,
                             help='Name of output pickle file.')

    parser.add_argument('--summary_h5',
                             dest='summary_h5',
                             default=None,
                             help='Name of output summary h5 file.')


    parser.add_argument('--ms_csv',
                             dest='ms_csv',
                             default=None,
                             help='Multisite detection CSV file.')


    # Calculation options:


    parser.add_argument('--control',
                             dest='control',
                             default=None,
                             type=validateNoneOrFile,
                             help='cmph.h5 file containing a control sample. Tool will perform a case-control analysis')

    # Temporary addition to test LDA for Ca5C detection:
    parser.add_argument('--useLDA',
                             action="store_true",
                             dest='useLDA',
                             default=False,
                             help='Set this flag to debug LDA for m5C/Ca5C detection')



    # Parameter options:

    parser.add_argument('--paramsPath',
                             dest='paramsPath',
                             default=_getResourcePath(),
                             type=validateNoneOrDir,
                             help='Directory containing in-silico trained model for each chemistry')

    parser.add_argument('--minCoverage',
                             dest='minCoverage',
                             default=3,
                             type=int,
                             help='Minimum coverage required to call a modified base')

    parser.add_argument('--maxQueueSize',
                             dest='maxQueueSize',
                             default=20,
                             type=int,
                             help='Max Queue Size')

    parser.add_argument('--maxCoverage',
                             dest='maxCoverage',
                             type=int, default=-1,
                             help='Maximum coverage to use at each site')

    parser.add_argument('--mapQvThreshold',
                             dest='mapQvThreshold',
                             type=float,
                             default=-1.0)

    parser.add_argument('--ipdModel',
                             dest='ipdModel',
                             default=None,
                             help='Alternate synthetic IPD model HDF5 file')

    parser.add_argument('--modelIters',
                             dest='modelIters',
                             type=int,
                             default=-1,
                             help='[Internal] Number of GBM model iteration to use')

    parser.add_argument('--cap_percentile',
                             dest='cap_percentile',
                             type=float,
                             default=99.0,
                             help='Global IPD percentile to cap IPDs at')


    parser.add_argument("--methylMinCov",
                             type=int,
                             dest='methylMinCov',
                             default=10,
                             help="Do not try to estimate methylFraction unless coverage is at least this.")

    parser.add_argument("--identifyMinCov",
                             type=int,
                             dest='identifyMinCov',
                             default=5,
                             help="Do not try to identify the modification type unless coverage is at least this.")

    parser.add_argument("--maxAlignments",
                             type=int,
                             dest="maxAlignments",
                             default=1500,
                             help="Maximum number of alignments to use for a given window")


    # Computation management options:

    parser.add_argument("-w", "--referenceWindow", "--referenceWindows",
                             "--refContigs", # backwards compatibility
                             type=str,
                             dest='referenceWindowsAsString',
                             default=None,
                             help="The window (or multiple comma-delimited windows) of the reference to " + \
                                  "be processed, in the format refGroup[:refStart-refEnd] "               + \
                                  "(default: entire reference).")

    def slurpWindowFile(fname):
        return ",".join(map(str.strip, open(fname).readlines()))


    parser.add_argument("--refContigIndex", type=int, dest='refContigIndex', default=-1,
                             help="For debugging purposes only - rather than enter a reference contig name, simply enter an index" ) 

    parser.add_argument("-W", "--referenceWindowsFile",
                              "--refContigsFile", # backwards compatibility
                             type=slurpWindowFile,
                             dest='referenceWindowsAsString',
                             default=None,
                             help="A file containing reference window designations, one per line")

    parser.add_argument("--skipUnrecognizedContigs",
                             type=bool,
                             default=False,
                             help="Whether to skip, or abort, unrecognized contigs in the -w/-W flags")
    
    # Debugging help options:

    parser.add_argument("--threaded", "-T",
                             action="store_true",
                             dest="threaded",
                             default=False,
                             help="Run threads instead of processes (for debugging purposes only)")

    parser.add_argument("--profile",
                             action="store_true",
                             dest="doProfiling",
                             default=False,
                             help="Enable Python-level profiling (using cProfile).")

    parser.add_argument('--usePdb',
                             action='store_true',
                             dest="usePdb",
                             default=False,
                             help="Enable dropping down into pdb debugger if an Exception is raised.")

    parser.add_argument("--seed",
                             action="store",
                             dest="randomSeed",
                             type=int,
                             default=None,
                             help="Random seed (for development and debugging purposes only)")

    # Verbosity
    parser.add_argument("--verbose",
                        action="store_true",
                        default=False)
    return parser


class KineticsToolsRunner(object):
    def __init__(self, args):
        self.args = args
        self._sharedAlignmentSet = None

    def start(self):
        self.validateArgs()
        return self.run()

    def getVersion(self):
        return __version__

    def validateArgs(self):
        if not os.path.exists(self.args.infile):
            self.parser.error('input.cmp.h5 file provided does not exist')

        if self.args.identify and self.args.control:
            self.parser.error('--control and --identify are mutally exclusive. Please choose one or the other')

        if self.args.useLDA:
            if self.args.m5Cclassifier is None:
                self.parser.error('Please specify a folder containing forward.csv and reverse.csv classifiers in --m5Cclassifier.')

        if self.args.m5Cgff:
            if not self.args.useLDA:
                self.parser.error('m5Cgff file can only be generated in --useLDA mode.')

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

        # set random seed
        # XXX note that this is *not* guaranteed to yield reproducible results
        # indepenently of the number of processing cores used!
        if self.options.randomSeed is not None:
            np.random.seed(self.options.randomSeed)

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
            p = WorkerType(self.options, self._workQueue, self._resultsQueue,
                self.ipdModel,
                sharedAlignmentSet=self._sharedAlignmentSet)
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

    def _queueChunksForWindow(self, refWindow):
        """
        Compute the chunk extents and queue up the work for a single reference
        """
        winId = refWindow.refId
        winStart = refWindow.start
        winEnd = refWindow.end
        pass

    def loadReferenceAndModel(self, referencePath):
        assert self._sharedAlignmentSet is not None
        # Load the reference contigs - annotated with their refID from the cmp.h5
        logging.info("Loading reference contigs %s" % referencePath)
        contigs = ReferenceUtils.loadReferenceContigs(referencePath,
            alignmentSet=self._sharedAlignmentSet)

        # There are three different ways the ipdModel can be loaded.
        # In order of precedence they are:
        # 1. Explicit path passed to --ipdModel
        # 2. Path to parameter bundle, model selected using the cmp.h5's sequencingChemistry tags
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

            majorityChem = ReferenceUtils.loadAlignmentChemistry(
                self._sharedAlignmentSet)
            ipdModel = os.path.join(self.args.paramsPath, majorityChem + ".h5")
            if majorityChem == 'unknown':
                logging.error("Chemistry cannot be identified---cannot perform kinetic analysis")
                sys.exit(1)
            elif not os.path.exists(ipdModel):
                logging.error("Aborting, no kinetics model available for this chemistry: %s" % ipdModel)
                sys.exit(1)
            else:
                logging.info("Using Chemistry matched IPD model: %s" % ipdModel)

        self.ipdModel = IpdModel(contigs, ipdModel, self.args.modelIters)

    def loadSharedAlignmentSet(self, cmpH5Filename):
        """
        Read the input AlignmentSet so the indices can be shared with the
        slaves.  This is also used to pass to ReferenceUtils for setting up
        the ipdModel object.
        """
        logging.info("Reading AlignmentSet: %s" % cmpH5Filename)
        logging.info("           reference: %s" % self.args.reference)
        self._sharedAlignmentSet = AlignmentSet(cmpH5Filename,
            referenceFastaFname=self.args.reference)
        # XXX this should ensure that the file(s) get opened, including any
        # .pbi indices - but need to confirm this
        self.refInfo = self._sharedAlignmentSet.referenceInfoTable

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

        # Load a copy of the cmpH5 alignment index to share with the slaves
        self.loadSharedAlignmentSet(self.args.infile)

        # Load reference and IpdModel
        self.loadReferenceAndModel(self.args.reference)
        
        # Spawn workers
        self._launchSlaveProcesses()

        # WARNING -- cmp.h5 file must be opened AFTER worker processes have been spawned
        # cmp.h5 we're using -- use this to orchestrate the work
        self.cmph5 = self._sharedAlignmentSet
        logging.info('Generating kinetics summary for [%s]' % self.args.infile)

        #self.referenceMap = self.cmph5['/RefGroup'].asDict('RefInfoID', 'ID')
        #self.alnInfo = self.cmph5['/AlnInfo'].asRecArray()

        # Resolve the windows that will be visited.
        if self.args.referenceWindowsAsString is not None:
            self.referenceWindows = []
            for s in self.args.referenceWindowsAsString.split(","):
                try:
                    win = ReferenceUtils.parseReferenceWindow(s, self.cmph5.referenceInfo)
                    self.referenceWindows.append(win)
                except:
                    if self.args.skipUnrecognizedContigs:
                        continue
                    else:
                        raise Exception, "Unrecognized contig!"
        else:
            self.referenceWindows = ReferenceUtils.createReferenceWindows(
                self.refInfo)

        # Main loop -- we loop over ReferenceGroups in the cmp.h5.  For each contig we will:
        # 1. Load the sequence into the main memory of the parent process
        # 2. Fork the workers
        # 3. chunk up the contig and

        self.workChunkCounter = 0

        # Iterate over references
        for window in self.referenceWindows:
            logging.info('Processing window/contig: %s' % (window,))
            for chunk in ReferenceUtils.enumerateChunks(1000, window):
                self._workQueue.put((self.workChunkCounter, chunk))
                self.workChunkCounter += 1

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

def args_runner(args):
    log = logging.getLogger()
    if args.verbose:
        log.setLevel(logging.INFO)
    else:
        log.setLevel(logging.WARN)
    kt = KineticsToolsRunner(args)
    return kt.start()

def resolved_tool_contract_runner(resolved_contract):
    """
    Run ipdSummary from a resolved tool contract.  This basically just
    translates the contract into arguments that can be passed to the argparse
    parser and then args_runner.

    :param resolved_contract:
    :type resolved_contract: ResolvedToolContract
    :return: Exit code
    """
    alignment_path = resolved_contract.task.input_files[0]
    reference_path = resolved_contract.task.input_files[1]
    gff_path = resolved_contract.task.output_files[0]
    csv_path = resolved_contract.task.output_files[1]
    args = [
        alignment_path,
        "--reference", reference_path,
        "--gff", gff_path,
        "--csv", csv_path,
        "--numWorkers", str(resolved_contract.task.nproc),
        "--pvalue", str(resolved_contract.task.options["basemods.pvalue"]),
    ]
    if resolved_contract.task.options["basemods.max_length"]:
        args.extend([
            "--maxLength",
            str(resolved_contract.task.options["basemods.max_length"]),
        ])
    if resolved_contract.task.options["basemods.compute_methyl_fraction"]:
        args.append("--methylFraction")
    if resolved_contract.task.options["basemods.identify"]:
        args.extend([
            "--identify",
            resolved_contract.task.options["basemods.identify"],
        ])
    parser = get_argument_parser()
    args_ = parser.parse_args(args)
    return args_runner(args_)

def main(argv=sys.argv, out=sys.stdout):
    # Log generously
    logFormat = '%(asctime)s [%(levelname)s] %(message)s'
    logging.basicConfig(format=logFormat, level=logging.WARN)
    stdOutHandler = logging.StreamHandler(sys.stdout)
    log = logging.getLogger()
    try:
        mp = get_argument_parser()
        return pacbio_args_or_contract_runner(argv[1:],
                                              mp,
                                              args_runner,
                                              resolved_tool_contract_runner,
                                              log,
                                              lambda *args: log)
    # FIXME is there a more central place to deal with this?
    except Exception as e:
        type, value, tb = sys.exc_info()
        traceback.print_exc(file=sys.stderr)
        # Note: if kt.args.usePdb
        # This won't work. If an exception is raised in parseArgs,
        # then kt.args is not defined yet.
        if '--pdb' in argv:
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
