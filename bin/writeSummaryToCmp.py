#!/usr/bin/env python

import cProfile
from pbcore.io import GffReader, Gff3Record
import os
import logging
import sys

from pbcore.util.ToolRunner import PBToolRunner

__version__ = "1.0"


class IpdRatioSummaryWriter(PBToolRunner):

    def __init__(self):
        desc = ['Summarizes kinetic modifications in the alignment_summary.gff file',
                'Notes: For all command-line arguments, default values are listed in [].']
        super(IpdRatioSummaryWriter, self).__init__('\n'.join(desc))

        self.parser.add_argument('--pickle',
                                 dest="pickle",
                                 help='Name of input GFF file [%(default)s]')

        self.parser.add_argument('--alignmentSummary',
                                 dest="alignmentSummary",
                                 help='Name alignment summary file [%(default)s]')

        self.parser.add_argument('--outfile',
                                 dest="outfile",
                                 help='Name of modified alignment summary GFF file [%(default)s]')

        self.parser.add_argument("--profile",
                                 action="store_true",
                                 dest="doProfiling",
                                 default=False,
                                 help="Enable Python-level profiling (using cProfile).")

    def getVersion(self):
        return __version__

    def validateArgs(self):
        if not os.path.exists(self.args.modifications):
            self.parser.error('input modifications gff file provided does not exist')

        if not os.path.exists(self.args.alignmentSummary):
            self.parser.error('input alignment summary gff file provided does not exist')

    def run(self):
        self.options = self.args

        # Log generously
        logFormat = '%(asctime)s [%(levelname)s] %(message)s'
        logging.basicConfig(level=logging.INFO, format=logFormat)
        stdOutHandler = logging.StreamHandler(sys.stdout)
        logging.Logger.root.addHandler(stdOutHandler)
        logging.info("t1")

        if self.args.doProfiling:
            cProfile.runctx("self._mainLoop()",
                            globals=globals(),
                            locals=locals(),
                            filename="profile-main4.out")

        else:
            return self._mainLoop()

    def _mainLoop(self):

        # Read in the existing modifications.gff
        modReader = GffReader(self.args.modifications)

        # Set up some additional headers to be injected
        headers = [
            ('source', 'kineticModificationCaller 1.3.1'),
            ('source-commandline', " ".join(sys.argv)),
            ('attribute-description', 'modsfwd - count of detected DNA modifications on forward strand'),
            ('attribute-description', 'modsrev - count of detected DNA modifications on reverse strand')
        ]

        # Get modification calls
        hits = [{"pos": x.start, "strand": x.strand} for x in modReader if x.type == 'modified_base']

        # Summary reader
        summaryFile = file(self.args.alignmentSummary)

        # Modified gff file
        summaryWriter = file(self.args.outfile, "w")

        self.seqMap = {}
        inHeader = True

        # Loop through
        for line in summaryFile:
            # Pass any metadata line straight through
            if line[0] == "#":

                # Parse headers
                splitFields = line.replace('#', '').split(' ')
                field = splitFields[0]
                value = " ".join(splitFields[1:])
                if field == 'sequence-header':
                    [internalTag, delim, externalTag] = value.strip().partition(' ')
                    self.seqMap[internalTag] = externalTag
                print(line.strip(), file=summaryWriter)
                continue

            if inHeader:
                # We are at the end of the header -- write the tool-specific headers
                for field in headers:
                    print(("##%s %s" % field), file=summaryWriter)
                inHeader = False

            # Parse the line
            rec = Gff3Record.fromString(line)

            if rec.type == 'region':
                # Get the hits in this interval, add them to the gff record
                intervalHits = [h for h in hits if rec.start <= h['pos'] <= rec.end]
                strand0Hits = len([h for h in intervalHits if h['strand'] == '+'])
                strand1Hits = len([h for h in intervalHits if h['strand'] == '-'])

                rec.modsfwd = strand0Hits
                rec.modsrev = strand1Hits

                print(str(rec), file=summaryWriter)

if __name__ == "__main__":
    kt = ModificationSummary()
    kt.start()
