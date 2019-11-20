#!/usr/bin/env python

import os, sys

from pbcore.util.ToolRunner import PBToolRunner
from operator import xor
import h5py

# Version info
__version__ = "1.1"


def validateFile(p):
    if os.path.isfile(p):
        return os.path.abspath(p)
    else:
        raise IOError("Unable to find {p}.".format(p=p))


class CopyIpdSummaryDatasets(PBToolRunner):
    """
    Copy IpdRatio datasets from infile:/ref0000x to outfile:/ref000x/Kinetics/IpdRatio
     or                    from infile:/ref0000x to mergefile:/ref000x (Dataset)
    """

    def __init__(self):
        super(CopyIpdSummaryDatasets, self).__init__(CopyIpdSummaryDatasets.__doc__)

        self.parser.add_argument('--infile',
                                 required=True,
                                 type=validateFile,
                                 dest='infile',
                                 help='Input cmp.h5 filename')

        self.parser.add_argument('--outfile',
                                 type=validateFile,
                                 required=False,
                                 help='Output cmp.h5 filename')

        self.parser.add_argument('--mergefile',
                                 type=validateFile,
                                 required=False,
                                 help='Filename of output h5 file for merging')

    def validateArgs(self):
        if not xor(bool(self.args.mergefile), bool(self.args.outfile)):
            raise Exception("Exactly one of --outfile, --mergefile is required")

    def getVersion(self):
        return __version__

    def copyToCmpH5(self):
        inFile = h5py.File(self.args.infile, mode='r')
        outFile = h5py.File(self.args.outfile, mode='r+')

        for refDataset in inFile.items():
            (name, ds) = refDataset

            if '/' + name in outFile:
                targetGroup = outFile['/' + name]

                if 'Kinetics' in targetGroup:
                    kinGroup = targetGroup['Kinetics']
                else:
                    kinGroup = targetGroup.create_group('Kinetics')

                if 'IpdRatio' in kinGroup:
                    del kinGroup['IpdRatio']

                h5py.h5o.copy(inFile.id, name, kinGroup.id, 'IpdRatio')


    def copyToMergeFile(self):
        inFile = h5py.File(self.args.infile, mode='r')
        mergeFile = h5py.File(self.args.mergefile, mode='r+')

        for refDataset in inFile.items():
            (name, ds) = refDataset
            h5py.h5o.copy(inFile.id, name, mergeFile.id, name)



    def run(self):
        if self.args.outfile: self.copyToCmpH5()
        else:                 self.copyToMergeFile()
        return 0


def main():
    kt = CopyIpdSummaryDatasets()
    rcode = kt.start()
    return rcode


if __name__ == "__main__":
    sys.exit(main())
