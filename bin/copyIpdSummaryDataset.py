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
