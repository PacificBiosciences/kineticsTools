"""
Gather task for BigWig output files
"""

from collections import namedtuple, OrderedDict, defaultdict
import itertools
import logging
import math
import os.path as op
import sys

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)
__version__ = "0.2"


def gather_bigwig(input_files, output_file):
    import pyBigWig
    chr_lengths = {}
    FileInfo = namedtuple("FileInfo", ("file_name", "file_id", "seqids"))
    files_info = []
    for i, file_name in enumerate(input_files):
        log.info("Reading header info from {f}...".format(f=file_name))
        if op.getsize(file_name) == 0:
            continue
        try:
            file_id = int(op.dirname(file_name).split("-")[-1])
        except ValueError:
            file_id = i
        bw_chunk = pyBigWig.open(file_name)
        seqids = []
        for (seqid, length) in bw_chunk.chroms().items():
            chr_lengths.setdefault(seqid, 0)
            chr_lengths[seqid] = max(length, chr_lengths[seqid])
            seqids.append(seqid)
        files_info.append(FileInfo(file_name, file_id, seqids))
        bw_chunk.close()
    if len(files_info) == 0:
        with open(output_file, "wb") as f:
            return output_file
    bw = pyBigWig.open(output_file, "w")
    files_info.sort(key=lambda f: f.file_id)
    regions = OrderedDict()
    seqid_files = defaultdict(list)
    for f in files_info:
        for seqid in f.seqids:
            log.debug("{f} ({i}): {s} {l}".format(f=f.file_name,
                                                  i=f.file_id,
                                                  s=seqid,
                                                  l=chr_lengths[seqid]))
            regions[seqid] = chr_lengths[seqid]
            seqid_files[seqid].append(f)
    bw.addHeader([(k, v) for k, v in regions.items()])
    seq_chunk = namedtuple("SeqChunk", ("starts", "ends", "values"))
    for (seqid, length) in regions.items():
        log.info("Collecting values for {i}...".format(i=seqid))
        chunks = []
        k = 0
        for file_info in seqid_files[seqid]:
            log.info("Reading values from {f}".format(f=file_info.file_name))
            bw_chunk = pyBigWig.open(file_info.file_name)
            starts, ends, values = [], [], []
            chr_max = bw_chunk.chroms()[seqid]
            for i, val in enumerate(bw_chunk.values(seqid, 0, chr_max)):
                if not math.isnan(val):
                    starts.append(i)
                    ends.append(i + 1)
                    values.append(val)
                    k += 1
            chunks.append(seq_chunk(starts, ends, values))
            bw_chunk.close()
        chunks.sort(key=lambda c: c.starts[0])
        starts = list(itertools.chain(*[x.starts for x in chunks]))
        ends = list(itertools.chain(*[x.ends for x in chunks]))
        values = list(itertools.chain(*[x.values for x in chunks]))
        seqids = [seqid] * len(starts)
        log.info("Adding {i}:{s}-{e}".format(i=seqid, s=starts[0], e=ends[-1]))
        bw.addEntries(seqids, starts, ends=ends, values=values)
    bw.close()
    return 0


def _get_parser():
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="INFO")
    p.add_argument("merged", help="Name of merged BigWig file")
    p.add_argument("chunks", nargs="+", help="Chunked BigWig files")
    return p


def run_args(args):
    return gather_bigwig(args.chunks, args.merged)


def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=run_args,
        alog=log,
        setup_log_func=setup_log)


if __name__ == '__main__':
    sys.exit(main())
