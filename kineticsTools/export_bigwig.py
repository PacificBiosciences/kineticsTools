
"""
Tool to export basemods annotations (from either CSV or GFF, depending on
reference size) to UCSC's BigWig binary format for display in SMRTView.
"""

from collections import namedtuple
import logging
import csv
import sys

from pbcore.io import ReferenceSet, GffReader
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "kineticstools.tasks.export_bigwig"
    DRIVER_EXE = "python -m kineticsTools.export_bigwig"
    MAX_CSV_ID = "kineticstools.task_options.genome_size_cutoff"
    MAX_CSV_DEFAULT = 10000000
    VERSION = "0.1"

BaseInfo = namedtuple("BaseInfo", ("seqid", "pos", "sense", "frames"))


def csv2bigwig(csv_file, bigwig_file):
    import pyBigWig
    with open(csv_file) as csv_in:
        reader = csv.reader(csv_in)
        ranges = {}
        records = []
        for rec in reader:
            if rec[0] == "refName":
                continue
            seqid = rec[0]
            if seqid[0] == seqid[-1] == '"':
                seqid = seqid[1:-1]
            ranges.setdefault(seqid, (sys.maxint, 0))
            pos = int(rec[1])
            ranges[seqid] = (min(ranges[seqid][0], pos),
                             max(ranges[seqid][0], pos))
            strand = int(rec[2])
            tMean = float(rec[5])
            frames = min(65535, int(tMean * 80))
            records.append(BaseInfo(seqid, pos, not bool(strand), frames))
        log.info("done reading CSV")
        records.sort(lambda a, b: cmp(a.pos, b.pos))
        records.sort(lambda a, b: cmp(a.seqid, b.seqid))
        bw = pyBigWig.open(bigwig_file, "w")
        regions = [(s, ranges[s][1]) for s in sorted(ranges.keys())]
        bw.addHeader(regions)
        k = 0
        seqids = []
        starts = []
        ends = []
        frames_enc = []
        while k < len(records):
            rec_plus = records[k] if records[k].sense else records[k + 1]
            rec_minus = records[k + 1] if records[k].sense else records[k]
            assert rec_plus.pos == rec_minus.pos, (rec_plus, rec_minus)
            seqids.append(rec_plus.seqid)
            starts.append(rec_plus.pos)
            ends.append(rec_plus.pos + 1)
            frames_enc.append(rec_minus.frames + 65536.0 * rec_plus.frames)
            k += 2
        log.info("Writing records for {n} bases".format(n=len(seqids)))
        bw.addEntries(seqids, starts, ends=ends, values=frames_enc)
        bw.close()
    return 0


def gff2bigwig(gff_file, bigwig_file):
    import pyBigWig
    with GffReader(gff_file) as gff:
        regions = []
        valid_seqids = set()
        for header in gff.headers:
            if header.startswith("##sequence-region"):
                _, seqid, start, end = header.split()
                regions.append((seqid, int(end) - int(start) + 1))
                valid_seqids.add(seqid)
        records = [rec for rec in gff]
        records.sort(lambda a, b: cmp(a.start, b.start))
        records.sort(lambda a, b: cmp(a.seqid, b.seqid))
        for strand, label in zip(["+", "-"], ["plus", "minus"]):
            bw = pyBigWig.open("{p}_{l}.bw".format(p=bw_file_prefix,
                                                   l=label), "w")
            bw.addHeader(regions)
            seqids = []
            starts = []
            ends = []
            ipds = []
            for rec in records:
                if rec.strand != strand:
                    continue
                assert rec.seqid in valid_seqids, rec.seqid
                seqids.append(rec.seqid)
                starts.append(rec.start)
                ends.append(rec.start + 1)
                ipds.append(float(rec.IPDRatio))
            bw.addEntries(seqids, starts, ends=ends, values=ipds)
            bw.close()
    return 0


def run(gff_file, csv_file, bigwig_file, reference,
        genome_size_cutoff=Constants.MAX_CSV_DEFAULT):
    with ReferenceSet(reference) as ref_ds:
        if ref_ds.totalLength > genome_size_cutoff:
            log.warn("Reference size above cutoff ({l} > {c}), will only include data for modified bases".format(
                l=ref_ds.totalLength, c=genome_size_cutoff))
            return gff2bigwig(gff_file, bigwig_file)
        else:
            return csv2bigwig(csv_file, bigwig_file)


def get_contract_parser():
    p = get_pbparser(
        Constants.TOOL_ID,
        Constants.VERSION,
        "Summarize Consensus",
        __doc__,
        Constants.DRIVER_EXE,
        default_level="ERROR")
    p.add_input_file_type(FileTypes.GFF, "basemods_gff",
                          "Basemods GFF", "Base modifications GFF file")
    p.add_input_file_type(FileTypes.CSV, "basemods_csv"
                          "Basemods CSV", "Base modifications CSV file")
    p.add_input_file_type(FileTypes.DS_REF, "reference",
                          "ReferenceSet XML", "ReferenceSet XML file")
    p.add_output_file_type(FileTypes.BIGWIG, "bigwig_file",
                           name="Bigwig file",
                           description="Bigwig file for SMRTView",
                           default_name="basemods")
    p.add_int(Constants.MAX_CSV_ID, "csvCutoff",
              default=Constants.MAX_CSV_DEFAULT,
              name="Genome size cutoff for displaying per-base information"
              description="Genome size cutoff for displaying per-base information")
    return p


def args_runner(args):
    return run(gff_file=args.basemods_gff,
               csv_file=args.basemods_csv,
               bigwig_file=args.bigwig_file,
               reference=args.reference,
               genome_size_cutoff=args.csvCutoff)


def resolved_tool_contract_runner(rtc):
    return run(gff_file=rtc.task.input_files[0],
               csv_file=rtc.task.input_files[1],
               bigwig_file=rtc.task.output_files[0],
               reference=rtc.task.input_files[2],
               genome_size_cutoff=rtc.task.options[Constants.MAX_CSV_ID])


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_contract_parser(),
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main())
