
"""
Tests for end-to-end tool contract support in kineticsTools.ipdSummary,
including consistency checks on output files.
"""

import unittest
import logging
import os.path
import csv

import h5py

import pbcommand.testkit

os.environ["PACBIO_TEST_ENV"] = "1" # turns off --verbose

DATA_DIR = "/pbi/dept/secondary/siv/testdata/kineticsTools"
REF_DIR = "/pbi/dept/secondary/siv/references/Helicobacter_pylori_J99"


class Constants(object):
    N_LINES_GFF = 338
    INITIAL_LINES_GFF = """\
gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t35\t35\t187\t-\t.\tcoverage=118;context=TTTAAGGGCGTTTTATGCCTAAATTTAAAAAATGATGCTGT;IPDRatio=5.68;identificationQv=196
gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t60\t60\t49\t-\t.\tcoverage=112;context=AAAAAGCTCGCTCAAAAACCCTTGATTTAAGGGCGTTTTAT;IPDRatio=2.58;identificationQv=33
gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t89\t89\t223\t+\t.\tcoverage=139;context=AGCGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAA;IPDRatio=5.69;identificationQv=187"""

# FIXME this is much too slow
@unittest.skipUnless(os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))
class TestIpdSummary(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m kineticsTools.ipdSummary "
    REQUIRES_PBCORE = True
    MAX_NPROC = 8
    RESOLVED_NPROC = 8
    INPUT_FILES = [
        os.path.join(DATA_DIR, "Hpyl_1_5000.xml"),
        os.path.join(REF_DIR, "sequence", "Helicobacter_pylori_J99.fasta"),
    ]
    TASK_OPTIONS = {
      "kinetics_tools.task_options.identify": "m6A,m4C",
      "kinetics_tools.task_options.max_length": 3000000000,
      "kinetics_tools.task_options.compute_methyl_fraction": False,
      "kinetics_tools.task_options.pvalue": 0.001
    }

    def run_after(self, rtc, output_dir):
        gff_file = os.path.join(output_dir, rtc.task.output_files[0])
        def lc(fn): return len(open(fn).readlines())
        self.assertEqual(lc(gff_file), Constants.N_LINES_GFF)
        def head2(fn,n):
            out = []
            for line in open(fn).read().splitlines():
                if line[0] != '#':
                    out.append(line)
                if len(out) == n:
                    break
            return "\n".join(out)
        self.assertEqual(head2(gff_file, 3), Constants.INITIAL_LINES_GFF)
        f = h5py.File(rtc.task.output_files[2])
        g = f[f.keys()[0]]
        seqh_fwd = ''.join([g['base'][x*2] for x in range(5000)])
        seqh_rev = ''.join([g['base'][(x*2)+1] for x in range(5000)])
        self.assertEqual(len(seqh_fwd), 5000)
        self.assertEqual(len(seqh_rev), 5000)


# FIXME this is covered by integration tests but the output is twitchy enough
# that it would be good to test here as well
#@unittest.skipUnless(os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
#    "%s or %s not available" % (DATA_DIR, REF_DIR))
@unittest.skip("FIXME")
class TestIpdSummaryChunk(TestIpdSummary):
    """
    This test is identical to the above except for using an AlignmentSet with
    filters applied as input.  We want the output to actually be within the
    range specified by the filters, not a seemingly arbitirary larger range
    around it.
    """
    INPUT_FILES = [
        os.path.join(DATA_DIR, "Hpyl_1_5000_chunk.xml"),
        os.path.join(REF_DIR, "sequence", "Helicobacter_pylori_J99.fasta"),
    ]

    def run_after(self, rtc, output_dir):
        gff_file = os.path.join(output_dir, rtc.task.output_files[0])
        logging.critical(gff_file)
        logging.critical("%s %s" % (csv_file, os.path.getsize(csv_file)))
        with open(csv_file) as f:
            records = [ r for r in csv.DictReader(f) ]
            logging.critical("start=%s end=%s" % (records[0]['tpl'],
                records[-1]["tpl"]))
            self.assertEqual(records[0]["tpl"], "1001")
            self.assertEqual(records[-1]["tpl"], "1050")


@unittest.skipUnless(os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))
class TestSummarizeModifications(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m kineticsTools.summarizeModifications"
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR, "Hpyl_1_5000_modifications.gff"),
        os.path.join(DATA_DIR, "Hpyl_1_5000_alignment_summary.gff"),
    ]


if __name__ == "__main__":
    unittest.main()
