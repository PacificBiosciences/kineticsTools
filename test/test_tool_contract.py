
"""
Tests for end-to-end tool contract support in kineticsTools.ipdSummary,
including consistency checks on output files.
"""

import unittest
import logging
import os.path

import pbcommand.testkit

DATA_DIR = "/mnt/secondary-siv/testdata/kineticsTools"
REF_DIR = "/mnt/secondary-siv/references/Helicobacter_pylori_J99"

EXPECTED_FILES = [
    'file.csv', # FIXME
    'file.gff',
]

class Constants(object):
    N_LINES_GFF = 338
    N_LINES_CSV = 13357
    INITIAL_LINES_CSV = """\
refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
"gi|12057207|gb|AE001439.1|",1,0,A,10,2.387,0.464,1.710,1.396,29
"gi|12057207|gb|AE001439.1|",1,1,T,1,0.492,0.075,0.602,0.817,57"""
    INITIAL_LINES_GFF = """\
gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t35\t35\t187\t-\t.\tcoverage=118;context=TTTAAGGGCGTTTTATGCCTAAATTTAAAAAATGATGCTGT;IPDRatio=5.68;identificationQv=196
gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t60\t60\t49\t-\t.\tcoverage=112;context=AAAAAGCTCGCTCAAAAACCCTTGATTTAAGGGCGTTTTAT;IPDRatio=2.58;identificationQv=33
gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t89\t89\t223\t+\t.\tcoverage=139;context=AGCGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAA;IPDRatio=5.69;identificationQv=187"""

@unittest.skipUnless(os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))
class TestIpdSummary(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m kineticsTools.ipdSummary "
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
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
        file_names = set(os.listdir(output_dir))
        for file_name in EXPECTED_FILES:
            self.assertTrue(file_name in file_names, "Missing %s" % file_name)
        csv_file = os.path.join(output_dir, EXPECTED_FILES[0])
        gff_file = os.path.join(output_dir, EXPECTED_FILES[1])
        def lc(fn): return len(open(fn).readlines())
        self.assertEqual(lc(gff_file), Constants.N_LINES_GFF)
        self.assertEqual(lc(csv_file), Constants.N_LINES_CSV)
        def head(fn,n): return "\n".join( open(fn).read().splitlines()[0:n] )
        self.assertEqual(head(csv_file, 3), Constants.INITIAL_LINES_CSV)
        def head2(fn,n):
            out = []
            i = 0
            for line in open(fn).read().splitlines():
                if line[0] != '#':
                    out.append(line)
                    i += 1
                if i == n:
                    break
            return "\n".join(out)
        self.assertEqual(head2(gff_file, 3), Constants.INITIAL_LINES_GFF)


@unittest.skipUnless(os.path.isdir(DATA_DIR) and os.path.isdir(REF_DIR),
    "%s or %s not available" % (DATA_DIR, REF_DIR))
class TestSummarizeModifications(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "python -m kineticsTools.summarizeModifications"
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR, "Hpyl_1_5000_modifications.gff"),
        os.path.join(DATA_DIR, "Hpyl_1_5000_alignment_summary.gff"),
    ]


if __name__ == "__main__":
    unittest.main()
