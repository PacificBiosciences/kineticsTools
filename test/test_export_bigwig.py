
from tempfile import NamedTemporaryFile
import unittest

from pbcommand.testkit import PbTestApp


try:
    import pyBigWig
except ImportError:
    pyBigWig = None

SKIP_IF_NOT_AVAILABLE = unittest.skipUnless(pyBigWig is not None,
                                            "Can't import pyBigWig")
DATA = op.join(op.dirname(__file__), "data")
GFF = op.join(DATA, "basemods.gff")
CSV = op.join(DATA, "basemods.csv")

#def _tmp_gff(): return NamedTemporaryFile(suffix=".gff").name
#def _tmp_csv(): return NamedTemporaryFile(suffix=".csv").name

@SKIP_IF_NOT_AVAILABLE
class TestCsvInput(PbTestApp):
    TASK_ID = "kineticstools.tasks.export_bigwig"
    DRIVER_BASE = "python -m kineticsTools.export_bigwig"
    INPUT_FILES = [GFF, CSV, NamedTemporaryFile(suffix=".fasta").name]
    TASK_OPTIONS = {"kineticstools.task_options.genome_size_cutoff": 100}

    @classmethod
    def setUpClass(cls):
        with open(cls.INPUT_FILES[-1], "w") as ref_out:
            ref_out.write(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")


@SKIP_IF_NOT_AVAILABLE
class TestGffInput(TestCsvInput):
    TASK_OPTIONS = {"kineticstools.task_options.genome_size_cutoff":10}
