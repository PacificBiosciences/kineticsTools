import logging
import os.path as op
import subprocess
import sys
import tempfile

import pytest

from pbcommand.pb_io.common import load_pipeline_chunks_from_json, \
    write_pipeline_chunks
from pbcommand.pb_io.report import load_report_from_json
from pbcommand.models import PipelineChunk

from kineticsTools.tasks.gather_bigwig import gather_bigwig

log = logging.getLogger(__name__)


def _generate_chunk_output_file(i=None):
    import pyBigWig
    records = [
        ("chr1", 4, 5, 0.45),
        ("chr2", 8, 9, 1.0),
        ("chr2", 9, 10, 6.7),
        ("chr1", 1, 2, 1.5),
        ("chr1", 2, 3, 4.5),
        ("chr1", 3, 4, 1.9)
    ]
    fn = tempfile.NamedTemporaryFile(suffix=".bw").name
    _records = records[(i * 3):(i * 3) + 3]
    assert len(_records) == 3
    ranges = {}
    for rec in _records:
        seqid = rec[0]
        pos = rec[1]
        ranges.setdefault(seqid, (sys.maxsize, 0))
        ranges[seqid] = (min(ranges[seqid][0], pos),
                         max(ranges[seqid][1], pos))
    bw = pyBigWig.open(fn, "w")
    regions = [(s, ranges[s][1] + 1) for s in sorted(ranges.keys())]
    bw.addHeader(regions)
    bw.addEntries([rec[0] for rec in _records],
                  [rec[1] - 1 for rec in _records],
                  ends=[rec[2] - 1 for rec in _records],
                  values=[rec[3] for rec in _records])
    bw.close()
    return fn


@pytest.mark.pybigwig
class TestGatherBigwig:
    NCHUNKS = 2

    @classmethod
    def setup_class(cls):
        cls._data_files = [_generate_chunk_output_file(i=i)
                           for i in range(cls.NCHUNKS)]

    def validate_output(self, output_file):
        import pyBigWig
        bw = pyBigWig.open(output_file)
        nrec = bw.header()["nBasesCovered"]
        assert nrec == 6
        assert pytest.approx(bw.stats("chr1", 2, 3)[0]) == 1.9
        assert pytest.approx(bw.stats("chr2", 7, 8)[0]) == 1.0

    def test_gather_bigwig(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".bw").name
        rc = gather_bigwig(self._data_files, ofn)
        self.validate_output(ofn)

    def test_gather_bigwig_cli(self):
        tmp_dir = tempfile.mkdtemp()
        ofn = op.join(tmp_dir, "gathered.bw")
        args = ["python3", "-m", "kineticsTools.tasks.gather_bigwig",
                ofn] + self._data_files
        log.info("Output will be in %s", tmp_dir)
        with open(op.join(tmp_dir, "stdout"), "w") as stdout:
            with open(op.join(tmp_dir, "stderr"), "w") as stderr:
                subprocess.check_call(args, stdout=stdout, stderr=stderr)
        self.validate_output(ofn)
