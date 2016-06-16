
import unittest
import tempfile

import numpy as np
import h5py

import pbcommand.testkit.core
from pbcommand.models import PipelineChunk
from pbcommand.pb_io.common import write_pipeline_chunks

from kineticsTools.tasks.gather_kinetics_h5 import gather_kinetics_h5

class SetUpHDF5(object):
    DATALENGTH = 10000
    CHUNKSIZE = 2048
    CHUNKED_FILES = [
        tempfile.NamedTemporaryFile(suffix=".h5").name,
        tempfile.NamedTemporaryFile(suffix=".h5").name
    ]

    @classmethod
    def makeInputs(cls):
        for k, ifn in enumerate(cls.CHUNKED_FILES):
            f = h5py.File(cls.CHUNKED_FILES[k], "w")
            a = f.create_dataset("base", (cls.DATALENGTH,),
                                 dtype="a1", compression="gzip",
                                 chunks=(cls.CHUNKSIZE,), compression_opts=2)
            b = f.create_dataset("score", (cls.DATALENGTH,),
                                 dtype="u4", compression="gzip",
                                 chunks=(cls.CHUNKSIZE,), compression_opts=2)
            c = f.create_dataset("tMean", (cls.DATALENGTH,),
                                 dtype="f4", compression="gzip",
                                 chunks=(cls.CHUNKSIZE,), compression_opts=2)
            start = k * (cls.DATALENGTH / 2)
            end = (k + 1) * (cls.DATALENGTH / 2)
            a[start:end] = np.array(['A' for x in range(start, end)])
            b[start:end] = np.array([2*(k+1) for x in range(start, end)])
            c[start:end] = np.sqrt(np.array(range(start+1, end+1)))
            f.close()


class TestGatherHDF5(unittest.TestCase, SetUpHDF5):

    @classmethod
    def setUpClass(cls):
        cls.makeInputs()

    def test_gather_kinetics_h5(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".h5").name
        gather_kinetics_h5(self.CHUNKED_FILES, ofn)
        f = h5py.File(ofn)
        self.assertTrue(all(f['base'].__array__() == 'A'))
        self.assertTrue(all(f['score'].__array__() > 0))
        self.assertEqual(f['score'].__array__().mean(), 3)
        self.assertTrue(all(f['tMean'].__array__() > 0))
        d = np.round(f['tMean'].__array__() ** 2).astype("u4")
        self.assertTrue(all(d == np.array(range(1, self.DATALENGTH+1))))


class TestGatherH5ToolContract(pbcommand.testkit.core.PbTestGatherApp, SetUpHDF5):
    DRIVER_BASE = "python -m kineticsTools.tasks.gather_kinetics_h5 "
    INPUT_FILES = [tempfile.NamedTemporaryFile(suffix=".json").name]
    CHUNK_KEY = "$chunk.h5_id"

    @classmethod
    def setUpClass(cls):
        super(TestGatherH5ToolContract, cls).setUpClass()
        cls.makeInputs()
        chunks = [PipelineChunk(chunk_id="chunk_data_{i}".format(i=i),
                                **({cls.CHUNK_KEY:fn}))
                  for i, fn in enumerate(cls.CHUNKED_FILES)]
        write_pipeline_chunks(chunks, cls.INPUT_FILES[0], None)


if __name__ == "__main__":
    unittest.main()
