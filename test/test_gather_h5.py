
# TODO test with more than one simulated reference contig

import unittest
import tempfile

import numpy as np
import h5py

import pbcommand.testkit.core
from pbcommand.models import PipelineChunk
from pbcommand.pb_io.common import write_pipeline_chunks

from kineticsTools.tasks.gather_kinetics_h5 import (gather_kinetics_h5,
                                                    gather_kinetics_h5_byref)

class SetUpHDF5(object):
    DATALENGTH = 10000
    CHUNKSIZE = 2048
    CHUNKED_FILES = [
        tempfile.NamedTemporaryFile(suffix=".h5").name,
        tempfile.NamedTemporaryFile(suffix=".h5").name
    ]

    @classmethod
    def get_base_group_input(cls, f):
        return f.create_group("chr1")

    def get_base_group_output(self, f):
        return f[f.keys()[0]]

    @classmethod
    def makeInputs(cls):
        for k, ifn in enumerate(cls.CHUNKED_FILES):
            f = h5py.File(cls.CHUNKED_FILES[k], "w")
            g = cls.get_base_group_input(f)
            a = g.create_dataset("base", (cls.DATALENGTH,),
                                 dtype="a1", compression="gzip",
                                 chunks=(cls.CHUNKSIZE,), compression_opts=2)
            b = g.create_dataset("score", (cls.DATALENGTH,),
                                 dtype="u4", compression="gzip",
                                 chunks=(cls.CHUNKSIZE,), compression_opts=2)
            c = g.create_dataset("tMean", (cls.DATALENGTH,),
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

    def _gather(self, ofn):
        gather_kinetics_h5_byref(self.CHUNKED_FILES, ofn)

    def test_gather_kinetics_h5(self):
        ofn = tempfile.NamedTemporaryFile(suffix=".h5").name
        self._gather(ofn)
        f = h5py.File(ofn)
        g = self.get_base_group_output(f)
        self.assertTrue(all(g['base'].__array__() == 'A'))
        self.assertTrue(all(g['score'].__array__() > 0))
        self.assertEqual(g['score'].__array__().mean(), 3)
        self.assertTrue(all(g['tMean'].__array__() > 0))
        d = np.round(g['tMean'].__array__() ** 2).astype("u4")
        self.assertTrue(all(d == np.array(range(1, self.DATALENGTH+1))))


class TestGatherHDF5Flat(TestGatherHDF5):

    @classmethod
    def get_base_group_input(cls, f):
        return f

    def get_base_group_output(self, f):
        return f

    def _gather(self, ofn):
        gather_kinetics_h5(self.CHUNKED_FILES, ofn)


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
