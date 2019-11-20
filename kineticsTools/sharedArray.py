from multiprocessing.sharedctypes import RawArray
import warnings
import numpy as np


class SharedArray:

    """
    Very simple wrapper for a chunk of shared memory that can be accessed across processes
    """

    def __init__(self, dtype, shape):
        self._rawArray = RawArray(dtype, shape)

    def getNumpyWrapper(self):
        """
        Construct a numpy array that wraps the raw shared memory array
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.ctypeslib.as_array(self._rawArray)
