#!/usr/bin/env python3

"""
One-time utility to convert our ipdModel HDF5 files to NumPy's .npz format
"""

import argparse
import re
import sys

import h5py
import numpy as np

def run(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("file_name", help="Path to hdf5 file")
    args = p.parse_args(argv)
    assert args.file_name.endswith(".h5")
    npz_file = re.sub(".h5", ".npz", args.file_name)
    with h5py.File(args.file_name, 'r') as f:
        k = list(f.keys())[0]
        d = {k2:f[k][k2][:] for k2 in f[k].keys()}
        np.savez(npz_file, **d)
        print("Wrote {o}".format(o=npz_file))
    return 0

if __name__ == "__main__":
    sys.exit(run(sys.argv[1:]))
