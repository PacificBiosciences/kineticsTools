#!/bin/bash -ex

source /mnt/software/Modules/current/init/bash
module load openblas
module load hdf5-tools/1.8.11

mkdir -p tmp
/opt/python-2.7.9/bin/python /mnt/software/v/virtualenv/13.0.1/virtualenv.py tmp/venv
source tmp/venv/bin/activate

(cd repos/pbcommand && make install)
export HDF5_DIR=/mnt/software/h/hdf5-tools/1.8.11
export PIP_CACHE=$PWD/.pip_cache
find $PIP_CACHE -name '*-linux_x86_64.whl' -delete
pip install -r requirements-ci.txt
pip install -r requirements-dev.txt
(cd repos/pbcore && make install)

pip install --no-index --install-option="--install-scripts=$PWD/tmp/venv/bin" ./
make test
