#!/bin/bash -ex

source /mnt/software/Modules/current/init/bash
mkdir -p tmp
/opt/python-2.7.9/bin/python /mnt/software/v/virtualenv/13.0.1/virtualenv.py tmp/venv
PIP_CACHE=$PWD/.pip_cache
find $PIP_CACHE -name '*-linux_x86_64.whl' -delete || true

source tmp/venv/bin/activate

rsync -avx /mnt/software/a/anaconda2/4.2.0/pkgs/mkl-11.3.3-0/lib/             tmp/venv/lib/
rsync -avx /mnt/software/a/anaconda2/4.2.0/pkgs/numpy-1.11.1-py27_0/bin/      tmp/venv/bin/
rsync -avx /mnt/software/a/anaconda2/4.2.0/pkgs/numpy-1.11.1-py27_0/lib/      tmp/venv/lib/
rsync -avx /mnt/software/a/anaconda2/4.2.0/pkgs/scipy-0.18.1-np111py27_0/lib/ tmp/venv/lib/

module load hdf5-tools/1.8.11
export HDF5_DIR=/mnt/software/h/hdf5-tools/1.8.11
(cd repos/pbcommand && make install)
pip --cache-dir=$PIP_CACHE install -r requirements-ci.txt
pip --cache-dir=$PIP_CACHE install -r requirements-dev.txt
(cd repos/pbcore && make install)

pip install --no-index --install-option="--install-scripts=$PWD/tmp/venv/bin" ./
make test
