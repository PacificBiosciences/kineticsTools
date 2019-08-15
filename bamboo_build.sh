#!/bin/bash

set +vx
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module purge
module load gcc
module load ccache
module load python/2
module load hdf5-tools  # for h5ls
set -vex

which gcc
gcc --version
which python
python --version

export PYTHONUSERBASE=$PWD/build
export PATH=${PYTHONUSERBASE}/bin:${PATH}

#PIP="pip --cache-dir=$bamboo_build_working_directory/.pip"
PIP=pip
if [[ -z ${bamboo_repository_branch_name+x} ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
elif [[ ${bamboo_repository_branch_name} == develop ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
elif [[ ${bamboo_repository_branch_name} == master ]]; then
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/master
else
  WHEELHOUSE=/mnt/software/p/python/wheelhouse/develop
fi
export WHEELHOUSE

rm -rf   build
mkdir -p build/{bin,lib,include,share}
PIP_INSTALL="${PIP} install --no-index --find-links=${WHEELHOUSE}"
#${PIP} install --user scipy # TODO: Delete this line when in our wheelhouse.
$PIP_INSTALL --user -r requirements-ci.txt
$PIP_INSTALL --user -r requirements-dev.txt
#iso8601 xmlbuilder tabulate pysam avro?
$PIP install --user ./

make test
