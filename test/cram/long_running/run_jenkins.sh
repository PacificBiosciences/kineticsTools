#!/bin/bash

BASE_PATH=$1
XML_DEST=$2

if [ -z "${BASE_PATH}" ] || [ ! -d "${BASE_PATH}" ]; then
  echo "Base path required as first argument"
  exit 1
fi
if [ -z "${XML_DEST}" ] || [ ! -d "${XML_DEST}" ]; then
  echo "XML output file required as first argument"
  exit 1
fi


cd ${BASE_PATH}
virtualenv ${BASE_PATH}/venv
${BASE_PATH}/venv/bin/pip install CramUnit
${BASE_PATH}/venv/bin/python ${BASE_PATH}/venv/bin/run_cram_unit.py -x ${XML_DEST} ${BASE_PATH}/tests/cram/long_running
