SHELL = /bin/bash -e

utest:
	PYTHONPATH=.:${PYTHONPATH} py.test -s -v test/test_internal.py

all: build install

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

clean:
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	find . -name "*.pyc" | xargs rm -rf;\
	rm -rf dist/

test: tests
check: tests
tests: cram-tests py-tests extra-tests

cram-tests:
	cram --xunit-file=cramtests.xml test/cram/*.t

long-tests:
	cram test/cram/long_running/*.t

py-tests:
	# pytest --cov=kineticsTools  # does not quite work since we run in test/ dir.
	#cd test/; pytest -s -v -n auto --dist=loadscope --durations=20 --junitxml=../nosetests.xml --cov-report=xml:../coverage.xml test_*.py
	cd test/; pytest -s -v --durations=20 --junitxml=../nosetests.xml --cov-report=xml:../coverage.xml test_*.py

extra-tests:
	cram --xunit-file=cramtests-extra.xml test/cram/extra/*.t

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'kineticsTools=='>/dev/null \
      && ( pip uninstall -y kineticsTools \
        || pip uninstall -y pbtools.kineticsTools ) \
      || true
	@pip install --no-index \
          --install-option="--install-data=$(DATA)" \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: tests test clean cram-tests unit-tests
.PHONY: clean
