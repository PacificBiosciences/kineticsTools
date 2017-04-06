SHELL = /bin/bash -e

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
tests: cram-tests unit-tests extra-tests

cram-tests:
	cram --xunit-file=cramtests.xml test/cram/*.t

long-tests:
	cram test/cram/long_running/*.t

unit-tests:
	#nosetests -s -v --with-xunit test/*.py
	py.test -s -v --junit-xml=nosetests.xml test/*.py

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
