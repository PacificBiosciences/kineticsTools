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

test:
	nosetests -s -v test/*.py 

pip-install: 
	@which pip > /dev/null
	@pip freeze|grep 'kineticsTools=='>/dev/null \
      && pip uninstall -y pbtools.kineticsTools \
      || pip uninstall -y kineticsTools \
      || true
	@pip install --no-index \
          --install-option="--install-data=$(DATA)" \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: test
.PHONY: clean
