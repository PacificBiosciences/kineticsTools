import os

import pytest


def _internal_data():
    return os.path.exists("/pbi/dept/secondary/siv/testdata")


def pytest_runtest_setup(item):
    for mark in item.iter_markers():
        if mark.name == 'internal_data':
            if not _internal_data():
                pytest.skip(
                    "need access to '/pbi/dept/secondary/siv/testdata'")
        elif mark.name == 'pybigwig':
            try:
                import pyBigWig
            except ImportError:
                pytest.skip("requires pyBigWig to be installed")
        else:
            raise LookupError("Unknown pytest mark: '{}'".format(mark.name))
