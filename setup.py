#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages

test_deps = [
    'pytest',
    'pytest-cov',
    'pytest-xdist',
]

setup(
    name='kineticsTools',
    version='0.7.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license='BSD-3-Clause-Clear',
    packages=find_packages(),
    include_package_data=True,
    exclude_package_data={'kineticsTools': ['tree_predict.c']},
    ext_modules=[Extension('kineticsTools/tree_predict', ['kineticsTools/tree_predict.c'],
                           extra_compile_args=['-O3', '-shared', '-std=c99'],
                           export_symbols=['innerPredict', 'innerPredictCtx', 'init_native'])],
    zip_safe=False,
    entry_points={'console_scripts': [
        'ipdSummary = kineticsTools.ipdSummary:main',
        'summarizeModifications = kineticsTools.summarizeModifications:main',
    ]},
    setup_requires=[
        'pytest-runner',
    ],
    install_requires=[
        'numpy >= 1.17',
        'pbcommand >= 2.0.0',
        'pbcore >= 2.0.0',
        'pyBigWig',
        'scipy >= 1.3',
    ],
    tests_require=test_deps,
    extras_require={'test': test_deps},
    python_requires='>=3.7',
)
