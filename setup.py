from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name = 'pbtools.kineticsTools',
    version='0.3.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license='LICENSE.txt',
    scripts = ['src/ipdSummary.py', 'src/summarizeModifications.py', 'src/copyIpdSummaryDataset.py'],
    packages = find_packages('src'),  
    package_dir = {'':'src'},
    namespace_packages = ['pbtools'],
    data_files= [('pbtools/kineticsTools/', ['src/pbtools/kineticsTools/kineticLut.h5'])],
    ext_modules=[Extension('tree_predict', ['src/C/tree_predict.c'], extra_compile_args=["-O3","-shared", "-std=c99"], export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    zip_safe = False,
    install_requires=[
        'pbcore >= 0.2.0',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0',
        'scipy >= 0.9.0'
        ]
    )
