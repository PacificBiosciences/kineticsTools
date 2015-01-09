from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='kineticsTools',
    version='0.5.1',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license=open('LICENSES.txt').read(),
    scripts=['bin/ipdSummary.py', 'bin/summarizeModifications.py', 'bin/copyIpdSummaryDataset.py'],
    packages=find_packages('.'),
    package_dir={'': '.'},
    package_data={'kineticsTools': ['resources/*.h5']},
    ext_modules=[Extension('kineticsTools/tree_predict', ['kineticsTools/tree_predict.c'],
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    zip_safe=False,
    install_requires=[
        'pbcore >= 0.8.0',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0',
        'scipy >= 0.9.0'
    ]
)
