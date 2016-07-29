from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='kineticsTools',
    version='0.6.0',
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license=open('LICENSES.txt').read(),
    packages=find_packages("."),
    package_data={'kineticsTools': ['resources/*.h5']},
    ext_modules=[Extension('kineticsTools/tree_predict', ['kineticsTools/tree_predict.c'],
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    zip_safe=False,
    install_requires=[
        'pbcore >= 1.2.8',
        'numpy >= 1.6.0',
        'h5py >= 1.3.0',
        'scipy >= 0.9.0',
        'pbcommand >= 0.3.22',
        #'pyBigWig'
    ],
    entry_points={'console_scripts': [
        "ipdSummary = kineticsTools.ipdSummary:main",
        "summarizeModifications = kineticsTools.summarizeModifications:main",
    ]},
)
