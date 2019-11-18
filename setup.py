from setuptools import setup, Extension, find_packages

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
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    zip_safe=False,
    install_requires=[
        'pbcore >= 2.0.0',
        'numpy >= 1.6.0',
        'scipy >= 0.9.0',
        'pbcommand >= 2.0.0',
        'pyBigWig',
    ],
    entry_points={'console_scripts': [
        "ipdSummary = kineticsTools.ipdSummary:main",
        "summarizeModifications = kineticsTools.summarizeModifications:main",
    ]},
    python_requires='>=3.7',
)
