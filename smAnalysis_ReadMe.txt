This version is forked from PacificBiosciences supported master version.

Contains additional files that support base modification analysis to be performed per individual single 
molecules rather than the multi-molecule approach in the standard kineticsTools.

Relevant files are analagous to those in the standard kineticsTools but have prefix 'sm' in front of the file
name, eg smIpdSummary.py vs ipdSummary.py.

The only files needed if you already have PacificBiosciences/kineticsTools installed are all files named sm*.py 
and the setup.py. Place them in the same directories as the analogous standard modules...eg place smIpdSummary.py 
in the same location as ipdSummary.py These are the files you'll need

setup.py
smIpdSummary.py
smKineticWorker.py
smResultWorker.py

To run the analysis, the command is analogous to running ipdSummary.py.

smIpdSummary.py --reference MNdeI.fasta --csv smMods.csv --ipdModel ~/analysis/etc/algorithm_parameters/2013-09/kineticsTools/C2.h5 --smBaseMod=True --errWin=True --minCoverage=5 MNdeI_case.cmp.h5


