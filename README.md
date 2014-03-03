kineticsTools
=============

Tools for detecting DNA modifications from single molecule, real-time (SMRT&reg;) sequencing data. This tool implements the P_ModificationDetection module in SMRT&reg; Portal, used by the RS_Modification_Detection and RS_Modifications_and_Motif_Detection protocol. Researchers interested in understanding or extending the modification detection algorthims can use these tools as a starting point. 

The current fork has an additional feature to analyze the base modifications for individual single molecules and all analogous code have a prefix "sm" in the python script name, eg smIpdSummary.py. 

The only files needed if you already have PacificBiosciences/kineticsTools installed are all files named sm*.py 
and the setup.py. Place them in the same directories as the analogous standard modules...eg place smIpdSummary.py 
in the same location as ipdSummary.py These are the files you'll need

setup.py
smIpdSummary.py
smKineticWorker.py
smResultWorker.py

To run the analysis, the command is analogous to running ipdSummary.py.

smIpdSummary.py --reference MNdeI.fasta --csv smMods.csv --ipdModel ~/analysis/etc/algorithm_parameters/2013-09/kineticsTools/C2.h5 --smBaseMod=True --errWin=True --minCoverage=5 MNdeI_case.cmp.h5

The output file is currently a csv file simimlar to the standard multi-molecule approach except this mode will generate an extra column named moleculeID. Hence, each molecule will have it's own set of base modification information at every position and strand that meets the user specified minCoverage.

Academic Publications:
 * [Characterization of DNA methyltransferase specificities using single-molecule, real-time DNA sequencing](http://nar.oxfordjournals.org/content/40/4/e29)
 * [The methylomes of six bacteria](http://nar.oxfordjournals.org/content/early/2012/10/02/nar.gks891.full)

Documentation:
 * [Tool documentation](http://github.com/PacificBiosciences/kineticsTools/blob/master/doc/manual.rst)
 * [Methods description](http://github.com/PacificBiosciences/kineticsTools/blob/master/doc/whitepaper/kinetics.pdf)
