

======================
Regional m5C Detection
======================

This document describes a simple workflow for obtaining basic estimates of hypo-methylated regions in eukaryotic data.   This application is based on and borrows heavily from the many original ideas of [1].

Following [1], we define a hypo-methylated region here to be a genomic region of any length containing at least 50 CG sites that are all (or mostly) judged by a bisulfite sequencing caller to be less than 50% methylated.


Requirements
------------

To estimate hypo-methylated regions as described here, the following two scripts are required:

- MaximumScoringSubsequences.py
- runMaxScoringSubsequences.py

In addition, the following two m5C classifiers are available to choose from:

- p4_c2_medaka_2_binary_classifier.csv
- p4_c2_arabidopsis_2_2_binary_classifier.csv


The workflow is currently two steps.  The first step requires running ipdSummary.py, the basic modification detection tool, from the command line with special arguments.   The second step requires running one of the additional scripts provided.

The scripts listed above have no additional dependencies beyond those required for kineticsTools.

An additional requirement relates to the chemistry used to collect data to which this workflow can be applied.  Since kinetics generally vary from chemistry to chemistry, models must be applied to data collected using the same chemistry as the training data.  Both models provided were trained using P4-C2 data.   **These P4-C2 models are not expected to work on P5-C3 chemistry data.**



Step 1. Application of LDA classifier
-------------------------------------


The first step involves running ipdSummary.py from the command line with the following arguments::

    ipdSummary.py  --useLDA
                   --refContigIndex <reference index>
                   --m5Cclassifier <csv file containing m5C binary classifier weights>
                   --m5Cgff  <m5C scores gff>
                   --reference <FASTA file>
                   <cmp.h5 file>

The current scheme is best suited to be run on one contig of the reference at a time: <reference index> specifies a contig number indexed from 1.   This is due to a limitation on the regional detection script, runMaxScoringSubsequences.py, which assumes that the input GFF file contains scores corresponding to only one reference.

The score is computed using a weights from an LDA classifier.  The classifier should have two columns stored in a csv format.  Each column should contain 127 values:  126 weights and one offset.  The first column contains weights for data mapped to the forward strand, and the second column contains weights for the reverse strand.  There are two classifiers provided.

The file <m5C scores gff> specifies the location of the output file, which will assign one score to each CG site in the specified contig of the reference.



Training features
~~~~~~~~~~~~~~~~~

The authors of [1] suggest training a classifier to help assign rough scores to each individual CG site.  Following their approach, we select training features that can be derived using the standard modifications.csv output of ipdSummary.py:

The modifications.csv file contains the following statistics for each site and each strand:  tMean (mean IPD), tErr (standard error of IPDs), cov (coverage), modelPrediction (predicted mean IPD for this sequence context at that site).   From these we can compute:

1. standard deviation
2. difference between tMean and modelPrediction
3. alternative error estimate obtained by adding standard deviation and model error in quadrature
4. t-statistic scaled to remove dependence on coverage

The predictors we chose are tMean, modelPrediction, standard deviation, scaled t-statistic, difference between tMean and modelPrediction, alternative error estimate values for every site in a window [-10, +10] around the site of interest.   Prior to training, a log-transformation is applied to all of these, except for the scaled t-statistics.

Following [1], we assume that there is generally some concordance in the methylation status of CG sites on the two strands, and that kinetic information for the two strands can be combined for the two strands.  Our method is to apply a different classifier to each strand, and then add the corresponding classification scores of the two strands for each CG site.



Training labels and data
~~~~~~~~~~~~~~~~~~~~~~~~

The training labels were obtained from bisulfite sequencing data.  A CG site was labeled un-methylated if bisulfite sequencing estimated that less that 0.5 of the molecules contained m5C.

The authors of [1] kindly provided tables summarizing their bisulfite sequencing data for the Medaka genome, as well as RSII sequencing data that was used for training the model stored in the file p4_c2_medaka_2_binary_classifier.csv.  In addition, we provide an alternative model for which the data was collected internally: p4_c2_arabidopsis_2_2_binary_classifier.csv.   The bisulfite sequencing results for this data were provided by Chongyuan Luo in the Ecker Lab at the Salk Institute for Biological Studies.



Output containing m5C scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GFF file < m5C scores gff > will have one row for every CG site in the specified contig of the reference.

1. The first column will list the reference index, and the second and third columns will be 'kinModCall' and 'CG', respectively.
2. The fourth and fifth columns should be identical, listing the template position of the CG site.
3. The sixth column will contain a score that is obtained by summing the LDA scores for forward and reverse strands.   A negative score is expected to loosely correspond to un-methylated sites.
4. The seventh column should only contain '+', for the forward strand.
5. The ninth column contains some extra information:  the coverage and IPD ratio on the forward strand at that CG site.




Step 2.  Estimation of hypo-methylated regions
----------------------------------------------

The second step involves taking the output of step 1 as the input to the additional scripts provided::

    runMaxScoringSubsequences.py   --infile <m5C scores gff>	
                                   --outfile <m5C regions output gff>

We apply the general method of [2] to the individual CG site scores obtained in Step 1 to estimate regions of hypo-methylation.  

The authors of [1] have developed an new method of boundary estimation that is specialized to this application and may yield superior results.   Their implementation is available here: https://github.com/hacone/AgIn



Output containing estimates of hypo-methylated regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The GFF file <m5C regions output gff> will have one row for hypo-methylated region in the specified contig of the reference.

Following [1], we assume that hypo- and hyper-methylated regions alternate.

1. The first column will list the reference index, and the second and third columns will be 'region' and 'hypomethylated', respectively.
2. The fourth and fifth columns contain start and stop positions of the hypo-methylated region.
3. The sixth column contains the negative sum of scores of CG sites in that region.   
4. The seventh column should only contain '+', for the forward strand.   
5. The ninth column contains some extra information:  the number of CG sites in that region, as well as the average coverage and IPD ratio of CG sites in that region on the forward strand.

Once again, we follow [1] and report only regions containing at least 50 CG sites.   


References
----------

1. Suzuki, Yuta, Wei Qu, Tatsuya Tsukahara, Stephen W. Turner, Jonas Korlach, Hideaki Yurino, Jun Yoshimura, Hiroyuki Takeda, and Shinichi Morishita, Completing CpG methylation statuses in a vertebrate genome by integrating SMRT sequencing kinetic data, to appear.
2. Ruzzo, Walter L. and Martin Tompa, A Linear Time Algorithm for Finding All Maximal Scoring Subsequences, 7th International Conference on Intelligent Systems for Molecular Biology, Heidelberg, Germany, August 1999.
