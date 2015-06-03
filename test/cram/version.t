A simple test of the version and help options:

  $ ipdSummary.py --version
  2.2

  $ ipdSummary.py
  usage: ipdSummary.py [-h] [--reference REFERENCE] [--outfile OUTFILE]
                       [--gff GFF] [--m5Cgff M5CGFF]
                       [--m5Cclassifer M5CCLASSIFIER] [--csv CSV]
                       [--csv_h5 CSV_H5] [--pickle PICKLE]
                       [--summary_h5 SUMMARY_H5] [--ms_csv MS_CSV]
                       [--control CONTROL] [--identify IDENTIFY]
                       [--methylFraction] [--useLDA] [--paramsPath PARAMSPATH]
                       [--maxLength MAXLENGTH] [--minCoverage MINCOVERAGE]
                       [--maxQueueSize MAXQUEUESIZE] [--maxCoverage MAXCOVERAGE]
                       [--mapQvThreshold MAPQVTHRESHOLD] [--pvalue PVALUE]
                       [--ipdModel IPDMODEL] [--modelIters MODELITERS]
                       [--cap_percentile CAP_PERCENTILE]
                       [--methylMinCov METHYLMINCOV]
                       [--identifyMinCov IDENTIFYMINCOV]
                       [-w REFERENCEWINDOWSASSTRING]
                       [--refContigIndex REFCONTIGINDEX]
                       [-W REFERENCEWINDOWSASSTRING]
                       [--skipUnrecognizedContigs SKIPUNRECOGNIZEDCONTIGS]
                       [--numWorkers NUMWORKERS] [--threaded] [--profile]
                       [--debug] [--verbose] [--version]
                       input.cmp.h5
  ipdSummary.py: error: too few arguments
  [2]
