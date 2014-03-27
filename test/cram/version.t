A simple test of the version and help options:

  $ ipdSummary.py -v
  2.2

  $ ipdSummary.py -h > tmp.txt
  $ head -19 tmp.txt
  usage: ipdSummary.py [-h] [-v] [--reference REFERENCE] [--outfile OUTFILE]
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
                       [--refContigs REFCONTIGS]
                       [--refContigIndex REFCONTIGINDEX]
                       [--refContigsFile REFCONTIGS] [--numWorkers NUMWORKERS]
                       [--threaded] [--profile] [--pdb]
                       input.cmp.h5
