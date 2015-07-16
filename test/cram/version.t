A simple test of the version and help options:

  $ ipdSummary.py --version
  2.2

  $ ipdSummary.py
  usage: ipdSummary.py [-h] [-v] [--reference REFERENCE] [--gff GFF] [--csv CSV]
                       [--identify IDENTIFY] [--methylFraction]
                       [--maxLength MAXLENGTH] [--pvalue PVALUE]
                       [--numWorkers NUMWORKERS] [--outfile OUTFILE]
                       [--m5Cgff M5CGFF] [--m5Cclassifer M5CCLASSIFIER]
                       [--csv_h5 CSV_H5] [--pickle PICKLE]
                       [--summary_h5 SUMMARY_H5] [--ms_csv MS_CSV]
                       [--control CONTROL] [--useLDA] [--paramsPath PARAMSPATH]
                       [--minCoverage MINCOVERAGE] [--maxQueueSize MAXQUEUESIZE]
                       [--maxCoverage MAXCOVERAGE]
                       [--mapQvThreshold MAPQVTHRESHOLD] [--ipdModel IPDMODEL]
                       [--modelIters MODELITERS]
                       [--cap_percentile CAP_PERCENTILE]
                       [--methylMinCov METHYLMINCOV]
                       [--identifyMinCov IDENTIFYMINCOV]
                       [--maxAlignments MAXALIGNMENTS]
                       [-w REFERENCEWINDOWSASSTRING]
                       [--refContigIndex REFCONTIGINDEX]
                       [-W REFERENCEWINDOWSASSTRING]
                       [--skipUnrecognizedContigs SKIPUNRECOGNIZEDCONTIGS]
                       [--threaded] [--profile] [--debug] [--seed RANDOMSEED]
                       [--verbose]
                       [--resolved-tool-contract RESOLVED_TOOL_CONTRACT]
                       aligned.subreads.xml
  ipdSummary.py: error: too few arguments
  [2]
