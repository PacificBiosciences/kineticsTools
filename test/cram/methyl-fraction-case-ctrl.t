Test basic mode of ipdSummary.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary:

  $ ipdSummary --numWorkers 1 --methylFraction --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --control $INPUT --reference $REFERENCE $INPUT >/dev/null

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,pvalue,caseMean,controlMean,caseStd,controlStd,ipdRatio,testStatistic,coverage,controlCoverage,caseCoverage,frac,fracLow,fracUp
  "lambda_NEB3011",*,?,?,*,*,*,*,*,*,1.000,0.0,*,*,*,,, (glob)
  "lambda_NEB3011",*,?,?,*,*,*,*,*,*,1.000,0.0,*,*,*,,, (glob)

Look at output gff file:

  $ cat tmp.gff
  ##gff-version 3
  ##source ipdSummary * (glob)
  ##source-commandline * (glob)
  ##sequence-region lambda_NEB3011 1 48502


What about the H5 file?

  $ h5ls -r tmp.h5
  /                        Group
  /lambda_NEB3011          Dataset {48502}
