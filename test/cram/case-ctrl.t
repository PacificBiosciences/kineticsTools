Case-control of a job against itself---shouldn't find any differences.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary.py:

  $ ipdSummary.py --numWorkers 1 --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --control $INPUT --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,pvalue,caseMean,controlMean,caseStd,controlStd,ipdRatio,testStatistic,coverage,controlCoverage,caseCoverage
  "lambda_NEB3011",*,?,?,*,*,*,*,*,*,1.000,0.0,3,3,3 (glob)
  "lambda_NEB3011",*,?,?,*,*,*,*,*,*,1.000,0.0,3,3,3 (glob)

  $ linecount tmp.csv
  7603

Look at output gff file:

  $ sed 's/\t/ /g' tmp.gff
  ##gff-version 3
  ##source ipdSummary.py v2.0
  ##source-commandline /home/UNIXHOME/dalexander/.virtualenvs/VE/bin/ipdSummary.py --numWorkers 1 --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --control /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/p4-c2-lambda-mod-decode.cmp.h5 --reference /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/lambda/sequence/lambda.fasta /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/p4-c2-lambda-mod-decode.cmp.h5
  ##sequence-region ref000001 1 48502

What about the IPD ratio H5 file?

  $ [ -f tmp.h5 ] && echo "IPD ratio H5 file found" || echo "IPD ratio H5 file not found"
  IPD ratio H5 file found

  $ h5ls -r tmp.h5
  /                        Group
  /ref000001               Dataset {48502}
