Test basic mode of ipdSummary.py.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Can we find the input and reference files?

  $ [ -f $INPUT ] && echo "cmpH5 file exists" || echo "cmpH5 file not found"
  cmpH5 file exists

  $ [ -f $REFERENCE ] && echo "fasta file exists" || echo "fasta file not found"
  fasta file exists

Run basic ipdSummary.py:

  $ ipdSummary.py --numWorkers 1 --identify m6A,m4C --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference $REFERENCE $INPUT

Look at output csv file:

  $ [ -f tmp.csv ] && echo "CSV file found" || echo "CSV file not found"
  CSV file found

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "lambda_NEB3011",*,*,*,*,*,*,*,*,* (glob)
  "lambda_NEB3011",*,*,*,*,*,*,*,*,* (glob)


Look at output gff file:

  $ [ -f tmp.gff ] && echo "GFF file found" || echo "GFF file not found"
  GFF file found

  $ head -2 tmp.gff
  ??gff-version * (glob)
  ??source ipdSummary.py v* (glob)

  $ tail -2 tmp.gff
  lambda_NEB3011\tkinModCall\tm6A\t*\t*\t*\t?\t.\tcoverage=*;context=*;IPDRatio=*;identificationQv=* (glob)
  lambda_NEB3011\tkinModCall\tm6A\t*\t*\t*\t?\t.\tcoverage=*;context=*;IPDRatio=*;identificationQv=* (glob)

  $  linecount tmp.gff
  188

What about the H5 file?

  $ [ -f tmp.h5 ] && echo "IPD ratio H5 file found" || echo "IPD ratio H5 file not found"
  IPD ratio H5 file found

  $ h5ls -r tmp.h5
  /                        Group
  /ref000001               Dataset {48502}
