Test basic mode of ipdSummary.py.

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
  * [INFO] Available CPUs: * (glob)
  * [INFO] Requested worker processes: 1 (glob)
  * [INFO] Launched worker processes. (glob)
  * [INFO] Worker KineticWorkerProcess-? (PID=*) started running (glob)
  * [INFO] Launched result collector process. (glob)
  * [INFO] Process KineticsWriter-? (PID=*) started running (glob)
  * [INFO] Creating IpdRatio dataset w/ name: *, Size: * (glob)
  * [INFO] Generating kinetics summary for [*] (glob)
  * [INFO] Processing reference entry: [?] (glob)
  * [INFO] Queueing chunks for ref: ?.  NumReads: *, Block Size: * (glob)
  * [INFO] Got chunk: (*, (*, *, *)) -- Process: <KineticWorkerProcess(KineticWorkerProcess-?, started daemon)> (glob)
  * [INFO] Making summary: * to * (glob)
  * [INFO] Got chunk: (*, (*, *, *)) -- Process: <KineticWorkerProcess(KineticWorkerProcess-?, started daemon)> (glob)
  * [INFO] Making summary: * to * (glob)
  * [INFO] Process KineticWorkerProcess-? (PID=*) done; exiting. (glob)
  * [INFO] Result thread shutting down... (glob)
  * [INFO] ipdSummary.py finished. Exiting. (glob)

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

  $  wc -l tmp.gff
  190 tmp.gff

What about the H5 file?

  $ [ -f tmp.h5 ] && echo "IPD ratio H5 file found" || echo "IPD ratio H5 file not found"
  IPD ratio H5 file found

  $ wc -l tmp.h5
  46 tmp.h5
