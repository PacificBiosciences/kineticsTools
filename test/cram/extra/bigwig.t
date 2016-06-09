Test support for BigWig output of IpdRatio.

  $ . $TESTDIR/../portability.sh

Load in data:

  $ DATA=/pbi/dept/secondary/siv/testdata/kineticsTools
  $ INPUT=$DATA/Hpyl_1_5000.xml
  $ REFERENCE=/pbi/dept/secondary/siv/references/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

Run basic ipdSummary:

  $ ipdSummary --log-level=WARNING --bigwig tmp1.bw --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --reference $REFERENCE --referenceWindows="gi|12057207|gb|AE001439.1|:0-5000" $INPUT
  $ ls tmp1.bw
  tmp1.bw
