Test detection and identification modes of ipdSummary using .bam file as input.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/mnt/secondary-siv/testdata/kineticsTools
  $ INPUT=$DATA/Hpyl_1_5000.bam
  $ REFERENCE=/mnt/secondary-siv/references/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

Run basic ipdSummary:

  $ ipdSummary --gff tmp1.gff --csv tmp1.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --reference $REFERENCE --referenceWindows="gi|12057207|gb|AE001439.1|:0-5000" $INPUT

Look at output csv file:

  $ head -3 tmp1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|12057207|gb|AE001439.1|",1,0,A,9,2.308,0.471,1.710,1.350,29
  "gi|12057207|gb|AE001439.1|",1,1,T,1,0.506,0.077,0.602,0.841,57

  $ linecount tmp1.csv
  10001

Look at output gff file:

  $ linecount tmp1.gff
  276
  $ cat tmp1.gff | head -20
  ##gff-version 3
  ##source ipdSummary v2.0
  ##source-commandline /home/nechols/py/bin/ipdSummary --gff tmp1.gff --csv tmp1.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --reference /mnt/secondary-siv/references/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta --referenceWindows=gi|12057207|gb|AE001439.1|:0-5000 /mnt/secondary-siv/testdata/kineticsTools/Hpyl_1_5000.bam
  ##sequence-region gi|12057207|gb|AE001439.1| 1 1643831
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t35\t35\t186\t-\t.\tcoverage=118;context=TTTAAGGGCGTTTTATGCCTAAATTTAAAAAATGATGCTGT;IPDRatio=5.67;identificationQv=195 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t60\t60\t51\t-\t.\tcoverage=112;context=AAAAAGCTCGCTCAAAAACCCTTGATTTAAGGGCGTTTTAT;IPDRatio=2.64;identificationQv=32 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t89\t89\t222\t+\t.\tcoverage=139;context=AGCGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAA;IPDRatio=5.58;identificationQv=187 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t91\t91\t218\t-\t.\tcoverage=143;context=ATTTTTAAACGCTATCTTGGATTCTTTGAGCAAAAAGCTCG;IPDRatio=6.36;identificationQv=217 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t113\t113\t39\t+\t.\tcoverage=132;context=CAAGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGT;IPDRatio=1.68 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t115\t115\t34\t+\t.\tcoverage=147;context=AGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGTAG;IPDRatio=1.88 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t122\t122\t220\t-\t.\tcoverage=158;context=GCAAACTCTACGCTGAGCCTAACACCCCTAAATTTTTAAAC;IPDRatio=6.48;identificationQv=202 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t232\t232\t225\t+\t.\tcoverage=173;context=AGCGTAAAATCGCCTTTTCCATGCTCCCTAATCGCTTGAAA;IPDRatio=5.94;identificationQv=212 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t233\t233\t281\t-\t.\tcoverage=183;context=ATTTCAAGCGATTAGGGAGCATGGAAAAGGCGATTTTACGC;IPDRatio=6.42;identificationQv=260 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t241\t241\t268\t+\t.\tcoverage=178;context=TCGCCTTTTCCATGCTCCCTAATCGCTTGAAATCCCAGTCT;IPDRatio=5.58;identificationQv=234 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t248\t248\t241\t-\t.\tcoverage=185;context=ATTTAAAAGACTGGGATTTCAAGCGATTAGGGAGCATGGAA;IPDRatio=6.17;identificationQv=225 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t274\t274\t231\t-\t.\tcoverage=190;context=TGAGATTGACGCTCTCATCGAACCGCATTTAAAAGACTGGG;IPDRatio=6.85;identificationQv=223 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t277\t277\t270\t+\t.\tcoverage=188;context=AGTCTTTTAAATGCGGTTCGATGAGAGCGTCAATCTCATTG;IPDRatio=7.48;identificationQv=254 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t312\t312\t38\t-\t.\tcoverage=204;context=GCTTTAAGCCTTTTTAATGGCGTGTTAGAAAAAATCAATGA;IPDRatio=1.89;identificationQv=4 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t373\t373\t388\t+\t.\tcoverage=219;context=TAATCTTTTTTTCTTCTAACATGCTGGAAGCGATTTTTTTA;IPDRatio=7.07;identificationQv=345 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t374\t374\t338\t-\t.\tcoverage=221;context=TTAAAAAAATCGCTTCCAGCATGTTAGAAGAAAAAAAGATT;IPDRatio=6.03;identificationQv=324 (esc)

Now try limiting the number of alignments:

  $ ipdSummary --gff tmp2.gff --csv tmp2.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --maxAlignments 100 --reference $REFERENCE --referenceWindows="gi|12057207|gb|AE001439.1|:0-5000" $INPUT

  $ N_DIFF=`diff tmp1.gff tmp2.gff | wc --lines`
  $ python -c "assert 100 < ${N_DIFF}, ${N_DIFF}"
