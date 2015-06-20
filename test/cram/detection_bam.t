Test detection and identification modes of ipdSummary.py using .bam file as input.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/mnt/secondary-siv/testdata/kineticsTools
  $ INPUT=$DATA/Hpyl_1_5000.bam
  $ REFERENCE=/mnt/secondary/Smrtanalysis/current/common/references/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

Run basic ipdSummary.py:

  $ ipdSummary.py --gff tmp1.gff --csv tmp1.csv --numWorkers 4 --pvalue 0.001 --identify m6A,m4C --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|12057207|gb|AE001439.1|",1,0,A,10,2.798,0.785,1.710,1.636,29
  "gi|12057207|gb|AE001439.1|",1,1,T,0,0.430,0.065,0.602,0.715,57

  $ linecount tmp1.csv
  13357

Look at output gff file:

  $ linecount tmp1.gff
  364
  $ cat tmp1.gff | head -20
  ##gff-version 3
  ##source ipdSummary.py * (glob)
  ##source-commandline * (glob)
  ##sequence-region gi|12057207|gb|AE001439.1| 1 1643831
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t18\t18\t36\t-\t.\tcoverage=110;context=CCTAAATTTAAAAAATGATGCTGTAATCACTAATCACTNNN;IPDRatio=2.14 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t35\t35\t138\t-\t.\tcoverage=118;context=TTTAAGGGCGTTTTATGCCTAAATTTAAAAAATGATGCTGT;IPDRatio=9.73;identificationQv=143 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t60\t60\t58\t-\t.\tcoverage=112;context=AAAAAGCTCGCTCAAAAACCCTTGATTTAAGGGCGTTTTAT;IPDRatio=3.82;identificationQv=38 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t89\t89\t160\t+\t.\tcoverage=139;context=AGCGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAA;IPDRatio=7.69;identificationQv=134 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t91\t91\t160\t-\t.\tcoverage=143;context=ATTTTTAAACGCTATCTTGGATTCTTTGAGCAAAAAGCTCG;IPDRatio=11.85;identificationQv=138 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t113\t113\t31\t+\t.\tcoverage=132;context=CAAGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGT;IPDRatio=1.72 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t115\t115\t41\t+\t.\tcoverage=147;context=AGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGTAG;IPDRatio=2.54 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t122\t122\t194\t-\t.\tcoverage=158;context=GCAAACTCTACGCTGAGCCTAACACCCCTAAATTTTTAAAC;IPDRatio=10.48;identificationQv=175 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t143\t143\t23\t+\t.\tcoverage=150;context=AGGCTCAGCGTAGAGTTTGCCAAGCTCTATGCATTCATTGA;IPDRatio=1.65;identificationQv=6 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t232\t232\t162\t+\t.\tcoverage=173;context=AGCGTAAAATCGCCTTTTCCATGCTCCCTAATCGCTTGAAA;IPDRatio=9.14;identificationQv=155 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t233\t233\t201\t-\t.\tcoverage=183;context=ATTTCAAGCGATTAGGGAGCATGGAAAAGGCGATTTTACGC;IPDRatio=9.40;identificationQv=185 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t241\t241\t208\t+\t.\tcoverage=178;context=TCGCCTTTTCCATGCTCCCTAATCGCTTGAAATCCCAGTCT;IPDRatio=7.09;identificationQv=188 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t248\t248\t173\t-\t.\tcoverage=185;context=ATTTAAAAGACTGGGATTTCAAGCGATTAGGGAGCATGGAA;IPDRatio=10.86;identificationQv=172 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t274\t274\t192\t-\t.\tcoverage=190;context=TGAGATTGACGCTCTCATCGAACCGCATTTAAAAGACTGGG;IPDRatio=11.44;identificationQv=189 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t277\t277\t200\t+\t.\tcoverage=188;context=AGTCTTTTAAATGCGGTTCGATGAGAGCGTCAATCTCATTG;IPDRatio=12.27;identificationQv=166 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t278\t278\t34\t+\t.\tcoverage=188;context=GTCTTTTAAATGCGGTTCGATGAGAGCGTCAATCTCATTGA;IPDRatio=1.87 (esc)

Now try limiting the number of alignments:

  $ ipdSummary.py --gff tmp2.gff --csv tmp2.csv --numWorkers 4 --pvalue 0.001 --identify m6A,m4C --maxAlignments 100 --reference $REFERENCE $INPUT

  $ N_DIFF=`diff tmp1.gff tmp2.gff | wc --lines`
  $ python -c "assert 100 < ${N_DIFF}, ${N_DIFF}"
