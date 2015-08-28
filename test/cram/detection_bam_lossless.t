Test detection and identification modes of ipdSummary using .bam file as input, with lossless encoding of pulse information.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/mnt/secondary-siv/testdata/kineticsTools
  $ INPUT=$DATA/Mjan_1_5000_lossless.bam
  $ export REF_DIR=/mnt/secondary-siv/references
  $ export REF_SEQ=${REF_DIR}/Methanocaldococcus_jannaschii_DSM2661/sequence/Methanocaldococcus_jannaschii_DSM2661.fasta

Run basic ipdSummary:

  $ ipdSummary --gff tmp1.gff --csv tmp1.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --reference $REF_SEQ $INPUT

Look at output csv file:

  $ head -3 tmp1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|6626255|gb|L77117.1|",1,0,T,3,0.844,0.348,0.870,0.970,14
  "gi|6626255|gb|L77117.1|",1,1,A,7,1.204,0.476,0.825,1.458,15

  $ linecount tmp1.csv
  12182

Look at output gff file:

  $ linecount tmp1.gff
  65
  $ cat tmp1.gff | head -20
  ##gff-version 3
  ##source ipdSummary v2.0
  ##source-commandline * (glob)
  ##sequence-region gi|6626255|gb|L77117.1| 1 1664970
  ##sequence-region gi|1500644|gb|L77118.1|MII2CG 1 58407
  ##sequence-region gi|1500688|gb|L77119.1|MII1CG 1 16550
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t424\t424\t170\t-\t.\tcoverage=116;context=TACGCTACATACGCTCCTCCATCTAATGCAGGAGCAAATTT;IPDRatio=4.81;identificationQv=153 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t738\t738\t33\t+\t.\tcoverage=169;context=AATTAAAATCAGACCGTTTCGGAATGGAAAGATTAATGTAA;IPDRatio=1.63 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t845\t845\t264\t-\t.\tcoverage=186;context=TGATTTTAATTTTGATTTCCATCGTGAAGTAATCCAAGTCG;IPDRatio=7.79;identificationQv=244 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t916\t916\t39\t+\t.\tcoverage=210;context=TTGATGCTTTATTTGGATGTTTGAAGAATTAAAATCAGACC;IPDRatio=1.76 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t929\t929\t31\t-\t.\tcoverage=213;context=TCCATTCCGAAACGGTCTGATTTTAATTCTTCAAACATCCA;IPDRatio=1.47 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t988\t988\t43\t+\t.\tcoverage=197;context=CCTTCAAGTTCAAAATTTTTCCCCATAATTAAAATCAGACC;IPDRatio=1.67 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t994\t994\t31\t+\t.\tcoverage=197;context=AGTTCAAAATTTTTCCCCATAATTAAAATCAGACCGTTTCG;IPDRatio=1.56 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t1126\t1126\t300\t-\t.\tcoverage=187;context=GGTCTGATTTTAATTGTTCCATCAACAAACGCAGTTATTGA;IPDRatio=5.92;identificationQv=284 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1186\t1186\t35\t+\t.\tcoverage=153;context=ACAGCATTATAACAATTGATAGACAGAAATTCAGAATTAAA;IPDRatio=1.82 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1188\t1188\t36\t+\t.\tcoverage=155;context=AGCATTATAACAATTGATAGACAGAAATTCAGAATTAAAAT;IPDRatio=1.65 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1450\t1450\t36\t-\t.\tcoverage=193;context=TTTACCTGTGAGTGTTGTAGTTCCAAGTAGATATTTCCATT;IPDRatio=1.55 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1611\t1611\t37\t+\t.\tcoverage=190;context=ATAATAATACTTTAACGCTTCTTTTAAATTAAAATCAGACC;IPDRatio=1.63 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1633\t1633\t35\t-\t.\tcoverage=175;context=GTAATAATTTCCATTCCGAAACGGTCTGATTTTAATTTAAA;IPDRatio=1.58 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1781\t1781\t33\t+\t.\tcoverage=257;context=AAATCAGACCGTTTCGGAATGGAAATTTTTTATCGAAACCT;IPDRatio=1.41 (esc)
