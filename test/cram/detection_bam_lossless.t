Test detection and identification modes of ipdSummary using .bam file as input, with lossless encoding of pulse information.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/pbi/dept/secondary/siv/testdata/kineticsTools
  $ INPUT=$DATA/Mjan_1_5000_lossless.bam
  $ export REF_DIR=/pbi/dept/secondary/siv/references
  $ export REF_SEQ=${REF_DIR}/Methanocaldococcus_jannaschii_DSM2661/sequence/Methanocaldococcus_jannaschii_DSM2661.fasta

Run basic ipdSummary:

  $ ipdSummary --log-level=WARNING --gff tmp1.gff --csv tmp1.csv --numWorkers 4 --pvalue 0.001 --identify m6A,m4C --reference $REF_SEQ --useChemistry P6-C4 $INPUT

Look at output csv file:

  $ head -3 tmp1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|6626255|gb|L77117.1|",1,0,T,3,0.865,0.347,0.870,0.995,14
  "gi|6626255|gb|L77117.1|",1,1,A,8,1.338,0.475,0.825,1.621,15

  $ linecount tmp1.csv
  12182

Look at output gff file:

  $ linecount tmp1.gff
  66
  $ cat tmp1.gff | head -20
  ##gff-version 3
  ##source ipdSummary v2.0
  ##source-commandline * (glob)
  ##sequence-region gi|6626255|gb|L77117.1| 1 1664970
  ##sequence-region gi|1500644|gb|L77118.1|MII2CG 1 58407
  ##sequence-region gi|1500688|gb|L77119.1|MII1CG 1 16550
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t424\t424\t171\t-\t.\tcoverage=116;context=TACGCTACATACGCTCCTCCATCTAATGCAGGAGCAAATTT;IPDRatio=4.81;identificationQv=153 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t738\t738\t33\t+\t.\tcoverage=169;context=AATTAAAATCAGACCGTTTCGGAATGGAAAGATTAATGTAA;IPDRatio=1.63 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t845\t845\t263\t-\t.\tcoverage=186;context=TGATTTTAATTTTGATTTCCATCGTGAAGTAATCCAAGTCG;IPDRatio=7.76;identificationQv=242 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t916\t916\t37\t+\t.\tcoverage=210;context=TTGATGCTTTATTTGGATGTTTGAAGAATTAAAATCAGACC;IPDRatio=1.73 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t929\t929\t31\t-\t.\tcoverage=213;context=TCCATTCCGAAACGGTCTGATTTTAATTCTTCAAACATCCA;IPDRatio=1.46 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t988\t988\t43\t+\t.\tcoverage=197;context=CCTTCAAGTTCAAAATTTTTCCCCATAATTAAAATCAGACC;IPDRatio=1.67 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t994\t994\t32\t+\t.\tcoverage=197;context=AGTTCAAAATTTTTCCCCATAATTAAAATCAGACCGTTTCG;IPDRatio=1.57 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t1126\t1126\t305\t-\t.\tcoverage=187;context=GGTCTGATTTTAATTGTTCCATCAACAAACGCAGTTATTGA;IPDRatio=5.99;identificationQv=289 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1186\t1186\t34\t+\t.\tcoverage=153;context=ACAGCATTATAACAATTGATAGACAGAAATTCAGAATTAAA;IPDRatio=1.81 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1188\t1188\t37\t+\t.\tcoverage=155;context=AGCATTATAACAATTGATAGACAGAAATTCAGAATTAAAAT;IPDRatio=1.65 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1450\t1450\t36\t-\t.\tcoverage=193;context=TTTACCTGTGAGTGTTGTAGTTCCAAGTAGATATTTCCATT;IPDRatio=1.55 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1611\t1611\t38\t+\t.\tcoverage=190;context=ATAATAATACTTTAACGCTTCTTTTAAATTAAAATCAGACC;IPDRatio=1.64 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1633\t1633\t33\t-\t.\tcoverage=175;context=GTAATAATTTCCATTCCGAAACGGTCTGATTTTAATTTAAA;IPDRatio=1.56 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1781\t1781\t33\t+\t.\tcoverage=257;context=AAATCAGACCGTTTCGGAATGGAAATTTTTTATCGAAACCT;IPDRatio=1.41 (esc)
