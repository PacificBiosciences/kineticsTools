Test detection and identification modes of ipdSummary.py using .bam file as input, with lossless encoding of pulse information.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/mnt/secondary-siv/testdata/kineticsTools
  $ INPUT=$DATA/Mjan_1_5000_lossless.bam
  $ export REF_DIR=/mnt/secondary-siv/references
  $ export REF_SEQ=${REF_DIR}/Methanocaldococcus_jannaschii_DSM2661/sequence/Methanocaldococcus_jannaschii_DSM2661.fasta
  $ export PBCORE_BAM_LOSSLESS_KINETICS=1

Run basic ipdSummary.py:

  $ ipdSummary.py --gff tmp1.gff --csv tmp1.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C --reference $REF_SEQ $INPUT

Look at output csv file:

  $ head -3 tmp1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|6626255|gb|L77117.1|",1,0,T,2,0.737,0.334,0.870,0.848,14
  "gi|6626255|gb|L77117.1|",1,1,A,8,1.338,0.475,0.825,1.621,15

  $ linecount tmp1.csv
  12182

Look at output gff file:

  $ linecount tmp1.gff
  64
  $ cat tmp1.gff | head -20
  ##gff-version 3
  ##source ipdSummary.py * (glob)
  ##source-commandline * (glob)
  ##sequence-region gi|6626255|gb|L77117.1| 1 1664970
  ##sequence-region gi|1500644|gb|L77118.1|MII2CG 1 58407
  ##sequence-region gi|1500688|gb|L77119.1|MII1CG 1 16550
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t424\t424\t172\t-\t.\tcoverage=116;context=TACGCTACATACGCTCCTCCATCTAATGCAGGAGCAAATTT;IPDRatio=4.82;identificationQv=154 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t738\t738\t33\t+\t.\tcoverage=169;context=AATTAAAATCAGACCGTTTCGGAATGGAAAGATTAATGTAA;IPDRatio=1.62 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t845\t845\t267\t-\t.\tcoverage=186;context=TGATTTTAATTTTGATTTCCATCGTGAAGTAATCCAAGTCG;IPDRatio=7.82;identificationQv=246 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t916\t916\t38\t+\t.\tcoverage=210;context=TTGATGCTTTATTTGGATGTTTGAAGAATTAAAATCAGACC;IPDRatio=1.75 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t929\t929\t31\t-\t.\tcoverage=213;context=TCCATTCCGAAACGGTCTGATTTTAATTCTTCAAACATCCA;IPDRatio=1.46 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t988\t988\t41\t+\t.\tcoverage=197;context=CCTTCAAGTTCAAAATTTTTCCCCATAATTAAAATCAGACC;IPDRatio=1.64 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tm6A\t1126\t1126\t300\t-\t.\tcoverage=187;context=GGTCTGATTTTAATTGTTCCATCAACAAACGCAGTTATTGA;IPDRatio=5.91;identificationQv=284 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1186\t1186\t35\t+\t.\tcoverage=153;context=ACAGCATTATAACAATTGATAGACAGAAATTCAGAATTAAA;IPDRatio=1.82 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1188\t1188\t37\t+\t.\tcoverage=155;context=AGCATTATAACAATTGATAGACAGAAATTCAGAATTAAAAT;IPDRatio=1.65 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1450\t1450\t36\t-\t.\tcoverage=193;context=TTTACCTGTGAGTGTTGTAGTTCCAAGTAGATATTTCCATT;IPDRatio=1.55 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1611\t1611\t37\t+\t.\tcoverage=190;context=ATAATAATACTTTAACGCTTCTTTTAAATTAAAATCAGACC;IPDRatio=1.63 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1633\t1633\t35\t-\t.\tcoverage=175;context=GTAATAATTTCCATTCCGAAACGGTCTGATTTTAATTTAAA;IPDRatio=1.58 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1781\t1781\t33\t+\t.\tcoverage=257;context=AAATCAGACCGTTTCGGAATGGAAATTTTTTATCGAAACCT;IPDRatio=1.41 (esc)
  gi|6626255|gb|L77117.1|\tkinModCall\tmodified_base\t1947\t1947\t33\t+\t.\tcoverage=167;context=TATTTTAAGTTACAAAAAAGTTTATAGTGGTTAATGATTAA;IPDRatio=1.62 (esc)
