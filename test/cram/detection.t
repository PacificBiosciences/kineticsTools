Test basic detection mode of ipdSummary.py.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary.py:

  $ ipdSummary.py --pvalue 0.001 --numWorkers 1 --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "lambda_NEB3011",13190,1,T,2,0.952,0.256,1.126,0.846,3
  "lambda_NEB3011",13190,0,A,3,1.119,0.096,1.125,0.995,3

  $ linecount tmp.csv
  7603

Look at output gff file:

  $ sed 's/\t/ /g' tmp.gff
  ##gff-version 3
  ##source ipdSummary.py v2.0
  ##source-commandline /home/UNIXHOME/dalexander/.virtualenvs/VE/bin/ipdSummary.py --pvalue 0.001 --numWorkers 1 --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/lambda/sequence/lambda.fasta /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/p4-c2-lambda-mod-decode.cmp.h5
  ##sequence-region ref000001 1 48502
  lambda_NEB3011 kinModCall modified_base 14060 14060 31 - . coverage=49;context=ACGTTATTGCGGAACTTACAACCGCTCAGGCATTTGCTGCA;IPDRatio=2.16
  lambda_NEB3011 kinModCall modified_base 14657 14657 35 + . coverage=155;context=CGGCACAGCCGGGCGATGTGCTGCTGTGCTGTTTTGGTTCA;IPDRatio=1.56
  lambda_NEB3011 kinModCall modified_base 14743 14743 32 - . coverage=173;context=TACCTCTCTCGTTTGCTCAGTTGTTCAGGAATATGGTGCAG;IPDRatio=1.51
  lambda_NEB3011 kinModCall modified_base 14769 14769 32 - . coverage=168;context=GTGTGCGTCGCTGCCATTTGTCGGTGTACCTCTCTCGTTTG;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 14800 14800 34 - . coverage=173;context=GCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTGCGTCG;IPDRatio=1.70
  lambda_NEB3011 kinModCall modified_base 14831 14831 33 + . coverage=161;context=CATGGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTG;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 14834 14834 31 + . coverage=166;context=GGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTGGTC;IPDRatio=1.69
  lambda_NEB3011 kinModCall modified_base 14841 14841 32 + . coverage=162;context=ATCTGCCTTTACGGGGATTTACAACGATTTGGTCGCCGCAT;IPDRatio=1.58
  lambda_NEB3011 kinModCall modified_base 14847 14847 35 - . coverage=172;context=AGGTCGATGCGGCGACCAAATCGTTGTAAATCCCCGTAAAG;IPDRatio=1.86
  lambda_NEB3011 kinModCall modified_base 14864 14864 46 - . coverage=166;context=CCCCCGTTTTCACACGAAGGTCGATGCGGCGACCAAATCGT;IPDRatio=1.71
  lambda_NEB3011 kinModCall modified_base 14909 14909 35 + . coverage=162;context=AGCCATCCGGGCACTGGCCACACAGCTCCCGGCGTTTCGTC;IPDRatio=1.77
  lambda_NEB3011 kinModCall modified_base 14983 14983 210 + . coverage=160;context=TTGCCGGGCGGGACGTCAGCACGTCCGGGTTAACGGCGCAG;IPDRatio=6.76
  lambda_NEB3011 kinModCall modified_base 14991 14991 34 - . coverage=172;context=TCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCCG;IPDRatio=1.63
  lambda_NEB3011 kinModCall modified_base 14992 14992 208 - . coverage=162;context=CTCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCC;IPDRatio=6.67
  lambda_NEB3011 kinModCall modified_base 14997 14997 61 - . coverage=140;context=AGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGAC;IPDRatio=2.35
  lambda_NEB3011 kinModCall modified_base 15006 15006 39 - . coverage=170;context=CCATCAGGCAGAGTCTCATGTAACTGCGCCGTTAACCCGGA;IPDRatio=1.74
  lambda_NEB3011 kinModCall modified_base 15008 15008 47 + . coverage=162;context=CGGGTTAACGGCGCAGTTACATGAGACTCTGCCTGATGGCG;IPDRatio=1.73
  lambda_NEB3011 kinModCall modified_base 15039 15039 51 - . coverage=169;context=CCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCAGG;IPDRatio=1.83
  lambda_NEB3011 kinModCall modified_base 15066 15066 37 + . coverage=161;context=CCCAGAGTCGCCGGGGCCAAGTCAGGTGGCGTATTCCAGAT;IPDRatio=1.68
  lambda_NEB3011 kinModCall modified_base 15074 15074 33 - . coverage=171;context=CCAGGACAATCTGGAATACGCCACCTGACTTGGCCCCGGCG;IPDRatio=1.79
  lambda_NEB3011 kinModCall modified_base 15078 15078 33 - . coverage=168;context=GCCCCCAGGACAATCTGGAATACGCCACCTGACTTGGCCCC;IPDRatio=1.61
  lambda_NEB3011 kinModCall modified_base 15089 15089 32 - . coverage=167;context=CAATGGCGGCAGCCCCCAGGACAATCTGGAATACGCCACCT;IPDRatio=1.56
  lambda_NEB3011 kinModCall modified_base 15095 15095 33 - . coverage=167;context=ATCCGGCAATGGCGGCAGCCCCCAGGACAATCTGGAATACG;IPDRatio=1.59
  lambda_NEB3011 kinModCall modified_base 15118 15118 32 - . coverage=163;context=GGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGGCAG;IPDRatio=1.54
  lambda_NEB3011 kinModCall modified_base 15121 15121 31 - . coverage=149;context=AAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGG;IPDRatio=1.58
  lambda_NEB3011 kinModCall modified_base 15124 15124 35 - . coverage=170;context=TGCAAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGG;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 15195 15195 33 - . coverage=164;context=AGCACCATACTGGCACCGAGAGAAAACAGGATGCCGGTCAT;IPDRatio=1.60
  lambda_NEB3011 kinModCall modified_base 15239 15239 33 + . coverage=157;context=TGGTGTGGCGCAGATGCTGGCACCGAAAGCCAGAACTCCCC;IPDRatio=1.64
  lambda_NEB3011 kinModCall modified_base 15261 15261 38 - . coverage=167;context=CCGTTATCCGTTGTCTGTATACGGGGAGTTCTGGCTTTCGG;IPDRatio=1.73
  lambda_NEB3011 kinModCall modified_base 15291 15291 31 + . coverage=158;context=ACGGATAACGGTAAGCAGAACACCTATTTCTCCTCACTGGA;IPDRatio=1.67
  lambda_NEB3011 kinModCall modified_base 15292 15292 41 - . coverage=166;context=ATCCAGTGAGGAGAAATAGGTGTTCTGCTTACCGTTATCCG;IPDRatio=1.61
  lambda_NEB3011 kinModCall modified_base 15402 15402 33 - . coverage=169;context=TGACCACCGTCCCCTTCGTCTGCCGTGCTGATCTCCTGAGA;IPDRatio=1.48
  lambda_NEB3011 kinModCall modified_base 15566 15566 31 + . coverage=118;context=GAAGGACAACCTGAAGTCCACGCAGTTGCTGAGTGTGATCG;IPDRatio=1.78

  $  linecount tmp.gff
  37

What about the H5 file?

  $ h5ls -r tmp.h5
  /                        Group
  /ref000001               Dataset {48502}
