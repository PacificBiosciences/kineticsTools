Test basic mode of ipdSummary.py.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary.py:

  $ ipdSummary.py --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "lambda_NEB3011",13190,0,A,1,0.871,0.230,1.125,0.775,3
  "lambda_NEB3011",13190,1,T,1,0.908,0.251,1.126,0.807,3


  $ linecount tmp.csv
  7603

Look at output gff file:

  $ sed 's/\t/ /g' tmp.gff
  ##gff-version 3
  ##source ipdSummary.py v2.0
  ##source-commandline /home/UNIXHOME/dalexander/.virtualenvs/VE/bin/ipdSummary.py --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/lambda/sequence/lambda.fasta /home/UNIXHOME/dalexander/Projects/software/smrtanalysis/bioinformatics/tools/kineticsTools/test/cram/../data/p4-c2-lambda-mod-decode.cmp.h5
  ##sequence-region ref000001 1 48502
  lambda_NEB3011 kinModCall modified_base 14060 14060 31 - . coverage=49;context=ACGTTATTGCGGAACTTACAACCGCTCAGGCATTTGCTGCA;IPDRatio=2.15
  lambda_NEB3011 kinModCall modified_base 14657 14657 34 + . coverage=155;context=CGGCACAGCCGGGCGATGTGCTGCTGTGCTGTTTTGGTTCA;IPDRatio=1.55
  lambda_NEB3011 kinModCall modified_base 14743 14743 34 - . coverage=173;context=TACCTCTCTCGTTTGCTCAGTTGTTCAGGAATATGGTGCAG;IPDRatio=1.54
  lambda_NEB3011 kinModCall m4C 14756 14756 21 - . coverage=165;context=CCATTTGTCGGTGTACCTCTCTCGTTTGCTCAGTTGTTCAG;IPDRatio=1.32;identificationQv=23
  lambda_NEB3011 kinModCall modified_base 14769 14769 32 - . coverage=168;context=GTGTGCGTCGCTGCCATTTGTCGGTGTACCTCTCTCGTTTG;IPDRatio=1.56
  lambda_NEB3011 kinModCall modified_base 14800 14800 32 - . coverage=173;context=GCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTGCGTCG;IPDRatio=1.66
  lambda_NEB3011 kinModCall modified_base 14805 14805 31 - . coverage=167;context=CAGATGCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTG;IPDRatio=1.64
  lambda_NEB3011 kinModCall modified_base 14834 14834 32 + . coverage=166;context=GGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTGGTC;IPDRatio=1.71
  lambda_NEB3011 kinModCall modified_base 14841 14841 32 + . coverage=162;context=ATCTGCCTTTACGGGGATTTACAACGATTTGGTCGCCGCAT;IPDRatio=1.58
  lambda_NEB3011 kinModCall modified_base 14847 14847 36 - . coverage=172;context=AGGTCGATGCGGCGACCAAATCGTTGTAAATCCCCGTAAAG;IPDRatio=1.85
  lambda_NEB3011 kinModCall modified_base 14864 14864 45 - . coverage=166;context=CCCCCGTTTTCACACGAAGGTCGATGCGGCGACCAAATCGT;IPDRatio=1.70
  lambda_NEB3011 kinModCall m4C 14884 14884 30 - . coverage=173;context=CAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGGT;IPDRatio=2.11;identificationQv=19
  lambda_NEB3011 kinModCall modified_base 14885 14885 33 - . coverage=166;context=CCAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGG;IPDRatio=2.52
  lambda_NEB3011 kinModCall modified_base 14909 14909 35 + . coverage=162;context=AGCCATCCGGGCACTGGCCACACAGCTCCCGGCGTTTCGTC;IPDRatio=1.76
  lambda_NEB3011 kinModCall m6A 14983 14983 208 + . coverage=160;context=TTGCCGGGCGGGACGTCAGCACGTCCGGGTTAACGGCGCAG;IPDRatio=6.72;identificationQv=192
  lambda_NEB3011 kinModCall m6A 14992 14992 206 - . coverage=162;context=CTCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCC;IPDRatio=6.58;identificationQv=186
  lambda_NEB3011 kinModCall m4C 14997 14997 62 - . coverage=140;context=AGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGAC;IPDRatio=2.32;identificationQv=11
  lambda_NEB3011 kinModCall m4C 14998 14998 22 - . coverage=172;context=CAGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGA;IPDRatio=1.51;identificationQv=40
  lambda_NEB3011 kinModCall modified_base 15006 15006 40 - . coverage=170;context=CCATCAGGCAGAGTCTCATGTAACTGCGCCGTTAACCCGGA;IPDRatio=1.75
  lambda_NEB3011 kinModCall modified_base 15008 15008 47 + . coverage=162;context=CGGGTTAACGGCGCAGTTACATGAGACTCTGCCTGATGGCG;IPDRatio=1.73
  lambda_NEB3011 kinModCall modified_base 15039 15039 51 - . coverage=169;context=CCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCAGG;IPDRatio=1.83
  lambda_NEB3011 kinModCall m6A 15041 15041 29 - . coverage=161;context=CCCCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCA;IPDRatio=1.54;identificationQv=3
  lambda_NEB3011 kinModCall modified_base 15066 15066 37 + . coverage=161;context=CCCAGAGTCGCCGGGGCCAAGTCAGGTGGCGTATTCCAGAT;IPDRatio=1.70
  lambda_NEB3011 kinModCall modified_base 15074 15074 33 - . coverage=171;context=CCAGGACAATCTGGAATACGCCACCTGACTTGGCCCCGGCG;IPDRatio=1.81
  lambda_NEB3011 kinModCall modified_base 15078 15078 33 - . coverage=168;context=GCCCCCAGGACAATCTGGAATACGCCACCTGACTTGGCCCC;IPDRatio=1.62
  lambda_NEB3011 kinModCall modified_base 15089 15089 32 - . coverage=167;context=CAATGGCGGCAGCCCCCAGGACAATCTGGAATACGCCACCT;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 15095 15095 31 - . coverage=167;context=ATCCGGCAATGGCGGCAGCCCCCAGGACAATCTGGAATACG;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 15118 15118 32 - . coverage=163;context=GGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGGCAG;IPDRatio=1.55
  lambda_NEB3011 kinModCall modified_base 15124 15124 34 - . coverage=170;context=TGCAAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGG;IPDRatio=1.57
  lambda_NEB3011 kinModCall modified_base 15195 15195 32 - . coverage=164;context=AGCACCATACTGGCACCGAGAGAAAACAGGATGCCGGTCAT;IPDRatio=1.59
  lambda_NEB3011 kinModCall modified_base 15239 15239 33 + . coverage=157;context=TGGTGTGGCGCAGATGCTGGCACCGAAAGCCAGAACTCCCC;IPDRatio=1.65
  lambda_NEB3011 kinModCall modified_base 15261 15261 39 - . coverage=167;context=CCGTTATCCGTTGTCTGTATACGGGGAGTTCTGGCTTTCGG;IPDRatio=1.73
  lambda_NEB3011 kinModCall modified_base 15292 15292 39 - . coverage=166;context=ATCCAGTGAGGAGAAATAGGTGTTCTGCTTACCGTTATCCG;IPDRatio=1.57
  lambda_NEB3011 kinModCall m6A 15381 15381 20 - . coverage=169;context=GCCGTGCTGATCTCCTGAGAAACCACGCGTGACCCCACGCG;IPDRatio=1.40;identificationQv=13
  lambda_NEB3011 kinModCall modified_base 15402 15402 32 - . coverage=169;context=TGACCACCGTCCCCTTCGTCTGCCGTGCTGATCTCCTGAGA;IPDRatio=1.46
  lambda_NEB3011 kinModCall modified_base 15566 15566 31 + . coverage=118;context=GAAGGACAACCTGAAGTCCACGCAGTTGCTGAGTGTGATCG;IPDRatio=1.79
  lambda_NEB3011 kinModCall m6A 15704 15704 21 - . coverage=84;context=CCTGCTCACCAGCCCGGAACACCACCGTGACACCGGATATG;IPDRatio=1.48;identificationQv=3
  lambda_NEB3011 kinModCall modified_base 16035 16035 31 - . coverage=38;context=ATCCTGCGCATCCGGATATTAAACGGGCGCGGCGGCAGGTT;IPDRatio=1.94
  lambda_NEB3011 kinModCall m6A 16380 16380 21 - . coverage=10;context=CATTTATCCACATCCGCCGCACCAAGACGTTTCCCCATGCC;IPDRatio=6.02;identificationQv=5
  lambda_NEB3011 kinModCall m6A 16658 16658 21 - . coverage=6;context=GGGCGCTGAAGCTGTAGCGGAACGGCGCGCCATCATCCGGC;IPDRatio=2.78;identificationQv=20


What about the H5 file?

  $ h5ls -r tmp.h5
  /                        Group
  /ref000001               Dataset {48502}
