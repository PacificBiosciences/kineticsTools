Test basic detection mode of ipdSummary.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary:

  $ ipdSummary --log-level=WARNING --pvalue 0.001 --numWorkers 1 --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "lambda_NEB3011",13190,1,T,1,0.909,0.252,1.126,0.808,3
  "lambda_NEB3011",13190,0,A,3,1.120,0.096,1.125,0.996,3

  $ linecount tmp.csv
  7603

Look at output gff file:

  $ cat tmp.gff
  ##gff-version 3
  ##source ipdSummary * (glob)
  ##source-commandline * (glob)
  ##sequence-region lambda_NEB3011 1 48502
  lambda_NEB3011\tkinModCall\tmodified_base\t14060\t14060\t31\t-\t.\tcoverage=49;context=ACGTTATTGCGGAACTTACAACCGCTCAGGCATTTGCTGCA;IPDRatio=2.15 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14657\t14657\t34\t+\t.\tcoverage=155;context=CGGCACAGCCGGGCGATGTGCTGCTGTGCTGTTTTGGTTCA;IPDRatio=1.54 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14743\t14743\t35\t-\t.\tcoverage=173;context=TACCTCTCTCGTTTGCTCAGTTGTTCAGGAATATGGTGCAG;IPDRatio=1.54 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14769\t14769\t31\t-\t.\tcoverage=168;context=GTGTGCGTCGCTGCCATTTGTCGGTGTACCTCTCTCGTTTG;IPDRatio=1.56 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14800\t14800\t34\t-\t.\tcoverage=173;context=GCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTGCGTCG;IPDRatio=1.68 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14831\t14831\t32\t+\t.\tcoverage=161;context=CATGGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTG;IPDRatio=1.55 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14834\t14834\t32\t+\t.\tcoverage=166;context=GGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTGGTC;IPDRatio=1.70 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14841\t14841\t33\t+\t.\tcoverage=162;context=ATCTGCCTTTACGGGGATTTACAACGATTTGGTCGCCGCAT;IPDRatio=1.59 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14847\t14847\t35\t-\t.\tcoverage=172;context=AGGTCGATGCGGCGACCAAATCGTTGTAAATCCCCGTAAAG;IPDRatio=1.83 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14864\t14864\t46\t-\t.\tcoverage=166;context=CCCCCGTTTTCACACGAAGGTCGATGCGGCGACCAAATCGT;IPDRatio=1.71 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14884\t14884\t31\t-\t.\tcoverage=173;context=CAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGGT;IPDRatio=2.12 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14885\t14885\t33\t-\t.\tcoverage=166;context=CCAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGG;IPDRatio=2.51 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14909\t14909\t36\t+\t.\tcoverage=162;context=AGCCATCCGGGCACTGGCCACACAGCTCCCGGCGTTTCGTC;IPDRatio=1.76 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14983\t14983\t210\t+\t.\tcoverage=160;context=TTGCCGGGCGGGACGTCAGCACGTCCGGGTTAACGGCGCAG;IPDRatio=6.76 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14991\t14991\t32\t-\t.\tcoverage=172;context=TCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCCG;IPDRatio=1.59 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14992\t14992\t209\t-\t.\tcoverage=162;context=CTCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCC;IPDRatio=6.64 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14997\t14997\t60\t-\t.\tcoverage=140;context=AGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGAC;IPDRatio=2.29 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15006\t15006\t39\t-\t.\tcoverage=170;context=CCATCAGGCAGAGTCTCATGTAACTGCGCCGTTAACCCGGA;IPDRatio=1.74 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15008\t15008\t47\t+\t.\tcoverage=162;context=CGGGTTAACGGCGCAGTTACATGAGACTCTGCCTGATGGCG;IPDRatio=1.74 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15039\t15039\t50\t-\t.\tcoverage=169;context=CCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCAGG;IPDRatio=1.82 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15066\t15066\t37\t+\t.\tcoverage=161;context=CCCAGAGTCGCCGGGGCCAAGTCAGGTGGCGTATTCCAGAT;IPDRatio=1.69 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15074\t15074\t33\t-\t.\tcoverage=171;context=CCAGGACAATCTGGAATACGCCACCTGACTTGGCCCCGGCG;IPDRatio=1.80 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15078\t15078\t33\t-\t.\tcoverage=168;context=GCCCCCAGGACAATCTGGAATACGCCACCTGACTTGGCCCC;IPDRatio=1.62 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15095\t15095\t32\t-\t.\tcoverage=167;context=ATCCGGCAATGGCGGCAGCCCCCAGGACAATCTGGAATACG;IPDRatio=1.59 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15118\t15118\t32\t-\t.\tcoverage=163;context=GGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGGCAG;IPDRatio=1.54 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15121\t15121\t31\t-\t.\tcoverage=149;context=AAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGG;IPDRatio=1.58 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15124\t15124\t34\t-\t.\tcoverage=170;context=TGCAAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGG;IPDRatio=1.57 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15195\t15195\t33\t-\t.\tcoverage=164;context=AGCACCATACTGGCACCGAGAGAAAACAGGATGCCGGTCAT;IPDRatio=1.60 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15239\t15239\t33\t+\t.\tcoverage=157;context=TGGTGTGGCGCAGATGCTGGCACCGAAAGCCAGAACTCCCC;IPDRatio=1.65 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15261\t15261\t38\t-\t.\tcoverage=167;context=CCGTTATCCGTTGTCTGTATACGGGGAGTTCTGGCTTTCGG;IPDRatio=1.73 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15291\t15291\t31\t+\t.\tcoverage=158;context=ACGGATAACGGTAAGCAGAACACCTATTTCTCCTCACTGGA;IPDRatio=1.68 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15292\t15292\t41\t-\t.\tcoverage=166;context=ATCCAGTGAGGAGAAATAGGTGTTCTGCTTACCGTTATCCG;IPDRatio=1.61 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15311\t15311\t31\t-\t.\tcoverage=169;context=TGCCCTGGGCAACCATGTTATCCAGTGAGGAGAAATAGGTG;IPDRatio=1.47 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15402\t15402\t31\t-\t.\tcoverage=169;context=TGACCACCGTCCCCTTCGTCTGCCGTGCTGATCTCCTGAGA;IPDRatio=1.45 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15566\t15566\t32\t+\t.\tcoverage=118;context=GAAGGACAACCTGAAGTCCACGCAGTTGCTGAGTGTGATCG;IPDRatio=1.79 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t16035\t16035\t32\t-\t.\tcoverage=38;context=ATCCTGCGCATCCGGATATTAAACGGGCGCGGCGGCAGGTT;IPDRatio=1.96 (esc)

  $  linecount tmp.gff
  40

What about the H5 file?

  $ h5ls -r tmp.h5
  /                        Group
  /lambda_NEB3011          Dataset {48502}
