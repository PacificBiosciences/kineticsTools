Test basic mode of ipdSummary.

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=$TESTDIR/../data
  $ INPUT=$DATA/p4-c2-lambda-mod-decode.cmp.h5
  $ REFERENCE=$DATA/lambda/sequence/lambda.fasta

Run basic ipdSummary:

  $ ipdSummary --log-level=WARNING --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --csv tmp.csv --gff tmp.gff --summary_h5 tmp.h5 --reference $REFERENCE $INPUT

Look at output csv file:

  $ head -3 tmp.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "lambda_NEB3011",13190,0,A,1,0.872,0.230,1.125,0.776,3
  "lambda_NEB3011",13190,1,T,1,0.909,0.252,1.126,0.808,3


  $ linecount tmp.csv
  7603

Look at output gff file:

  $ cat tmp.gff
  ##gff-version 3
  ##source ipdSummary * (glob)
  ##source-commandline * (glob)
  ##sequence-region lambda_NEB3011 1 48502
  lambda_NEB3011\tkinModCall\tmodified_base\t14657\t14657\t34\t+\t.\tcoverage=155;context=CGGCACAGCCGGGCGATGTGCTGCTGTGCTGTTTTGGTTCA;IPDRatio=1.55 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14743\t14743\t34\t-\t.\tcoverage=173;context=TACCTCTCTCGTTTGCTCAGTTGTTCAGGAATATGGTGCAG;IPDRatio=1.54 (esc)
  lambda_NEB3011\tkinModCall\tm4C\t14756\t14756\t21\t-\t.\tcoverage=165;context=CCATTTGTCGGTGTACCTCTCTCGTTTGCTCAGTTGTTCAG;IPDRatio=1.32;identificationQv=24 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14769\t14769\t32\t-\t.\tcoverage=168;context=GTGTGCGTCGCTGCCATTTGTCGGTGTACCTCTCTCGTTTG;IPDRatio=1.56 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14800\t14800\t32\t-\t.\tcoverage=173;context=GCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTGCGTCG;IPDRatio=1.66 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14805\t14805\t31\t-\t.\tcoverage=167;context=CAGATGCGCGCCATGCCCGGTGACGCCAGAGGGAGTGTGTG;IPDRatio=1.64 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14831\t14831\t31\t+\t.\tcoverage=161;context=CATGGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTG;IPDRatio=1.52 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14834\t14834\t32\t+\t.\tcoverage=166;context=GGCGCGCATCTGCCTTTACGGGGATTTACAACGATTTGGTC;IPDRatio=1.71 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14841\t14841\t32\t+\t.\tcoverage=162;context=ATCTGCCTTTACGGGGATTTACAACGATTTGGTCGCCGCAT;IPDRatio=1.58 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14847\t14847\t36\t-\t.\tcoverage=172;context=AGGTCGATGCGGCGACCAAATCGTTGTAAATCCCCGTAAAG;IPDRatio=1.85 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14864\t14864\t45\t-\t.\tcoverage=166;context=CCCCCGTTTTCACACGAAGGTCGATGCGGCGACCAAATCGT;IPDRatio=1.70 (esc)
  lambda_NEB3011\tkinModCall\tm4C\t14884\t14884\t30\t-\t.\tcoverage=173;context=CAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGGT;IPDRatio=2.10;identificationQv=19 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14885\t14885\t33\t-\t.\tcoverage=166;context=CCAGTGCCCGGATGGCTTCAGCCCCCGTTTTCACACGAAGG;IPDRatio=2.52 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t14909\t14909\t35\t+\t.\tcoverage=162;context=AGCCATCCGGGCACTGGCCACACAGCTCCCGGCGTTTCGTC;IPDRatio=1.76 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t14983\t14983\t208\t+\t.\tcoverage=160;context=TTGCCGGGCGGGACGTCAGCACGTCCGGGTTAACGGCGCAG;IPDRatio=6.71;identificationQv=191 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t14992\t14992\t207\t-\t.\tcoverage=162;context=CTCATGTAACTGCGCCGTTAACCCGGACGTGCTGACGTCCC;IPDRatio=6.56;identificationQv=186 (esc)
  lambda_NEB3011\tkinModCall\tm4C\t14997\t14997\t62\t-\t.\tcoverage=140;context=AGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGAC;IPDRatio=2.32;identificationQv=11 (esc)
  lambda_NEB3011\tkinModCall\tm4C\t14998\t14998\t22\t-\t.\tcoverage=172;context=CAGAGTCTCATGTAACTGCGCCGTTAACCCGGACGTGCTGA;IPDRatio=1.51;identificationQv=41 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15006\t15006\t40\t-\t.\tcoverage=170;context=CCATCAGGCAGAGTCTCATGTAACTGCGCCGTTAACCCGGA;IPDRatio=1.75 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15008\t15008\t47\t+\t.\tcoverage=162;context=CGGGTTAACGGCGCAGTTACATGAGACTCTGCCTGATGGCG;IPDRatio=1.73 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15039\t15039\t51\t-\t.\tcoverage=169;context=CCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCAGG;IPDRatio=1.83 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t15041\t15041\t29\t-\t.\tcoverage=161;context=CCCCGGCGACTCTGGGAACAATATGAATTACAGCGCCATCA;IPDRatio=1.54;identificationQv=3 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15066\t15066\t38\t+\t.\tcoverage=161;context=CCCAGAGTCGCCGGGGCCAAGTCAGGTGGCGTATTCCAGAT;IPDRatio=1.70 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15074\t15074\t33\t-\t.\tcoverage=171;context=CCAGGACAATCTGGAATACGCCACCTGACTTGGCCCCGGCG;IPDRatio=1.81 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15078\t15078\t33\t-\t.\tcoverage=168;context=GCCCCCAGGACAATCTGGAATACGCCACCTGACTTGGCCCC;IPDRatio=1.62 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15089\t15089\t32\t-\t.\tcoverage=167;context=CAATGGCGGCAGCCCCCAGGACAATCTGGAATACGCCACCT;IPDRatio=1.57 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15095\t15095\t31\t-\t.\tcoverage=167;context=ATCCGGCAATGGCGGCAGCCCCCAGGACAATCTGGAATACG;IPDRatio=1.57 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15118\t15118\t32\t-\t.\tcoverage=163;context=GGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGGCGGCAG;IPDRatio=1.55 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15124\t15124\t35\t-\t.\tcoverage=170;context=TGCAAGGGTGGCTCCGGCGGTAAAGAATGATCCGGCAATGG;IPDRatio=1.57 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15195\t15195\t32\t-\t.\tcoverage=164;context=AGCACCATACTGGCACCGAGAGAAAACAGGATGCCGGTCAT;IPDRatio=1.59 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15239\t15239\t33\t+\t.\tcoverage=157;context=TGGTGTGGCGCAGATGCTGGCACCGAAAGCCAGAACTCCCC;IPDRatio=1.65 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15261\t15261\t39\t-\t.\tcoverage=167;context=CCGTTATCCGTTGTCTGTATACGGGGAGTTCTGGCTTTCGG;IPDRatio=1.73 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15292\t15292\t39\t-\t.\tcoverage=166;context=ATCCAGTGAGGAGAAATAGGTGTTCTGCTTACCGTTATCCG;IPDRatio=1.57 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t15381\t15381\t20\t-\t.\tcoverage=169;context=GCCGTGCTGATCTCCTGAGAAACCACGCGTGACCCCACGCG;IPDRatio=1.40;identificationQv=13 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15402\t15402\t32\t-\t.\tcoverage=169;context=TGACCACCGTCCCCTTCGTCTGCCGTGCTGATCTCCTGAGA;IPDRatio=1.47 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t15566\t15566\t32\t+\t.\tcoverage=118;context=GAAGGACAACCTGAAGTCCACGCAGTTGCTGAGTGTGATCG;IPDRatio=1.79 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t15704\t15704\t21\t-\t.\tcoverage=84;context=CCTGCTCACCAGCCCGGAACACCACCGTGACACCGGATATG;IPDRatio=1.47;identificationQv=3 (esc)
  lambda_NEB3011\tkinModCall\tmodified_base\t16035\t16035\t31\t-\t.\tcoverage=38;context=ATCCTGCGCATCCGGATATTAAACGGGCGCGGCGGCAGGTT;IPDRatio=1.93 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t16380\t16380\t21\t-\t.\tcoverage=10;context=CATTTATCCACATCCGCCGCACCAAGACGTTTCCCCATGCC;IPDRatio=6.03;identificationQv=5 (esc)
  lambda_NEB3011\tkinModCall\tm6A\t16658\t16658\t21\t-\t.\tcoverage=6;context=GGGCGCTGAAGCTGTAGCGGAACGGCGCGCCATCATCCGGC;IPDRatio=2.78;identificationQv=20 (esc)


What about the H5 file?

  $ h5ls -r tmp.h5
  /                        Group
  /lambda_NEB3011          Dataset {48502}
