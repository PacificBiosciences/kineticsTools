Test detection and identification modes of ipdSummary using .xml dataset file as input.  Results should be identical to those using the equivalent .bam file.  This will also be tested for a split dataset (two non-overlapping .bam files in one .xml).

  $ . $TESTDIR/portability.sh

Load in data:

  $ DATA=/pbi/dept/secondary/siv/testdata/kineticsTools
  $ INPUT=$DATA/Hpyl_1_5000.xml
  $ REFERENCE=/pbi/dept/secondary/siv/references/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

Run basic ipdSummary:

  $ ipdSummary --log-level=WARNING --outfile tmp_xml1 --numWorkers 4 --pvalue 0.001 --identify m6A,m4C --reference $REFERENCE --referenceWindows="gi|12057207|gb|AE001439.1|:0-5000" --useChemistry P6-C4 $INPUT

Look at output csv file:

  $ head -3 tmp_xml1.csv
  refName,tpl,strand,base,score,tMean,tErr,modelPrediction,ipdRatio,coverage
  "gi|12057207|gb|AE001439.1|",1,0,A,10,2.387,0.464,1.710,1.396,29
  "gi|12057207|gb|AE001439.1|",1,1,T,1,0.492,0.075,0.602,0.817,57

  $ linecount tmp_xml1.csv
  10001

Look at output gff file:

  $ linecount tmp_xml1.gff
  274
  $ cat tmp_xml1.gff | head -20
  ##gff-version 3
  ##source ipdSummary v2.0
  ##source-commandline * (glob)
  ##sequence-region gi|12057207|gb|AE001439.1| 1 1643831
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t35\t35\t187\t-\t.\tcoverage=118;context=TTTAAGGGCGTTTTATGCCTAAATTTAAAAAATGATGCTGT;IPDRatio=5.68;identificationQv=196 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t60\t60\t49\t-\t.\tcoverage=112;context=AAAAAGCTCGCTCAAAAACCCTTGATTTAAGGGCGTTTTAT;IPDRatio=2.58;identificationQv=33 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t89\t89\t223\t+\t.\tcoverage=139;context=AGCGAGCTTTTTGCTCAAAGAATCCAAGATAGCGTTTAAAA;IPDRatio=5.69;identificationQv=187 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t91\t91\t217\t-\t.\tcoverage=143;context=ATTTTTAAACGCTATCTTGGATTCTTTGAGCAAAAAGCTCG;IPDRatio=6.34;identificationQv=214 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t113\t113\t41\t+\t.\tcoverage=132;context=CAAGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGT;IPDRatio=1.69 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tmodified_base\t115\t115\t33\t+\t.\tcoverage=147;context=AGATAGCGTTTAAAAATTTAGGGGTGTTAGGCTCAGCGTAG;IPDRatio=1.88 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t122\t122\t222\t-\t.\tcoverage=158;context=GCAAACTCTACGCTGAGCCTAACACCCCTAAATTTTTAAAC;IPDRatio=6.51;identificationQv=204 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t232\t232\t221\t+\t.\tcoverage=173;context=AGCGTAAAATCGCCTTTTCCATGCTCCCTAATCGCTTGAAA;IPDRatio=5.90;identificationQv=209 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t233\t233\t282\t-\t.\tcoverage=183;context=ATTTCAAGCGATTAGGGAGCATGGAAAAGGCGATTTTACGC;IPDRatio=6.43;identificationQv=262 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t241\t241\t267\t+\t.\tcoverage=178;context=TCGCCTTTTCCATGCTCCCTAATCGCTTGAAATCCCAGTCT;IPDRatio=5.57;identificationQv=234 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t248\t248\t239\t-\t.\tcoverage=185;context=ATTTAAAAGACTGGGATTTCAAGCGATTAGGGAGCATGGAA;IPDRatio=6.25;identificationQv=223 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t274\t274\t229\t-\t.\tcoverage=190;context=TGAGATTGACGCTCTCATCGAACCGCATTTAAAAGACTGGG;IPDRatio=6.83;identificationQv=220 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t277\t277\t272\t+\t.\tcoverage=188;context=AGTCTTTTAAATGCGGTTCGATGAGAGCGTCAATCTCATTG;IPDRatio=7.50;identificationQv=257 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm4C\t312\t312\t37\t-\t.\tcoverage=204;context=GCTTTAAGCCTTTTTAATGGCGTGTTAGAAAAAATCAATGA;IPDRatio=1.88;identificationQv=3 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t373\t373\t393\t+\t.\tcoverage=219;context=TAATCTTTTTTTCTTCTAACATGCTGGAAGCGATTTTTTTA;IPDRatio=7.11;identificationQv=353 (esc)
  gi|12057207|gb|AE001439.1|\tkinModCall\tm6A\t374\t374\t337\t-\t.\tcoverage=221;context=TTAAAAAAATCGCTTCCAGCATGTTAGAAGAAAAAAAGATT;IPDRatio=6.06;identificationQv=323 (esc)

Now try with a split dataset:

  $ INPUT=$DATA/Hpyl_1_5000_split.xml
  $ ipdSummary --log-level=WARNING --gff tmp_xml2.gff --csv tmp_xml2.csv --numWorkers 4 --pvalue 0.001 --identify m6A,m4C --reference $REFERENCE --referenceWindows="gi|12057207|gb|AE001439.1|:0-5000" --useChemistry P6-C4 $INPUT
  $ linecount tmp_xml2.gff
  274
