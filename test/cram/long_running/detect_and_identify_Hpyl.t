
Run base modification detection on H. pylori P6 chemistry validation data.

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Hpyl_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary-siv/references
  $ export REF_SEQ=${REF_DIR}/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

  $ ipdSummary.py ${BAMFILE} --reference ${REF_SEQ} --gff tst_Hpyl.gff --csv tst_Hpyl.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C

  $ linecount tst_Hpyl.csv
  3287635

This one also has lots of modifications, mostly m6A.

  $ linecount tst_Hpyl.gff
  88452

  $ grep -c m4C tst_Hpyl.gff
  13723

  $ grep -c m6A tst_Hpyl.gff
  59721

  $ grep -c modified_base tst_Hpyl.gff
  15006
