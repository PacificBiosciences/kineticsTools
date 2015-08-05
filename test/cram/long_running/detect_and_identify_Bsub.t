

Run base modification detection on B. subtilis P6 chemistry validation data

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Bsub_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary-siv/references
  $ export REF_SEQ=${REF_DIR}/B_subtilis_strW23/sequence/B_subtilis_strW23.fasta

  $ ipdSummary ${BAMFILE} --reference ${REF_SEQ} --gff tst_Bsub.gff --csv tst_Bsub.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C

  $ linecount tst_Bsub.csv
  8055333

This is an amplified control but it will still find some "modifications".

  $ linecount tst_Bsub.gff
  25513

Most of these can't be positively identified, however.

  $ grep -c m4C tst_Bsub.gff
  3141

  $ grep -c m6A tst_Bsub.gff
  573

  $ grep -c modified_base tst_Bsub.gff
  21797
