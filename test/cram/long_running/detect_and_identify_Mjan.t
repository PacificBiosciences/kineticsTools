
Run base modification detection on M. jannaschii P6 chemistry validation data.

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/pbi/dept/secondary/siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Mjan_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary/Smrtanalysis/current/common/references
  $ export REF_SEQ=${REF_DIR}/Methanocaldococcus_jannaschii_DSM2661/sequence/Methanocaldococcus_jannaschii_DSM2661.fasta

  $ ipdSummary ${BAMFILE} --reference ${REF_SEQ} --gff tst_Mjan.gff --csv tst_Mjan.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C

  $ linecount tst_Mjan.csv
  3479855

This one has relatively few modificiations, 60% of which are identifiable.

  $ linecount tst_Mjan.gff
  18654

  $ grep -c m4C tst_Mjan.gff
  3328

  $ grep -c m6A tst_Mjan.gff
  9506

  $ grep -c modified_base tst_Mjan.gff
  5816
