

Run base modification detection on C. aggregans P6 chemistry validation data

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Cagg_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary/Smrtanalysis/current/common/references
  $ export REF_SEQ=${REF_DIR}/Chloroflexus_aggregans_DSM9485/sequence/Chloroflexus_aggregans_DSM9485.fasta

  $ ipdSummary.py ${BAMFILE} --reference ${REF_SEQ} --gff tst_Cagg.gff --csv tst_Cagg.csv --numWorkers 12 --pvalue 0.001 --identify m6A,m4C

  $ linecount tst_Cagg.csv
  9369863

This one actually has real modifications (and lots of them).

  $ linecount tst_Cagg.gff
  188678
  $ grep -c m4C tst_Cagg.gff
  40055
  $ grep -c m6A tst_Cagg.gff
  138566
  $ grep -c modified_base tst_Cagg.gff
  10055
