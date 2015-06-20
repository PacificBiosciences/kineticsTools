
Run base modification detection on M. jannaschii P6 chemistry validation data.

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Mjan_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary/Smrtanalysis/current/common/references
  $ export REF_SEQ=${REF_DIR}/Methanocaldococcus_jannaschii_DSM2661/sequence/Methanocaldococcus_jannaschii_DSM2661.fasta

  $ ipdSummary.py ${BAMFILE} --reference ${REF_SEQ} --gff tst_Mjan.gff --csv tst_Mjan.csv --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --seed 12345

  $ linecount tst_Mjan.csv
  3479855

This one has relatively few modificiations, 60% of which are identifiable.

  $ linecount tst_Mjan.gff
  26673
  $ NLINES="`linecount tst_Mjan.gff`"
  $ python -c "assert 26660 < ${NLINES} < 26675, ${NLINES}"

  $ grep -c m4C tst_Mjan.gff
  5277
  $ NM4="`grep -c m4C tst_Mjan.gff`"
  $ python -c "assert 5270 < ${NM4} < 5280, ${NM4}"

  $ grep -c m6A tst_Mjan.gff
  11568
  $ NM6="`grep -c m6A tst_Mjan.gff`"
  $ python -c "assert 11565 < ${NM6} < 11580, ${NM6}"

  $ grep -c modified_base tst_Mjan.gff
  9824
  $ N_OTHER="`grep -c modified_base tst_Mjan.gff`"
  $ python -c "assert 9805 < ${N_OTHER} < 9825, ${N_OTHER}"
