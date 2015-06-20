

Run base modification detection on B. subtilis P6 chemistry validation data

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Bsub_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary/Smrtanalysis/current/common/references
  $ export REF_SEQ=${REF_DIR}/B_subtilis_strW23/sequence/B_subtilis_strW23.fasta

  $ ipdSummary.py ${BAMFILE} --reference ${REF_SEQ} --gff tst_Bsub.gff --csv tst_Bsub.csv --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --seed 12345

  $ linecount tst_Bsub.csv
  8055333

This is an amplified control but it will still find some "modifications".

  $ linecount tst_Bsub.gff
  41460
  $ NLINES="`linecount tst_Bsub.gff`"
  $ python -c "assert 41450 < ${NLINES} < 41470, ${NLINES}"

Most of these can't be positively identified, however.

  $ grep -c m4C tst_Bsub.gff
  8390
  $ NM4="`grep -c m4C tst_Bsub.gff`"
  $ python -c "assert 8380 < ${NM4} < 8391, ${NM4}"

  $ grep -c m6A tst_Bsub.gff
  3278
  $ NM6="`grep -c m6A tst_Bsub.gff`"
  $ python -c "assert 3270 < ${NM6} < 3280, ${NM6}"

  $ grep -c modified_base tst_Bsub.gff
  29790
  $ N_OTHER="`grep -c modified_base tst_Bsub.gff`"
  $ python -c "assert 29780 < ${N_OTHER} < 29800, ${N_OTHER}"
