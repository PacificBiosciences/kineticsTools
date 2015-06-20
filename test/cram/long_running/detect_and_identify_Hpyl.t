
Run base modification detection on H. pylori P6 chemistry validation data.

  $ . $TESTDIR/../portability.sh

  $ export DATA_DIR=/mnt/secondary-siv/testdata/kineticsTools
  $ export BAMFILE=${DATA_DIR}/Hpyl_aligned.subreads.bam
  $ export REF_DIR=/mnt/secondary/Smrtanalysis/current/common/references
  $ export REF_SEQ=${REF_DIR}/Helicobacter_pylori_J99/sequence/Helicobacter_pylori_J99.fasta

  $ ipdSummary.py ${BAMFILE} --reference ${REF_SEQ} --gff tst_Hpyl.gff --csv tst_Hpyl.csv --numWorkers 1 --pvalue 0.001 --identify m6A,m4C --seed 12345

  $ linecount tst_Hpyl.csv
  3287635

This one also has lots of modifications, mostly m6A.

  $ linecount tst_Hpyl.gff
  88465
  $ NLINES="`linecount tst_Hpyl.gff`"
  $ python -c "assert 88440 < ${NLINES} < 88530, ${NLINES}"

  $ grep -c m4C tst_Hpyl.gff
  13721
  $ NM4="`grep -c m4C tst_Hpyl.gff`"
  $ python -c "assert 13715 < ${NM4} < 13755, ${NM4}"

  $ grep -c m6A tst_Hpyl.gff
  59745
  $ NM6="`grep -c m6A tst_Hpyl.gff`"
  $ python -c "assert 59710 < ${NM6} < 59750, ${NM6}"

  $ grep -c modified_base tst_Hpyl.gff
  14997
  $ N_OTHER="`grep -c modified_base tst_Hpyl.gff`"
  $ python -c "assert 14950 < ${N_OTHER} < 15060, ${N_OTHER}"
