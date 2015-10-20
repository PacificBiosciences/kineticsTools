=========================================================
README for /pbi/dept/secondary/siv/testdata/kineticsTools
=========================================================

Most of these files are derived from Tyson Clark's P6 chemistry validation
experiments.  Bsub is an amplified control.

Bsub:
/mnt/data3/vol56/2530923/0003
/mnt/data3/vol56/2530923/0004

Cagg:
/mnt/data3/vol56/2530926/0003
/mnt/data3/vol56/2530926/0004
/mnt/data3/vol56/2530928/0003
/mnt/data3/vol56/2530928/0004

Hpyl:
/mnt/data3/vol56/2530926/0007
/mnt/data3/vol56/2530926/0008
/mnt/data3/vol56/2530928/0006
/mnt/data3/vol56/2530928/0005

Mjan:
/mnt/data3/vol56/2530924/0007
/mnt/data3/vol56/2530924/0008
/mnt/data3/vol56/2530928/0001
/mnt/data3/vol56/2530928/0002

Method to generate alignment files (python-like pseudo-code):
  bam_files = []
  for i, file_name in enumerate(h5_files):
    call("bax2bam subreads.1.bax.h5 -o unaligned.1")
    call("pbalign unaligned." + str(i) + ".subreads.bam /ref/seq/dir aligned." + str(i) + "1.bam --nproc=12 --seed=1 --minAccuracy=0.75 --minLength=50 --concordant --algorithmOptions=' -minMatch 12 -bestn 10 -minPctIdentity 70.0 -useQuality -minRawSubreadScore 800'")
    bam_files.append("unaligned." + str(i) + ".bam")
  call("~nechols/bin/bamtools merge -out aligned3.bam -in " +
       " -in ".join(bam_files))
  call("samtools index aligned3.bam")

Running ipdSummary on these inputs takes 1-3 hours on our cluster.  This
can be significantly reduced by using parallelization, but this introduces
stochastic behavior for some jobs.
