###ALign_and_Filtering
THR=8
IND=hg19_sacCer3/bwa_index/hg19_sacCer.fa #genome index prefix of hg19_saccer genome
BLACKLIST=Genome_Reference/hg19/hg19-blacklist.v2.bed
WORK_DIR=/Work/
cd $WORK_DIR
mkdir trim_galore
mkdir Align
############################################Fastqc - Trimming
nohup sh -c 'parallel trim_galore --illumina --paired --fastqc -o trim_galore/ {} {=s/_1/_2/=} ::: *_1.fastq.gz'&
###Rename files
for i in *_val_1.fq.gz
do
  a=$(basename $i _val_1.fq.gz)
  mv $i $a.trim.fq.gz
done
for i in *_val_2.fq.gz
do
  a=$(basename $i _val_2.fq.gz)
  mv $i $a.trim.fq.gz
done
############################################Alignment
#create filelist.txt with filenames
nohup sh -c '
THR=8
PICARD="Tools/picard.jar"
for i in $(cat filelist.txt)
do
  EXT=".trim.fq.gz"
  HEAD=`basename $i _1.fastq.gz`
  FOR=$HEAD"_1"
  REV=$HEAD"_2"
  bwa mem -t $THR $IND trim_galore/$FOR$EXT trim_galore/$REV$EXT |samtools view -bS -@ $THR - | samtools sort -@ $THR -m 5G - -o Align/$HEAD.pre.bam
  samtools index -@ $THR Align/$HEAD.pre.bam
  samtools view -h -@ $THR Align/$HEAD.pre.bam $(cat chr_list.txt) |samtools view -@ 8 -q 20 -L $BLACKLIST -U Align/$HEAD.dup.bam -o /dev/null -
  samtools index -@ $THR Align/$HEAD.dup.bam
  java -Xmx40g -jar $PICARD MarkDuplicates I="Align/$HEAD.dup.bam" O="Align/$HEAD.bam" M="Align/$HEAD.md.txt" AS=true REMOVE_DUPLICATES=true
  samtools index -@ $THR Align/$HEAD.bam
  samtools stats -@ $THR Align/$HEAD.pre.bam > Align/$HEAD.pre.stats
  samtools flagstat -@ $THR Align/$HEAD.pre.bam > Align/$HEAD.pre.flagstat
  samtools idxstats -@ $THR Align/$HEAD.pre.bam > Align/$HEAD.pre.idxstat
  samtools stats -@ $THR Align/$HEAD.bam > Align/$HEAD.stats
  samtools flagstat -@ $THR Align/$HEAD.bam > Align/$HEAD.flagstat
  samtools idxstats -@ $THR Align/$HEAD.bam > Align/$HEAD.idxstat
done
'&
###########QUALITY FILTER
for i in *.bam
do
SAMPLE=$(basename $i .bam)
samtools view -q 10 -h $i  | awk 'substr($0,1,1)=="@" || ($9>= 60 && $9<=900) || ($9<=-60 && $9>=-900)' | samtools view -b > $SAMPLE.filtered.bam
samtools index -@ 8 $SAMPLE.filtered.bam
done
nohup sh -c '
for i in *filtered.bam
do
SAMPLE=$(basename $i .bam)
samtools stats -@ 8 $SAMPLE.bam > $SAMPLE.stats
samtools flagstat -@ 8 $SAMPLE.bam > $SAMPLE.flagstat
samtools idxstats -@ 8  $SAMPLE.bam > $SAMPLE.idxstat
done
' &
############################################ Genome splitting
nohup sh -c '
hg19_genome=Genome_Reference/hg19_sacCer3/hg19.genome.bed
saccer3_genome=Genome_Reference/hg19_sacCer3/saccer3.genome.bed
touch all_counts
for f in *.bam ; do echo $f;echo $f >> all_counts &&  samtools view -c $f >> all_counts; done
for f in *.bam
do
  #hg19
  samtools view -b -@8 $f -L $hg19_genome > ${f%%.bam}.hg19.bam
  #saccer3
  samtools view -b -@8 $f -L $saccer3_genome > ${f%%.bam}.saccer3.bam
  done
touch all_saccer3_counts
for f in *saccer3.bam; do echo $f;echo $f >> all_saccer3_counts &&  samtools view -c $f >> all_saccer3_counts; done
touch all_hg19_counts
for f in *hg19.bam; do echo $f;echo $f >> all_hg19_counts &&  samtools view -c $f >> all_hg19_counts; done
'&
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_hg19_counts > all_hg19_counts.tab.txt
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_saccer3_counts > all_saccer3_counts.tab.txt
############################################ === recover stats ===
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_counts > all_counts.tab.txt
#####MACS peak calling on human
mkdir MACS
nohup sh -c '
macs2 callpeak -f BAMPE -t 5-H-IP_S5.Rep1.filtered.hg19.bam -c 5-IN_S6.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/5-H-IP.Rep1 -n 5-H-IP.Rep1
macs2 callpeak -f BAMPE -t 5-H-IP_S5.Rep2.filtered.hg19.bam -c 5-IN_S6.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/5-H-IP.Rep2 -n 5-H-IP.Rep2
macs2 callpeak -f BAMPE -t 60-H-IP_S8.Rep1.filtered.hg19.bam -c 60-IN_S9.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/60-H-IP.Rep1 -n 60-H-IP.Rep1
macs2 callpeak -f BAMPE -t 60-H-IP_S8.Rep2.filtered.hg19.bam -c 60-IN_S9.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/60-H-IP.Rep2 -n 60-H-IP.Rep2
macs2 callpeak -f BAMPE -t 5-IP_S4.Rep1.filtered.hg19.bam -c 5-IN_S6.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/5-IP.Rep1 -n 5-IP.Rep1
macs2 callpeak -f BAMPE -t 5-IP_S4.Rep2.filtered.hg19.bam -c 5-IN_S6.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/5-IP.Rep2 -n 5-IP.Rep2
macs2 callpeak -f BAMPE -t 60-IP_S7.Rep1.filtered.hg19.bam -c 60-IN_S9.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/60-IP.Rep1 -n 60-IP.Rep1
macs2 callpeak -f BAMPE -t 60-IP_S7.Rep2.filtered.hg19.bam -c 60-IN_S9.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/60-IP.Rep2 -n 60-IP.Rep2
macs2 callpeak -f BAMPE -t CT-H-IP_S2.Rep1.filtered.hg19.bam -c CT-IN_S3.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/CT-H-IP.Rep1 -n CT-H-IP.Rep1
macs2 callpeak -f BAMPE -t CT-H-IP_S2.Rep2.filtered.hg19.bam -c CT-IN_S3.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/CT-H-IP.Rep2 -n CT-H-IP.Rep2
macs2 callpeak -f BAMPE -t CT-IP_S1.Rep1.filtered.hg19.bam -c CT-IN_S3.Rep1.filtered.hg19.bam -g hs -B --outdir MACS/CT-IP.Rep1 -n CT-IP.Rep1
macs2 callpeak -f BAMPE -t CT-IP_S1.Rep2.filtered.hg19.bam -c CT-IN_S3.Rep2.filtered.hg19.bam -g hs -B --outdir MACS/CT-IP.Rep2 -n CT-IP.Rep2
'&
######StepbyStep
#crea bedgraph
mkdir Signal_tracks
nohup sh -c '
for i in *hg19.bam
do
SAMPLE=$(basename $i .filtered.hg19.bam)
macs2 pileup -i $i -f BAMPE -o Signal_tracks/$SAMPLE.bdg
done
'&
#moltiplica scalfactor
macs2 bdgopt -i $SAMPLE.bdg -m multiply -p 0.65200919 -o $SAMPLE.scaled.bdg





nohup sh -c '
macs2 bdgopt -i 5-H-IP_S5.Rep1.bdg -m multiply -p 0.65200919 -o 5-H-IP_S5.Rep1.scaled.bdg
macs2 bdgopt -i 5-H-IP_S5.Rep2.bdg -m multiply -p 0.160254485 -o 5-H-IP_S5.Rep2.scaled.bdg
macs2 bdgopt -i 5-IP_S4.Rep1.bdg -m multiply -p 1.170921637 -o 5-IP_S4.Rep1.scaled.bdg
macs2 bdgopt -i 5-IP_S4.Rep2.bdg -m multiply -p 0.438903194 -o 5-IP_S4.Rep2.scaled.bdg
macs2 bdgopt -i 5-IN_S6.Rep1.bdg -m multiply -p 0.434874033 -o 5-IN_S6.Rep1.scaled.bdg
macs2 bdgopt -i 5-IN_S6.Rep2.bdg -m multiply -p 0.285427429 -o 5-IN_S6.Rep2.scaled.bdg
macs2 bdgcmp -t 5-H-IP_S5.Rep1.scaled.bdg -c 5-IN_S6.Rep1.scaled.bdg -m qpois -o 5-H-IP_S5.Rep1.topeakcall.bdg
macs2 bdgcmp -t 5-H-IP_S5.Rep2.scaled.bdg -c 5-IN_S6.Rep2.scaled.bdg -m qpois -o 5-H-IP_S5.Rep2.topeakcall.bdg
macs2 bdgcmp -t 5-IP_S4.Rep1.scaled.bdg -c 5-IN_S6.Rep1.scaled.bdg -m qpois -o 5-IP_S4.Rep1.topeakcall.bdg
macs2 bdgcmp -t 5-IP_S4.Rep2.scaled.bdg -c 5-IN_S6.Rep2.scaled.bdg  -m qpois -o 5-IP_S4.Rep2.topeakcall.bdg
macs2 bdgbroadcall -i 5-H-IP_S5.Rep1.topeakcall.bdg -c 3 --outdir 5-H-IP_S5.Rep1 --o-prefix 5-H-IP_S5.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i 5-H-IP_S5.Rep2.topeakcall.bdg -c 3 --outdir 5-H-IP_S5.Rep2 --o-prefix 5-H-IP_S5.Rep2 -g 200 -G 1000
macs2 bdgbroadcall -i 5-IP_S4.Rep1.topeakcall.bdg -c 3 --outdir 5-IP_S4.Rep1 --o-prefix 5-IP_S4.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i 5-IP_S4.Rep2.topeakcall.bdg -c 3 --outdir 5-IP_S4.Rep2 --o-prefix 5-IP_S4.Rep2 -g 200 -G 1000
macs2 bdgpeakcall -i 5-H-IP_S5.Rep1.topeakcall.bdg -c 3 --outdir 5-H-IP_S5.Rep1 --o-prefix 5-H-IP_S5.Rep1
macs2 bdgpeakcall -i 5-H-IP_S5.Rep2.topeakcall.bdg -c 3 --outdir 5-H-IP_S5.Rep2 --o-prefix 5-H-IP_S5.Rep2
macs2 bdgpeakcall -i 5-IP_S4.Rep1.topeakcall.bdg -c 3 --outdir 5-IP_S4.Rep1 --o-prefix 5-IP_S4.Rep1
macs2 bdgpeakcall -i 5-IP_S4.Rep2.topeakcall.bdg -c 3 --outdir 5-IP_S4.Rep2 --o-prefix 5-IP_S4.Rep2
'&
nohup sh -c '
macs2 bdgopt -i 60-H-IP_S8.Rep1.bdg -m multiply -p 1.456870091 -o 60-H-IP_S8.Rep1.scaled.bdg
macs2 bdgopt -i 60-H-IP_S8.Rep2.bdg -m multiply -p 0.178424036 -o 60-H-IP_S8.Rep2.scaled.bdg
macs2 bdgopt -i 60-IP_S7.Rep1.bdg -m multiply -p 1 -o 60-IP_S7.Rep1.scaled.bdg
macs2 bdgopt -i 60-IP_S7.Rep2.bdg -m multiply -p 0.395534427 -o 60-IP_S7.Rep2.scaled.bdg
macs2 bdgopt -i 60-IN_S9.Rep1.bdg -m multiply -p 0.373576459 -o 60-IN_S9.Rep1.scaled.bdg
macs2 bdgopt -i 60-IN_S9.Rep2.bdg -m multiply -p 0.259185941 -o 60-IN_S9.Rep2.scaled.bdg
macs2 bdgcmp -t 60-H-IP_S8.Rep1.scaled.bdg -c 60-IN_S9.Rep1.scaled.bdg -m qpois -o 60-H-IP_S8.Rep1.topeakcall.bdg
macs2 bdgcmp -t 60-H-IP_S8.Rep2.scaled.bdg -c 60-IN_S9.Rep2.scaled.bdg -m qpois -o 60-H-IP_S8.Rep2.topeakcall.bdg
macs2 bdgcmp -t 60-IP_S7.Rep1.scaled.bdg -c 60-IN_S9.Rep1.scaled.bdg -m qpois -o 60-IP_S7.Rep1.topeakcall.bdg
macs2 bdgcmp -t 60-IP_S7.Rep2.scaled.bdg -c 60-IN_S9.Rep2.scaled.bdg  -m qpois -o 60-IP_S7.Rep2.topeakcall.bdg
macs2 bdgbroadcall -i 60-H-IP_S8.Rep1.topeakcall.bdg -c 3 --outdir 60-H-IP_S8.Rep1 --o-prefix 60-H-IP_S8.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i 60-H-IP_S8.Rep2.topeakcall.bdg -c 3 --outdir 60-H-IP_S8.Rep2 --o-prefix 60-H-IP_S8.Rep2 -g 200 -G 1000
macs2 bdgbroadcall -i 60-IP_S7.Rep1.topeakcall.bdg -c 3 --outdir 60-IP_S7.Rep1 --o-prefix 60-IP_S7.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i 60-IP_S7.Rep2.topeakcall.bdg -c 3 --outdir 60-IP_S7.Rep2 --o-prefix 60-IP_S7.Rep2 -g 200 -G 1000
macs2 bdgpeakcall -i 60-H-IP_S8.Rep1.topeakcall.bdg -c 3 --outdir 60-H-IP_S8.Rep1 --o-prefix 60-H-IP_S8.Rep1
macs2 bdgpeakcall -i 60-H-IP_S8.Rep2.topeakcall.bdg -c 3 --outdir 60-H-IP_S8.Rep2 --o-prefix 60-H-IP_S8.Rep2
macs2 bdgpeakcall -i 60-IP_S7.Rep1.topeakcall.bdg -c 3 --outdir 60-IP_S7.Rep1 --o-prefix 60-IP_S7.Rep1
macs2 bdgpeakcall -i 60-IP_S7.Rep2.topeakcall.bdg -c 3 --outdir 60-IP_S7.Rep2 --o-prefix 60-IP_S7.Rep2
'&
nohup sh -c '
macs2 bdgopt -i CT-H-IP_S2.Rep1.bdg -m multiply -p 0.439454171 -o CT-H-IP_S2.Rep1.scaled.bdg
macs2 bdgopt -i CT-H-IP_S2.Rep2.bdg -m multiply -p 0.102144009 -o CT-H-IP_S2.Rep2.scaled.bdg
macs2 bdgopt -i CT-IP_S1.Rep1.bdg -m multiply -p 0.677663401 -o CT-IP_S1.Rep1.scaled.bdg
macs2 bdgopt -i CT-IP_S1.Rep2.bdg -m multiply -p 0.39306033 -o CT-IP_S1.Rep2.scaled.bdg
macs2 bdgopt -i CT-IN_S3.Rep1.bdg -m multiply -p 0.617261648 -o CT-IN_S3.Rep1.scaled.bdg
macs2 bdgopt -i CT-IN_S3.Rep2.bdg -m multiply -p 0.183994298 -o CT-IN_S3.Rep2.scaled.bdg
macs2 bdgcmp -t CT-H-IP_S2.Rep1.scaled.bdg -c CT-IN_S3.Rep1.scaled.bdg -m qpois -o CT-H-IP_S2.Rep1.topeakcall.bdg
macs2 bdgcmp -t CT-H-IP_S2.Rep2.scaled.bdg -c CT-IN_S3.Rep2.scaled.bdg -m qpois -o CT-H-IP_S2.Rep2.topeakcall.bdg
macs2 bdgcmp -t CT-IP_S1.Rep1.scaled.bdg -c CT-IN_S3.Rep1.scaled.bdg -m qpois -o CT-IP_S1.Rep1.topeakcall.bdg
macs2 bdgcmp -t CT-IP_S1.Rep2.scaled.bdg -c CT-IN_S3.Rep2.scaled.bdg  -m qpois -o CT-IP_S1.Rep2.topeakcall.bdg
macs2 bdgbroadcall -i CT-H-IP_S2.Rep1.topeakcall.bdg -c 3 --outdir CT-H-IP_S2.Rep1 --o-prefix CT-H-IP_S2.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i CT-H-IP_S2.Rep2.topeakcall.bdg -c 3 --outdir CT-H-IP_S2.Rep2 --o-prefix CT-H-IP_S2.Rep2 -g 200 -G 1000
macs2 bdgbroadcall -i CT-IP_S1.Rep1.topeakcall.bdg -c 3 --outdir CT-IP_S1.Rep1 --o-prefix CT-IP_S1.Rep1 -g 200 -G 1000
macs2 bdgbroadcall -i CT-IP_S1.Rep2.topeakcall.bdg -c 3 --outdir CT-IP_S1.Rep2 --o-prefix CT-IP_S1.Rep2 -g 200 -G 1000
macs2 bdgpeakcall -i CT-H-IP_S2.Rep1.topeakcall.bdg -c 3 --outdir CT-H-IP_S2.Rep1 --o-prefix CT-H-IP_S2.Rep1
macs2 bdgpeakcall -i CT-H-IP_S2.Rep2.topeakcall.bdg -c 3 --outdir CT-H-IP_S2.Rep2 --o-prefix CT-H-IP_S2.Rep2
macs2 bdgpeakcall -i CT-IP_S1.Rep1.topeakcall.bdg -c 3 --outdir CT-IP_S1.Rep1 --o-prefix CT-IP_S1.Rep1
macs2 bdgpeakcall -i CT-IP_S1.Rep2.topeakcall.bdg -c 3 --outdir CT-IP_S1.Rep2 --o-prefix CT-IP_S1.Rep2
'&
###Conversion
for i in *scaled.bdg
do
  SAMPLE= $(basename $i .bdg)
  bedGraphtoBigwig $i hg19genome $SAMPLE.bw
done
