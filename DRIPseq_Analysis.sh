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
  
for f in *.bam
do
java -Xmx10g -jar $PICARD MarkDuplicates I=$f O=${f%%.bam}.nodup.bam M=${f%%.bam}.md.txt AS=true REMOVE_DUPLICATES=true
  echo "==========="
  echo "          "
done
touch all_saccer3_counts
for f in *saccer3*nodup.bam; do echo $f;echo $f >> all_saccer3_counts &&  samtools view -c $f >> all_saccer3_counts; done
touch all_hg19_counts
for f in *hg19*nodup.bam; do echo $f;echo $f >> all_hg19_counts &&  samtools view -c $f >> all_hg19_counts; done
touch all_counts
for f in *.bam ; do echo $f;echo $f >> all_counts &&  samtools view -c $f >> all_counts; done
awk '{printf "%s\t%s",$0,(NR%2?FS:RS)}' all_counts.txt > all_couns.tab.txt'
&
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

###Create a consensus file of bed peaks.
#To perform differential analysis
#deeptools bigwig with scale factors
BLACKLIST=Genome_Reference/hg19/hg19-blacklist.v2.bed
bamCoverage -b 5-H-IP_S5.hg19.nodup.bam -o 5-H-IP_S5.hg19.nodup.bw --scaleFactor 42.30 -p max -bl $BLACKLIST
bamCoverage -b 5-IN_S6.hg19.nodup.bam -o 5-IN_S6.hg19.nodup.bw --scaleFactor 9.91 -p max -bl $BLACKLIST
bamCoverage -b 5-IP_S4.hg19.nodup.bam -o 5-IP_S4.hg19.nodup.bw --scaleFactor 24.18 -p max -bl $BLACKLIST
bamCoverage -b 60-H-IP_S8.hg19.nodup.bam -o 60-H-IP_S8.hg19.nodup.bw --scaleFactor 21 -p max -bl $BLACKLIST
bamCoverage -b 60-IN_S9.hg19.nodup.bam -o 60-IN_S9.hg19.nodup.bw --scaleFactor 11.44 -p max -bl $BLACKLIST
bamCoverage -b 60-IP_S7.hg19.nodup.bam -o 60-IP_S7.hg19.nodup.bw --scaleFactor 28.74 -p max -bl $BLACKLIST
bamCoverage -b CT-H-IP_S2.hg19.nodup.bam -o CT-H-IP_S2.hg19.nodup.bw --scaleFactor 62.32 -p max -bl $BLACKLIST
bamCoverage -b CT-IN_S3.hg19.nodup.bam -o CT-IN_S3.hg19.nodup.bw --scaleFactor 6.73 -p max -bl $BLACKLIST
bamCoverage -b CT-IP_S1.hg19.nodup.bam -o CT-IP_S1.hg19.nodup.bw --scaleFactor 41.55 -p max -bl $BLACKLIST
######DROPA
cd Tools/DROPA
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Gain_Fast_Stable -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Gain_Fast_Stable.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Gain_Fast_Trans -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Gain_Fast_Trans.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Gain_Slow -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Gain_Slow.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Loss_Fast_Stable -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Loss_Fast_Stable.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Loss_Fast_Trans -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Loss_Fast_Trans.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o Loss_Slow -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome Loss_Slow.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o NoChange -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome NoChange.sorted.bed
python3 DROPA_v1.0.0.py -ref GeneReference/hg19_Ensembl/ -o AllPeaks -ex GROtoDROPA.txt -shuffle 100 -gsize GeneReference/hg19.genome All.Peaks.sorted.bed
###Metaplots
####DRIP_plot centered on peak center 
BLACKLIST=Genome_Reference/hg19/hg19-blacklist.v2.bed
SAMPLELIST=$(ls Deeptools/*IP.bw)
Sample_bed=$(ls  Deeptools/*bed)
computeMatrix reference-point -S $SAMPLELIST -R $Sample_bed -o out.mtx -bl $BLACKLIST -bs 50 -b 3000 -a 3000 --smartLabels -p max --missingDataAsZero --referencePoint center
plotHeatmap -m out.mtx -out DRIP.pdf --plotType se --perGroup --colorList 'white,darkred' 
##################
####DRP_plot centered on RSS
BLACKLIST=Genome_Reference/hg19/hg19-blacklist.v2.bed
SAMPLELIST=$(ls Deeptools/*IP.bw)
Sample_bed=$(ls  Deeptools/*bed)
computeMatrix reference-point -S $SAMPLELIST -R $Sample_bed -o out.mtx -bl $BLACKLIST -bs 100 -b 3000 -a 8000 --smartLabels -p max --missingDataAsZero --referencePoint center
plotHeatmap -m out.mtx -out Peaks_DRIP.pdf --plotType se --perGroup  --colorList 'white,darkred' --outFileSortedRegions tmp.bed
##################Plots for ChIPseq data
for i in Top1seq pol2 GROseq
do
SAMPLELIST=$(ls Deeptools/$i/*.bw)
computeMatrix reference-point -S $SAMPLELIST -R tmp.bed -o out.mtx -bl $BLACKLIST -bs 100 -b 3000 -a 8000 --smartLabels -p max --missingDataAsZero --referencePoint center
plotHeatmap -m out.mtx -out Peaks.$i.pdf --plotType se --perGroup --colorList 'white,darkred' --sortRegions no
done
###Plot on genes
for i in DRIP
do
SAMPLELIST=$(ls Deeptools/$i/*.bw)
REGIONS=$(ls *.gtf)
computeMatrix scale-regions -S  $SAMPLELIST -R $REGIONS -o out.mtx -bl $BLACKLIST -bs 50 -b 1500 -a 1500 -m 6000 --smartLabels -p max --missingDataAsZero --transcriptID gene --unscaled5prime 1000   
plotHeatmap -m out.mtx.gz -out Allgenes.Rloop.anterior.pdf --linesAtTickMarks --plotType se --perGroup --heatmapHeight 5 --colorList 'white,darkred' --outFileSortedRegions genes_sorted.bed
done 
for i in  pol2 GROseq END
do
SAMPLELIST=$(ls Baranello/$i/*.bw)
computeMatrix scale-regions -S  $SAMPLELIST -R genes_sorted.bed -o out.mtx -bl $BLACKLIST -bs 50 -b 1500 -a 1500 -m 6000 --smartLabels -p max --missingDataAsZero --transcriptID gene --unscaled5prime 1000
plotHeatmap -m out.mtx -out Allgenes.$i.anterior.pdf --linesAtTickMarks --plotType se --perGroup --heatmapHeight 5 --colorList 'white,darkgreen' --sortRegions no
done 
