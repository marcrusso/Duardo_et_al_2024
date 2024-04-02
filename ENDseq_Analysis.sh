### FastQC of reads

$ for i in *.fastq.gz  
do  
fastqc -o ./FASTQC $i -t 8  
done  
#--------------------------------

### Trimming of adapters

$ for i in `cat samples.txt.txt` 
do 
SAMPLE=`basename $i _R1.fastq.gz` 
PAIR1=$SAMPLE"_R1.fastq.gz" 
PAIR2=$SAMPLE"_R2.fastq.gz" 
trim_galore $PAIR1 $PAIR2 --paired --illumina --fastqc -o ../TRIMMED/ 
done 

#--------------------------------

### Alignment

THR=8 
for i in `cat samples_trimmed.txt ` 
do  
SAMPLE=`basename $i _R1_val_1.fq.gz` 
PAIR1=$SAMPLE"_R1_val_1" 
PAIR2=$SAMPLE"_R2_val_2" 
echo $SAMPLE 
bowtie2 -p $THR -q --end-to-end -x ../Bowtie2_index/hg19/hg19.fa -1 $PAIR1 -2 $PAIR2 | samtools view -bS --threads 8 | samtools sort --threads 8 -o ../BAM_files_hg19/$SAMPLE".sorted.bam"  
done 

#--------------------------------

### Properly paired reads

for i in *.bam  
do 
SAMPLE=$(basename $i .sorted.bam) 
samtools view  -f 2 -F 512 -b -o ./PP_BAM_mm10/$SAMPLE"_pp".bam $i 
done 

#--------------------------------

### Filtering of aligned reads by MAPQ>30 

parallel --jobs 4 samtools view {} -q 30 -o ./MAPQ_filtered_BAM/{} ::: *sorted.bam
#--------------------------------

### Split properly paired and quality filtered reads in forward (F1R2) and reverse (F2R1), according to first in pair orientation

for i in *_pp.bam  
do 
SAMPLE=$(basename $i .bam) 
FLAG_FW1=$(echo _fwd99) 
FLAG_FW2=$(echo _fwd147) 
FLAG_REV1=$(echo _rev83) 
FLAG_REV2=$(echo _rev163) 
echo $i 
samtools view -f 99 $i -o $SAMPLE$FLAG_FW1.bam 
echo $SAMPLE$FLAG_FW1 
samtools view -f 147 $i -o $SAMPLE$FLAG_FW2.bam
echo $SAMPLE$FLAG_FW2 
samtools view -f 83 $i -o $SAMPLE$FLAG_REV1.bam 
echo $SAMPLE$FLAG_REV1 
samtools view -f 163 $i -o $SAMPLE$FLAG_REV2.bam 
echo $SAMPLE$FLAG_REV2 
samtools merge -f merged_$SAMPLE"_fwd.bam" $SAMPLE$FLAG_FW1.bam $SAMPLE$FLAG_FW2.bam 
samtools merge -f merged_$SAMPLE"_rev.bam" $SAMPLE$FLAG_REV1.bam $SAMPLE$FLAG_REV2.bam 
done 
# --------------------------------

### Peak calling
#double end DSB
parallel --jobs 4 "macs2 callpeak -t {} -c CTRL_a_pp.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ../../../Peak_calling_hg19_PP/no_lambda_PP/" ::: *0_a_pp.bam 
parallel --jobs 4 "macs2 callpeak -t {} -c CTRL_a_pp.bam -g 2652783500 --nomodel --nolambda -q 0.05 --keep-dup all -f BAMPE -n {} -B --outdir ../../../Peak_calling_mm10_PP/no_lambda_PP/" ::: *0_a_pp.bam
parallel --jobs 4 "macs2 callpeak -t {} -c CTRL_b_pp.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ../../../Peak_calling_hg19_PP/no_lambda_PP/" ::: *0_b_pp.bam 
parallel --jobs 4 "macs2 callpeak -t {} -c CTRL_b_pp.bam -g 2652783500 --nomodel --nolambda -q 0.05 --keep-dup all -f BAMPE -n {} -B --outdir ../../../Peak_calling_mm10_PP/no_lambda_PP/" ::: *0_b_pp.bam

#single end DSB
parallel --jobs 4 "macs2 callpeak -t {} -c merged_CTRL_a_pp_fwd.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ./Peaks_fwd.rev/peaks.fwd/" ::: *0_a_pp_fwd.bam
parallel --jobs 4 "macs2 callpeak -t {} -c merged_CTRL_a_pp_rev.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ./Peaks_fwd.rev/peaks.rev/" ::: *0_a_pp_rev.bam
parallel --jobs 4 "macs2 callpeak -t {} -c merged_CTRL_b_pp_fwd.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ./Peaks_fwd.rev/peaks.fwd/" ::: *0_b_pp_fwd.bam
parallel --jobs 4 "macs2 callpeak -t {} -c merged_CTRL_b_pp_rev.bam -g 2864785220 --nomodel --nolambda -q 0.05 --keep-dup all -n {} -f BAMPE -B --outdir ./Peaks_fwd.rev/peaks.rev/" ::: *0_b_pp_rev.bam




# BW with orientation 
PERCORSO_BAM="/mnt/d/END-Seq/BAM_files_hg19/MAPQ_filtered_BAM/PP_BAM_MAPQ/FW_REV/Merged_fw.rev/"
BLACKLIST="/mnt/d/END-Seq/"
OUTPUT="/mnt/d/END-Seq/BAM_files_hg19/MAPQ_filtered_BAM/PP_BAM_MAPQ/Bigwig_new"
bamCoverage --bam $PERCORSO_BAM/merged_CPT10_a_pp_fwd.bam -p 8 --scaleFactor 0.427412908 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT10.a.fw.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT10_a_pp_rev.bam -p 8 --scaleFactor 0.427412908 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT10.a.rev.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT20_a_pp_fwd.bam -p 8 --scaleFactor 0.702165267 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT20.a.fw.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT20_a_pp_rev.bam -p 8 --scaleFactor 0.702165267 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT20.a.rev.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CTRL_a_pp_fwd.bam -p 8 --scaleFactor 0.606158929 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CTRL.a.fwd.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CTRL_a_pp_rev.bam -p 8 --scaleFactor 0.606158929 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CTRL.a.rev.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT10_b_pp_fwd.bam -p 8 --scaleFactor 1 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT10.b.fwd.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT10_b_pp_rev.bam -p 8 --scaleFactor 1 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT10.b.rev.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT20_b_pp_fwd.bam -p 8 --scaleFactor 0.611025701 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT20.b.fwd.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CPT20_b_pp_rev.bam -p 8 --scaleFactor 0.611025701 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CPT20.b.rev.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CTRL_b_pp_fwd.bam -p 8 --scaleFactor 0.509311971 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CTRL.b.fwd.new.bw
bamCoverage --bam $PERCORSO_BAM/merged_CTRL_b_pp_rev.bam -p 8 --scaleFactor 0.509311971 -bl $BLACKLIST/hg19-blacklist.v2.bed -bs 10 -o $OUTPUT/CTRL.b.rev.new.bw
 


