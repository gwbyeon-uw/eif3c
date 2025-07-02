#!/bin/bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=4G

reads=$1
name=$2
outdir=$3
starIndex=index/star
gtf=index/mm10/gencode.vM18.annotation.gtf
win=index/mm10_all_w10_s5.bed
post="w10_s5"

mkdir -p $outdir
mkdir -p $outdir/fastqc
mkdir -p $outdir/STAR/$name
mkdir -p $outdir/STAR2/$name

#FASTQC
fastqc -t 24 --noextract -o $outdir/fastqc -q $reads

#Illumina adapter cutting
cutadapt -j 24 -m 18 -a CTGTAGGCACCATCAAT --overlap=5 --trimmed-only --quality-cutoff=20 -o $outdir/$name".trimmed.fastq.gz" $reads

#UMI tools; reformats fastq so that UMI is in the header
umi_tools extract -I $outdir/$name".trimmed.fastq.gz" --bc-pattern=NNNNXXXXXXNNNNN --stdout=$outdir/$name".umi.fastq.gz"

#clip 1 base from 5' end
cutadapt -j 24 -u 6 -o $outdir/$name".5pclip.fastq.gz" $outdir/$name".umi.fastq.gz"

#STAR alignment
STAR --runThreadN 24 --genomeDir $starIndex --readFilesIn $outdir/$name".5pclip.fastq.gz" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/STAR/$name/$name"." --outFilterMultimapNmax 1 --alignEndsType EndToEnd --alignIntronMax 1000000 --alignIntronMin 20 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterType BySJout --quantMode TranscriptomeSAM GeneCounts

#Index BAM
samtools index $outdir/STAR/$name/$name".Aligned.sortedByCoord.out.bam"

#UMI-tools deduplication
umi_tools dedup -I $outdir/STAR/$name/$name".Aligned.sortedByCoord.out.bam" -v 6 --buffer-whole-contig --output-stats=$outdir/$name".dedup" -S $outdir/$name".dedup.bam"

#Sort
samtools sort --threads 24 -n $outdir/$name".dedup.bam" > $outdir/$name".dedup.namesort.bam"

#Realign using UMI dedup'ed reads
STAR --runThreadN 24 --genomeDir $starIndex --readFilesType SAM SE --readFilesIn $outdir/$name".dedup.namesort.bam" --readFilesCommand samtools view --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/STAR2/$name/$name"." --outFilterMultimapNmax 1 --alignEndsType EndToEnd --alignIntronMax 1000000 --alignIntronMin 20 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterType BySJout --quantMode TranscriptomeSAM GeneCounts

#Coverage
bedtools coverage -counts -split -s -sorted -g index/star/chrNameLength.txt -b $outdir/STAR2/$name/$name".Aligned.sortedByCoord.out.bam" -a $win | gzip > $outdir/$name"_winCounts_"$post".gz"
