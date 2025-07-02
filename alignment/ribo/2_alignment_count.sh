#!/bin/bash
#SBATCH -n 24
#SBATCH --mem-per-cpu=4G

reads=$1
name=$2
outdir=$3
starIndex=index/star
bed=index/vm18_ens93_regions_sorted.bed

#Illumina adapter cutting
cutadapt -j 24 -m 18 -a CTGTAGGCACCATCAAT --quality-cutoff=20 -o $outdir/$name".trimmed.fastq.gz" $reads

#STAR
STAR --runThreadN 24 --genomeDir $starIndex --readFilesIn $outdir/$name".trimmed.fastq.gz" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $outdir/STAR/$name/$name"." --outFilterMultimapNmax 1 --alignEndsType EndToEnd --alignIntronMax 1000000 --alignIntronMin 20 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterType BySJout --quantMode TranscriptomeSAM GeneCounts

#Sort
samtools sort --threads 24 -o $outdir/$name".tx.bam" $outdir/STAR/$name/$name".Aligned.toTranscriptome.out.bam"

#Count
bedtools coverage -counts -s -sorted -g index/sortorder_sizes.txt -b $outdir/$name".tx.bam" -a $bed | gzip > $outdir/$name".counts.gz"
