#!/bin/bash
#SBATCH -n 16
#SBATCH --mem-per-cpu=4G

mmrDNAFasta=index/mm_rDNA/mouse_rDNA.fa
mm10Fasta=index/mm10/fasta
index=index
gtf=index/mm10/gencode.vM18.annotation.gtf

mkdir -p $index/star

#Make STAR index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $index/star --genomeFastaFiles $mmrDNAFasta $mm10Fasta/*.fa --sjdbGTFfile $gtf --sjdbOverhang 50

#Create BED file for coverage
bedtools makewindows -w 20 -s 10 -b mm10_all.bed | awk -F"\t" ' BEGIN { OFS="\t" } { print $1, $2, $3, $1"_"$2"_"$3"_+",(NR-1)*2,"+"; print $1, $2, $3, $1"_"$2"_"$3"_-",(NR-1)*2+1,"-" } ' > mm10_all_w20_s10.bed
