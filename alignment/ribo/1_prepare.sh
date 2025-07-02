#!/bin/bash
#SBATCH -n 16
#SBATCH --mem-per-cpu=4G

mmrDNAFasta=index/mm_rDNA/mouse_rDNA.fa
mm10Fasta=index/mm10/fasta
index=index
gtf=index/mm10/gencode.vM18.annotation.gtf

mkdir $index/star

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $index/star --genomeFastaFiles $mmrDNAFasta $mm10Fasta/*.fa --sjdbGTFfile $gtf --sjdbOverhang 50

