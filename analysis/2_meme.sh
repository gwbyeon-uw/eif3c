#!/bin/bash

#mm10 2bit available in http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit
twoBitToFa -bed=sharp.bed mm10.2bit sharp.fa
twoBitToFa -bed=bg.bed mm10.2bit bg.fa

meme sharp.fa -neg bg.fa -objfun de -oc sharp_meme -rna -minw 6 -maxw 20 -nmotifs 20 -evt 0.001 -test mrs -brief 1000000 -p 2
