#!/bin/bash
outdir=out

first=$(ls $outdir/*_winCounts_w10*.gz | head -n1)
header="chr\tstart\tend\tname\tid\tstrand"
par="<(zcat $first | cut -f1,2,3,4,5,6)"
for file in `ls $outdir/*_winCounts_w10*.gz`
do 
	name=$(basename $file)
	name=${name%_winCounts_w10*.gz}
	header=$header"\t"$name
	par=$par" <(zcat $file | cut -f7)"
done
cmd="paste $par"

echo -e $header | gzip > winCounts_w10_s5.gz
eval $cmd | awk -F"\t" ' BEGIN {OFS="\t"} { for(i=7;i<=NF && c==0; ++i) { c+=($i>0)}  ; if(c>0) { print $0 }; c=0 }' | gzip >> winCounts_w10_s5.gz
