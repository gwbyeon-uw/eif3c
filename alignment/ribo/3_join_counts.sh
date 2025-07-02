#!/bin/bash

outdir=out

first=$(ls $outdir/*counts.gz | head -n1)
header="chr\tstart\tend\tname\tid\tstrand"
par="<(zcat $first | cut -f1,2,3,4,5,6)"
for file in `ls $outdir/*counts.gz`
do 
	name=$(basename $file)
	name=${name%.counts.gz}
	header=$header"\t"$name
	par=$par" <(zcat $file | cut -f7)"
done
cmd="paste $par"
echo $cmd

echo -e $header | gzip > counts.gz
eval $cmd | gzip >> counts.gz
