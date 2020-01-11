#!/bin/bash

hic_srr=$1 #path to file containing SRR of HIC experiment
chip_srr=$2 #path to file containing SRR of real Chipseq experiment
input_srr=$3 #path to file containing SRR of negative control Chipseq experiment
dire=$4 #where the roots directories should be created, don't write / at the end
chr_fasta=$5 #path to all chr.fasta files, don't write / at the end

./pipeline_hic.sh $hic_srr $dire $chr_fasta

./pipeline_chipseq.sh $chip_srr $input_srr $dire
hic=()
chip=()
a=0
b=0

for i in $(cat $hic_srr)
do
	hic[$a]=$i
	let 'a+=1'
done

for j in $(cat $chip_srr)
do
	chip[$b]=$j
	let 'b+=1'
done

d=$a
let 'd-=1'

for c in `seq 0 $d`
do
    echo ${hic[$c]}
    echo ${chip[0]}
	python3 examples_codes/proportion_chip_peak_other_impl.py $dire/HIC/results/${hic[$c]} $dire/CHIPSEQ/results/IP/${chip[0]}/chipseq_peaks.txt ${hic[$c]}_${chip[0]}
done
echo 'All the experiment are terminated ! Enrichiment plot is in the HIC/results/ directory.'
