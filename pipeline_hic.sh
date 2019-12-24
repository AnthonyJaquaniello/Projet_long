#!/bin/bash

srr=$1 #path to file containing SRR of HIC experiment
dire=$2 #where the root directory (HIC) should be created, don't write / at the end
ref=$3 #path to all chr.fasta files, don't write / at the end

list_srr=$(cat $srr)

echo 'Constructing the arborescence...'
mkdir $dire/HIC/
mkdir $dire/HIC/data/ $dire/HIC/alignment/ $dire/HIC/results/
mkdir $dire/HIC/data/fastq/ $dire/HIC/data/ref/ $dire/HIC/data/obsolete/

cat $ref/*.fasta > $dire/HIC/data/ref/genome.fa
cp $ref/*.fasta $dire/HIC/data/ref/

echo 'Downloading fastq file...'
for i in $list_srr;do
	fasterq-dump $i -O $dire/HIC/data/fastq/
done
echo '--> Done !'

echo 'Indexing the genome...'
bowtie2-build $dire/HIC/data/ref/genome.fa genome
mv *.bt2 $dire/HIC/data/ref/
echo '--> Done !'

for i in $list_srr;do
	DIR="$dire/HIC" SRR="$i" snakemake -j 8
	chromosight detect $dire/HIC/results/$i/*.bg2 --output $dire/HIC/results/$i/
	python3 examples_codes/loops_visualisation_yeast.py $dire/HIC/results/$i $dire/HIC/results/$i/loops_out.txt
done
