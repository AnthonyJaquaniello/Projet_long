#!/bin/bash

srr=$1 #path to srr list file .txt
input=$2 #idem
dire=$3 #path to location where dir root will be created

list_real=$(cat $srr)
list_input=$(cat $input)

echo 'Constructing the arborescence...'
mkdir $dire/CHIPSEQ/
mkdir $dire/CHIPSEQ/data/ $dire/CHIPSEQ/results/ $dire/CHIPSEQ/alignment/
mkdir $dire/CHIPSEQ/data/fastq/ $dire/CHIPSEQ/data/ref/ $dire/CHIPSEQ/results/input/ $dire/CHIPSEQ/results/IP/
echo '--> done !'

echo 'Indexing of the refrence...'
cp $dire/HIC/data/obsolete/genome.fa  $dire/CHIPSEQ/data/ref/genome.fa
bowtie2-build $dire/CHIPSEQ/data/ref/genome.fa genome
mv *.bt2 $dire/CHIPSEQ/data/ref/
echo '--> done !'

echo 'Traitment of all real chip-seq experiment...'
for i in $list_real; do
	fasterq-dump "$i" -t $dire/CHIPSEQ -O $dire/CHIPSEQ/data/fastq/
	bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$i.fastq -S $dire/CHIPSEQ/alignment/$i.sam
	samtools view -Sb $dire/CHIPSEQ/alignment/$i.sam > $dire/CHIPSEQ/alignment/$i.bam
	samtools 'sort' $dire/CHIPSEQ/alignment/$i.bam > $dire/CHIPSEQ/alignment/$i.sorted.bam
	samtools index $dire/CHIPSEQ/alignment/$i.sorted.bam $dire/CHIPSEQ/alignment/$i.sorted.bam.bai
	rm $dire/CHIPSEQ/alignment/$i.sam $dire/CHIPSEQ/alignment/$i.bam
	echo 'Peak extraction...'
	mkdir $dire/CHIPSEQ/results/IP/$i/
	python3 examples_codes/peaks_extract.py $dire/CHIPSEQ/alignment/$i.sorted.bam $dire/CHIPSEQ/results/IP/$i/
	echo $i 'terminated'
done
echo '--> done !'
echo 'Traitment of all input controls...'
for i in $list_input; do
	echo 'Downloading of data...'
	fasterq-dump "$i" -t $dire/CHIPSEQ -O $dire/CHIPSEQ/data/fastq/
	mv $dire/CHIPSEQ/data/fastq/$i.fastq $dire/CHIPSEQ/data/fastq/$i.input.fastq
	echo 'Alignment...'
	bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$i.input.fastq -S $dire/CHIPSEQ/alignment/$i.input.sam
	echo 'Conversion in bam format...'
	samtools view -Sb $dire/CHIPSEQ/alignment/$i.input.sam > $dire/CHIPSEQ/alignment/$i.input.bam
	echo 'Sorting...'
	samtools 'sort' $dire/CHIPSEQ/alignment/$i.input.bam > $dire/CHIPSEQ/alignment/$i.input.sorted.bam
	echo 'Indexing...'
	samtools index $dire/CHIPSEQ/alignment/$i.input.sorted.bam $dire/CHIPSEQ/alignment/$i.input.sorted.bam.bai
	rm $dire/CHIPSEQ/alignment/$i.input.sam $dire/CHIPSEQ/alignment/$i.input.bam
	echo 'Peak extraction...'
	mkdir $dire/CHIPSEQ/results/input/$i/
	python3 examples_codes/peaks_extract.py $dire/CHIPSEQ/alignment/$i.input.sorted.bam $dire/CHIPSEQ/results/input/$i/
	echo $i 'terminated'
done
echo '--> done'
