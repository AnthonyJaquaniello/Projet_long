#!/bin/bash

srr=$1 #path to srr list file .txt
input=$2 #idem
dire=$3 #path to location where dir root will be created

list_real=$(cat $srr)
list_input=$(cat $input)

real_tab=()
input_tab=()
a=0
b=0


for i in $list_real
do
	real_tab[$a]=$i
	let 'a+=1'
done

for j in $list_input
do
	input_tab[$b]=$j
	let 'b+=1'
done


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
	echo 'Downloading of the data...'
    	fasterq-dump "$i" -t $dire/CHIPSEQ -O $dire/CHIPSEQ/data/fastq/
	if [ -e $dire/CHIPSEQ/data/fastq/$i'_1.fastq' ] && [ -e $dire/CHIPSEQ/data/fastq/$i'_2.fastq' ] #if paired-end
	then
	    echo '[PAIRED-END] Alignment...'
	    bowtie2 -p8 '-x' $dire/CHIPSEQ/data/ref/genome -1 $dire/CHIPSEQ/data/fastq/$i'_1.fastq' -2 $dire/CHIPSEQ/data/fastq/$i'_2.fastq' -S $dire/CHIPSEQ/alignment/$i'.sam'
	    echo '[PAIRED-END] Sorting...'
	    sort -V -k1 $dire/CHIPSEQ/alignment/$i'.sam' > $dire/CHIPSEQ/alignment/$i'.sorted.sam'
	    echo '[PAIRED-END] Conversion in bam format...'
	    samtools view -Sb $dire/CHIPSEQ/alignment/$i'.sorted.sam' > $dire/CHIPSEQ/alignment/$i'.sorted.bam'
	    rm $dire/CHIPSEQ/alignment/$i'.sam' 
	    echo '[PAIRED-END] Indexing...'
	    samtools index $dire/CHIPSEQ/alignment/$i'.sorted.bam' $dire/CHIPSEQ/alignment/$i'.sorted.bam.bai'
	    echo $i 'terminated'
	else
	    echo '[SINGLE-END] Alignment...'
	    bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$i.fastq -S $dire/CHIPSEQ/alignment/$i.sam
	    echo '[SINGLE-END] Conversion in bam format...'
	    samtools view -Sb $dire/CHIPSEQ/alignment/$i.sam > $dire/CHIPSEQ/alignment/$i.bam
	    rm $dire/CHIPSEQ/alignment/$i.sam    
	    echo '[SINGLE-END] Sorting...'
	    samtools 'sort' $dire/CHIPSEQ/alignment/$i.bam > $dire/CHIPSEQ/alignment/$i.sorted.bam
	    rm $dire/CHIPSEQ/alignment/$i.bam
	    echo '[SINGLE-END] Indexing...'
	    samtools index $dire/CHIPSEQ/alignment/$i.sorted.bam $dire/CHIPSEQ/alignment/$i.sorted.bam.bai
	    echo $i 'terminated'
	fi
done


echo 'Traitment of all input controls...'
for i in $list_input; do
	echo $i
	echo 'Downloading of data...'
	fasterq-dump "$i" -t $dire/CHIPSEQ -O $dire/CHIPSEQ/data/fastq/
	if [ -e $dire/CHIPSEQ/data/fastq/$i'_1.fastq' ] && [ -e $dire/CHIPSEQ/data/fastq/$i'_2.fastq' ] 
	then
	    mv $dire/CHIPSEQ/data/fastq/$i'_1.fastq' $dire/CHIPSEQ/data/fastq/$i'_1.input.fastq'
	    mv $dire/CHIPSEQ/data/fastq/$i'_2.fastq' $dire/CHIPSEQ/data/fastq/$i'_2.input.fastq'
	    echo '[PAIRED-END] Alignment...'
	    bowtie2 -p8 '-x' $dire/CHIPSEQ/data/ref/genome -1 $dire/CHIPSEQ/data/fastq/$i'_1.input.fastq' -2 $dire/CHIPSEQ/data/fastq/$i'_2.input.fastq' -S $dire/CHIPSEQ/alignment/$i'.input.sam'
	    echo '[PAIRED-END] Sorting...'
	    sort -V -k1 $dire/CHIPSEQ/alignment/$i'.input.sam' > $dire/CHIPSEQ/alignment/$i'.input.sorted.sam'
	    echo '[PAIRED-END] Conversion in bam format...'
	    samtools view -Sb $dire/CHIPSEQ/alignment/$i'.input.sorted.sam' > $dire/CHIPSEQ/alignment/$i'.input.sorted.bam'
	    rm $dire/CHIPSEQ/alignment/$i'.input.sam' 
	    echo '[PAIRED-END] Indexing...'
	    samtools index $dire/CHIPSEQ/alignment/$i'.input.sorted.bam' $dire/CHIPSEQ/alignment/$i'.input.sorted.bam.bai'
	    echo $i 'terminated'
	else
	    mv $dire/CHIPSEQ/data/fastq/$i'.fastq' $dire/CHIPSEQ/data/fastq/$i'.input.fastq' 
	    echo '[SINGLE-END] Alignment...'
	    bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$i.input.fastq -S $dire/CHIPSEQ/alignment/$i.input.sam
	    echo '[SINGLE-END] Conversion in bam format...'
	    samtools view -Sb $dire/CHIPSEQ/alignment/$i.input.sam > $dire/CHIPSEQ/alignment/$i.input.bam
	    rm $dire/CHIPSEQ/alignment/$i.input.sam
	    echo '[SINGLE-END] Sorting...'
	    samtools 'sort' $dire/CHIPSEQ/alignment/$i.input.bam > $dire/CHIPSEQ/alignment/$i.input.sorted.bam
	    rm $dire/CHIPSEQ/alignment/$i.input.bam
	    echo '[SINGLE-END] Indexing...'
	    samtools index $dire/CHIPSEQ/alignment/$i.input.sorted.bam $dire/CHIPSEQ/alignment/$i.input.sorted.bam.bai
	    echo $i 'terminated'
	fi
done

echo 'Extraction and cleaning of chip-seq peaks...'
d=$a
let 'd-=1'
for c in `seq 0 $d`
do
    echo ${real_tab[$c]}
    echo ${input_tab[$c]}
    macs2 callpeak -t $dire/CHIPSEQ/alignment/${real_tab[$c]}.sorted.bam -c $dire/CHIPSEQ/alignment/${input_tab[$c]}.input.sorted.bam -n chipseq_peaks --outdir $dire/CHIPSEQ/results/IP/${real_tab[$c]}/ --format BAM --nomodel
done
echo 'All chipseq experiment are terminated !'
