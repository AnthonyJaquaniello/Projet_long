srr=$1
ref=$2 #path to genome.fa
dire=$3

mkdir $dire/data/ $dire/results/ $dire/alignment/
mkdir $dire/data/fastq/ $dire/data/ref/

fasterq-dump $srr -O $dire/data/fastq/
bowtie2-build $ref genome
mv *.bt2 $dire/data/ref/
bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/data/ref/genome -q $dire/data/fastq/$srr.fastq -S $dire/alignment/$srr.sam
