srr=$1
input=$2
dire=$3
ref=$4 #path to genome.fa

mkdir $dire/CHIPSEQ/
mkdir $dire/CHIPSEQ/data/ $dire/CHIPSEQ/results/ $dire/CHIPSEQ/alignment/
mkdir $dire/CHIPSEQ/data/fastq/ $dire/CHIPSEQ/data/ref/ $dire/CHIPSEQ/results/input/ $dire/CHIPSEQ/results/IP/

fasterq-dump $srr -O $dire/CHIPSEQ/data/fastq/
fasterq-dump $input -O $dire/CHIPSEQ/data/fastq/
mv $dire/CHIPSEQ/data/fastq/$input.fastq $dire/CHIPSEQ/data/fastq/$input.input.fastq

cp $ref $dire/CHIPSEQ/data/ref/genome.fa
bowtie2-build $dire/CHIPSEQ/data/ref/genome.fa genome
mv *.bt2 $dire/CHIPSEQ/data/ref/

bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$srr.fastq -S $dire/CHIPSEQ/alignment/$srr.sam
bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/CHIPSEQ/data/ref/genome -q $dire/CHIPSEQ/data/fastq/$input.input.fastq $dire/CHIPSEQ/alignment/$input.input.sam

samtools view -Sb $dire/CHIPSEQ/alignment/$srr.sam > $dire/CHIPSEQ/alignment/$srr.bam
samtools view -Sb $dire/CHIPSEQ/alignment/$input.input.sam > $dire/CHIPSEQ/alignment/$input.input.bam
samtools 'sort' $dire/CHIPSEQ/alignment/$srr.bam > $dire/CHIPSEQ/alignment/$srr.sorted.bam
samtools 'sort' $dire/CHIPSEQ/alignment/$input.input.bam > $dire/CHIPSEQ/alignment/$input.input.sorted.bam
samtools index $dire/CHIPSEQ/alignment/$srr.sorted.bam $dire/CHIPSEQ/alignment/$srr.sorted.bam.bai
samtools index $dire/CHIPSEQ/alignment/$input.input.sorted.bam $dire/CHIPSEQ/alignment/$input.input.sorted.bam.bai

python3 examples_codes/peaks_extract.py $dire/CHIPSEQ/alignment/$srr.sorted.bam $dire/CHIPSEQ/results/IP/
python3 examples_codes/peaks_extract.py $dire/CHIPSEQ/alignment/$input.input.sorted.bam $dire/CHIPSEQ/results/input/
