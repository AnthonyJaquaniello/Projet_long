srr=$1
ref=$2 #path to chr
dire=$3

mkdir $dire/data/ $dire/results/ $dire/alignment/
mkdir $dire/data/fastq/ $dire/data/ref/

#fasterq-dump $srr -O $dire/data/fastq/

cat $ref/*.fasta > $dire/data/ref/genome.fa
chr=$(ls $ref/*.fasta) # contain all paths to chr.fasta
for i in $chr; do
    echo $i
    bowtie2-build $i "$i"
    mv *.bt2 $dire/data/ref/
done
bowtie2-build $dire/data/ref/genome.fa genome
mv *.bt2 $dire/data/ref/

bowtie2 -p8 '--local' '--very-sensitive-local' '-x' $dire/data/ref/genome -q $dire/data/fastq/$srr.fastq -S $dire/alignment/$srr.sam
