#!/bin/bash
#avant faut créer un dossier HIC
srr=$1
dire=$2 #chemin absolu vers dossier HIC
#Le script-shell sera exécuté dans Projet_long/

mkdir $dire/data/ $dire/alignment/ $dire/results/
mkdir $dire/data/fastq/ $dire/data/ref/ $dire/data/obsolete/

cat *.fasta > genome.fa
mv genome.fa $dire/data/ref/
mv *.fasta $dire/data/ref/

fasterq-dump $srr -O $dire/data/fastq/

bowtie2-build $dire/data/ref/genome.fa genome
mv genome* $dire/data/ref/

DIR="$dire" SRR="$srr" snakemake -j 8
