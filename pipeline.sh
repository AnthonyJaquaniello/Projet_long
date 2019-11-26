#!/bin/bash
#avant faut créer un dossier HIC
srr=$1
ref=$2
dire=$3 #chemin absolu vers dossier HIC
#Le script-shell sera exécuté dans Projet_long/

mkdir $3/data/ $3/alignment/ $3/results/
mkdir $3/data/fastq/ $3/data/ref/ $3/data/obsolete/
mv $2 genome.fa
mv genome.fa $3/data/ref

fasterq-dump $1 -O $3/data/fastq/

bowtie2-build $3/data/ref/genome.fa genome
mv genome* $3/data/ref/

DIR="$3" SRR="$1" snakemake -j 8
