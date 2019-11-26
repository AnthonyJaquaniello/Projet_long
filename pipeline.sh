#!/bin/bash
#avant faut créer un dossier HIC
srr=$1
ref=$2
dire=$3 #chemin absolu vers dossier HIC
#Le script-shell sera exécuté dans Projet_long/

conda activate projet_long

mkdir $3/data/ $3/alignment/ $3/results/
mkdir $3/data/fastq/ $3/data/ref/ $3/data/obsolete/
mv $2 genome.fa
mv $2 $3/data/ref

fastqer-dump $1 --split-3 -O $3/data/fastq/

bowtie2-build $3/data/ref/$2 genome
mv genome* $3/data/ref/

DIR="$3" SRR="$1" snakemake -j 8
