import os

INDEX = ['1','2','3','4','rev.1','rev.2']
SRR = os.environ.get("SRR").split()
DIR = os.environ.get("DIR").split()

rule all:
    input:
        expand('{dir}/results/{srr}.bed2', srr=SRR, dir=DIR)

rule alignment:
    input:
        left_read= '{DIR}/data/fastq/{sample}_1.fastq',
        right_read='{DIR}/data/fastq/{sample}_2.fastq',
    output:
        left_align='{DIR}/alignment/{sample}_1.sam',
        right_align='{DIR}/alignment/{sample}_2.sam'
    shell:
        """
        cat {DIR}/data/ref/NC*.fasta > {DIR}/data/ref/genome.fa
        bowtie2-build {DIR}/data/ref/genome.fa genome
        mv genome* {DIR}/data/ref/
        bowtie2  -p18 --local --very-sensitive-local --no-hd --no-sq -x {DIR}/data/ref/genome -q {input.left_read} -S {output.left_align}
        bowtie2 -p 18 --local --very-sensitive-local --no-hd --no-sq -x {DIR}/data/ref/genome -q {input.right_read} -S {output.right_align}
        mv {DIR}/data/ref/genome.fa {DIR}/data/obsolete/
        """

rule sam_extraction:
    input:
        left_old ='{DIR}/alignment/{sample}_1.sam',
        right_old='{DIR}/alignment/{sample}_2.sam'
    output:
        left_new='{DIR}/alignment/{sample}_1.sam.0',
        right_new='{DIR}/alignment/{sample}_2.sam.0'
    shell:
        """
        awk '{{print $1,$3,$4,$2,$5;}}' {input.left_old} > {output.left_new}
        awk '{{print $1,$3,$4,$2,$5;}}' {input.right_old} > {output.right_new}
        rm {input.left_old} {input.right_old}
        """

rule sam_sorting:
    input:
        left_unsorted='{DIR}/alignment/{sample}_1.sam.0',
        right_unsorted='{DIR}/alignment/{sample}_2.sam.0'
    output:
        left_sorted='{DIR}/alignment/{sample}_1.sam.0.sorted',
        right_sorted='{DIR}/alignment/{sample}_2.sam.0.sorted'
    shell:
        """
        sort -V -k1 {input.left_unsorted} > {output.left_sorted}
        sort -V -k1 {input.right_unsorted} > {output.right_sorted}
        rm {input.left_unsorted} {input.right_unsorted}
        """

rule pairing:
    input:
        left='{DIR}/alignment/{sample}_1.sam.0.sorted',
        right='{DIR}/alignment/{sample}_2.sam.0.sorted'
    output:
        '{DIR}/alignment/{sample}.merged.dat'
    shell:
        """
         paste {input.left} {input.right} > {output}
         rm {input.left} {input.right}
         """

rule quality_filtering:
    input:
        '{DIR}/alignment/{sample}.merged.dat'
    output:
        '{DIR}/alignment/{sample}.merged_qualfilt.dat'
    shell:
        """
        awk '{{if($1==$6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}}'  {input}  > {output}
        rm {input}
        """

rule fragment_restriction:
    input:
        genome='{DIR}/data/ref/',
        align='{DIR}/alignment/{sample}.merged_qualfilt.dat'
    params:
        enz='DpnII'
    output:
        '{DIR}/alignment/{sample}.merged_qualfilt.dat.indices'
    script:
        "examples_codes/fragment_attribution.py"
        
rule event_filtering:
    input:
        file='{DIR}/alignment/{sample}.merged_qualfilt.dat.indices'
    params:
        uncut_threshold='4',
        loop_threshold='2',
        srr='{sample}'
    output:
        '{DIR}/alignment/{sample}.merged_qualfilt.dat.indices.filtered'
    script:
        "examples_codes/library_events_ARG.py"
        
rule creating:
    input:
        '{DIR}/alignment/{sample}.merged_qualfilt.dat.indices.filtered'
    output:
        "{DIR}/alignment/{sample}.merged_qualfilt.dat.indices.filtered.bed2"
    script:
        "examples_codes/convert_pairs_bed2d.py"

rule post_processing:
    input:
        bad_name='{DIR}/alignment/{sample}.merged_qualfilt.dat.indices.filtered.bed2'
    output:
        good_name='{DIR}/results/{sample}.bed2'
    shell:
        """
        mv {input.bad_name} {output.good_name}
        mv *.npz {DIR}/results
        mv *.png {DIR}/results
        rm {DIR}/alignment/*.dat
        rm {DIR}/alignment/*.indices
        """

