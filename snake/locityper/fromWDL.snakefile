import pandas as pd
import os, sys

    # command <<<
    #     set -euxo pipefail
    #     nthreads=$(nproc)
    #     echo "using ${nthreads} threads"
    #     gunzip -c ~{reference} > reference.fa
    #     samtools faidx reference.fa
    #     locityper add \
    #         -@ ${nthreads} \
    #         -d vcf_db \
    #         -v ~{vcf} \
    #         -r reference.fa \
    #         -j ~{counts_jf} \
    #         -L ~{bed}
    #     echo "compressing DB"
    #     tar -czf ~{output_tar} vcf_db
    #     echo "done compressing DB"
    # >>>

configfile: "config.yaml"

MANIFEST = config.get("MANIFEST", "manifest.tab")
REF = config['REF']
locityper_db = config['db']
BED_FILE = config['bed']
JF_FILE = config['counts_jf']

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="SAMPLE")

def get_cram(wildcards):
    return manifest_df.at[wildcards.sample, "cram"]

def get_crai(wildcards):
    return manifest_df.at[wildcards.sample, "crai"]

wildcard_constraints:
    sample="|".join(manifest_df.index),

localrules:
    all,


rule all:
    input:
        expand(
            "results/gts.{sample}.csv",
            sample=manifest_df.index
        ),

rule LocityperPreprocessAndGenotype:
    input:
        cram=get_cram,
        crai=get_crai,
        ref=REF,
        counts=JF_FILE,
        db=locityper_db,
        bed=BED_FILE,
    output:
        output_file="results/{sample}/success.txt",
    threads: 8,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=24,
    singularity:
        "docker://us.gcr.io/broad-dsp-lrma/lr-locityper:1.0.0",
    shell:        """
        set -eux; mkdir -p results/{wildcards.sample}
        hostname; date; echo "preprocessing and genotyping {wildcards.sample}";

        mkdir -p results/{wildcards.sample}/locityper_preproc

        locityper preproc -a {input.cram} \
            -r {input.ref} \
            -j {input.counts} \
            -@ {threads} \
            --technology illumina \
            -o results/{wildcards.sample}/locityper_preproc

        mkdir -p results/{wildcards.sample}/out_dir
        mkdir -p results/{wildcards.sample}/out_dir/loci

        process_single_locus() {{
            line="$1"
            locus_name=$(echo "$line" | cut -f4)
            echo "Processing locus: $locus_name"

            # Ensure the locus-specific directory exists
            mkdir -p "results/{wildcards.sample}/out_dir/loci/$locus_name"

            locityper genotype -a {input.cram} \
                -r {input.ref} \
                -d {input.db} \
                -p results/{wildcards.sample}/locityper_preproc \
                -S greedy:i=5k,a=1 -S anneal:i=20,a=20 \
                --subset-loci $locus_name \
                -o results/{wildcards.sample}/out_dir
        }}
        export -f process_single_locus

        cat {input.bed} | parallel --line-buffer -j {threads} process_single_locus {{}}

        date
        touch {output.output_file}
        """

rule Summarize:
    input:
        input="results/{sample}/success.txt",
    output:
        summary="results/gts.{wildcards.sample}.csv",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: attempt * 8,
        hrs=24,
    singularity:
        "docker://us.gcr.io/broad-dsp-lrma/lr-hidive:0.1.122"
    shell:        """
        cd results
        /usr/local/bin/python3 /locityper/extra/into_csv.py -i ./{wildcards.sample} -o gts.{wildcards.sample}.csv
        """

