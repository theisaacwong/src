import pandas as pd
import os, sys

configfile: "config.yaml"

MANIFEST = config.get("MANIFEST", "manifest.tab")
DOCKER = config.get("DOCKER", "docker://iqwong/altina:1.2")
NTHREADS = config.get("NTHREADS", 1)

REF_DICT = config['REF']

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="SAMPLE")


def get_aln(wildcards):
    return manifest_df.at[wildcards.sample, wildcards.ref]

def get_ref(wildcards):
    return REF_DICT[wildcards.ref]


wildcard_constraints:
    sample="|".join(manifest_df.index),
    ref='|'.join(REF_DICT),


localrules:
    all,

rule all:
    input:
        expand(
            "{ref}/sawfish/{sample}/sawfish_call/genotyped.sv.vcf.gz",
            sample=manifest_df.index,
            ref=REF_DICT
        ),


rule sawfish_discover:
    input:
        bam = get_aln,
        ref = get_ref
    output:
        bcf = '{ref}/sawfish/{sample}/sawfish_disc/candidate.sv.bcf'
    resources:
        mem = lambda wildcards, attempt: attempt * 12 + 50,
        hrs = 24,
    threads: 4
    shell:        """
        sawfish discover --threads {threads} --ref {input.ref} --bam {input.bam} --output-dir $( dirname {output.bcf} ) --clobber
        """

rule sawfish_call:
    input:
        bcf = '{ref}/sawfish/{sample}/sawfish_disc/candidate.sv.bcf'
    output:
        vcf = '{ref}/sawfish/{sample}/sawfish_call/genotyped.sv.vcf.gz'
    resources:
        mem = lambda wildcards, attempt: attempt * 12 + 50,
        hrs = 24,
    threads: 4
    shell:        """
        sawfish joint-call --threads {threads} --sample $( dirname {input.bcf} ) --clobber --output-dir $( dirname {output.vcf} )
        """
