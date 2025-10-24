import pandas as pd
import os, sys

configfile: "config.yaml"

MANIFEST = config.get("MANIFEST", "manifest.tab")
DOCKER = config.get("DOCKER", "docker://iqwong/altina:1.2")
NTHREADS = config.get("NTHREADS", 1)

REF_DICT = config['REF']

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="SAMPLE")

def get_h1(wildcards):
    return manifest_df.at[wildcards.sample, "ASM_QC_HAP1"]

def get_h2(wildcards):
    return manifest_df.at[wildcards.sample, "ASM_QC_HAP2"]

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
            "{ref}/hapdiff/{sample}/hapdiff_phased.vcf.gz",
            sample=manifest_df.index,
            ref=REF_DICT,
        ),

rule hapdiff:
    input:
        h1 = get_h1,
        h2 = get_h2,
        ref = get_ref,
    output:
        vcf='{ref}/hapdiff/{sample}/hapdiff_phased.vcf.gz',
    resources:
        mem=20,
        hrs=24,
    threads: 6
    shell:        """
        hapdiff.py --reference {input.ref} --pat {input.h1} --mat {input.h2} --out-dir $( dirname {output.vcf} ) -t {threads} --sample {wildcards.sample} --sv-size 50
        """

