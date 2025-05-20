import pandas as pd
import os, sys

configfile: "config.yaml"

MANIFEST = config.get("MANIFEST", "manifest.tab")
DOCKER = config.get("DOCKER", "docker://iqwong/altina:1.2")
NTHREADS = config.get("NTHREADS", 1)
REF = config.get("REF", "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta")

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="SAMPLE")


def get_input(wildcards):
    return manifest_df.at[wildcards.sample, "INPUT"]

def get_h1(wildcards):
    return manifest_df.at[wildcards.sample, "ASM_QC_HAP1"]

def get_h2(wildcards):
    return manifest_df.at[wildcards.sample, "ASM_QC_HAP2"]

def get_hifi(wildcards):
    return manifest_df.at[wildcards.sample, "ALIGNMENT_T2T"]

def get_ref(wildcards):
    return manifest_df.at[wildcards.sample, "REF"]


wildcard_constraints:
    sample="|".join(manifest_df.index),

localrules:
    all,

rule all:
    input:
        expand(
            "hapdiff/{sample}/hapdiff_phased.vcf.gz",
            sample=manifest_df.index
        ),
        expand(
            "sawfish/{sample}/sawfish_call/genotyped.sv.vcf.gz",
            sample=manifest_df.index
        ),

rule hapdiff:
    input:
        h1 = get_h1,
        h2 = get_h2,
        ref = REF
    output:
        vcf="hapdiff/{sample}/hapdiff_phased.vcf.gz",
    resources:
        mem=20,
        hrs=24,
        disk_free=1,
    threads: 6
    shell:        """
        hapdiff.py --reference {input.ref} --pat {input.h1} --mat {input.h2} --out-dir $( dirname {output.vcf} ) -t {threads} --sample {wildcards.sample} --sv-size 50
        """


rule sawfish_discover:
    input:
        bam = get_hifi,
        ref = get_ref
    output:
        bcf = 'sawfish/{sample}/sawfish_disc/candidate.sv.bcf'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "sawfish/0.12.4",
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:        """
        sawfish discover --threads {threads} --ref {input.ref} --bam {input.bam} --output-dir $( dirname {output.bcf} ) --clobber
        """

rule sawfish_call:
    input:
        bcf = 'sawfish/{sample}/sawfish_disc/candidate.sv.bcf'
    output:
        vcf = 'sawfish/{sample}/sawfish_call/genotyped.sv.vcf.gz'
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:        """
        sawfish joint-call --threads {threads} --sample $( dirname {input.bcf} ) --clobber --output-dir $( dirname {output.vcf} )
        """