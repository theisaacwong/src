import pandas as pd


configfile: "config.yaml"


MANIFEST = config.get("MANIFEST", "manifest.tab")
REF = config["REF"]
TRF_BED = config["TRF_BED"]
PBSV_VERSION = config["PBSV_VERSION"]

manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col="SAMPLE")

pbsv_dict = {"sv": "DEL,INS,INV", "dup": "DUP", "bnd": "BND"}


def find_aln(wildcards):
    return manifest_df.at[wildcards.sample, "ALN"]


def find_svtype(wildcards):
    return pbsv_dict[wildcards.svtype]


wildcard_constraints:
    sample="|".join(manifest_df.index),


localrules:
    all,


rule all:
    input:
        expand(
            "results/{sample}/pbsv_{sample}_{svtype}.vcf.gz",
            sample=manifest_df.index,
            svtype=pbsv_dict,
        ),


rule pbsv_svsig:
    input:
        bam=find_aln,
        trf=TRF_BED,
    output:
        svsig="tmp/{sample}/{sample}.svsig.gz",
    singularity:
        f"docker://eichlerlab/pbsv:{PBSV_VERSION}"
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=24,
    shell:
        """
        pbsv discover --tandem-repeats {input.trf} {input.bam} {output.svsig}
        """


rule pbsv_call:
    input:
        ref=REF,
        svsig=rules.pbsv_svsig.output.svsig,
    output:
        vcf=temp("results/{sample}/pbsv_{sample}_{svtype}.vcf"),
    singularity:
        f"docker://eichlerlab/pbsv:{PBSV_VERSION}"
    threads: 4
    params:
        vcf_out=find_svtype,
    resources:
        mem=lambda wildcards, attempt: attempt * 32,
        hrs=24,
    shell:
        """
        pbsv call -j {threads} --ccs --types {params.vcf_out} {input.ref} {input.svsig} {output.vcf} 
        """


rule gzip_index:
    input:
        vcf=rules.pbsv_call.output.vcf,
    output:
        vcf="results/{sample}/pbsv_{sample}_{svtype}.vcf.gz",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 32,
        hrs=24,
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf}
        sleep 120s
        tabix -p vcf -f {output.vcf}
        """
