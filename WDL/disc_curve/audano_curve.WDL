version 1.0

workflow discovery_curve {

    input {
        String docker = "iqwong/disccurve:0.1"

        File callable_bed
        Array[File] input_bams
        Array[String] input_samples

        String samtools_n_cpu = 4
        String samtools_mem_gb = 8


        File truvari_vcf
        File sample_file
    }


    scatter (scatter_index in range(length(input_bams))) {
        call create_callable {
            input:
                bam_file = input_bams[scatter_index],
                sample_id = input_samples[scatter_index],
                callable_bed = callable_bed,
                docker = docker,
                samtools_mem_gb = samtools_mem_gb,
                samtools_n_cpu = samtools_n_cpu
        }
    }

    call create_table {
        input:
            truvari_vcf = truvari_vcf,
            sample_file = sample_file,
            out_bed = out_bed,
            out_tsv = out_tsv,
            out_haps = out_haps,
    }



    output {
        Array[File] output_depth_files = create_callable.output_depth_file
        Array[File] output_h1_files = create_callable.output_h1_file
        Array[File] output_h2_files = create_callable.output_h2_file
    }
}

task create_callable {
    input {
        File bam_file
        String sample_id
        String docker
        File callable_bed
        String samtools_mem_gb
        String samtools_n_cpu
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 50)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_depth = base_filename + ".depth.txt" 
    String output_h1 = base_filename + "_callable_regions_h1_500.bed.gz"
    String output_h2 = base_filename + "_callable_regions_h2_500.bed.gz"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h
        
        samtools depth -@ ${nthreads} -b ~{callable_bed} ~{bam_file} > ~{output_depth}

        python <<CODE
        import gzip
        import pandas as pd
        import pysam
        from pybedtools import BedTool

        df = pd.read_csv(
        "~{output_depth}",sep="\t",header=None,names=["#CHROM", "POS", "DEPTH"]
        )
        df["END"] = df["POS"] + 1
        df_h1 = df.loc[df["DEPTH"] > 0][["#CHROM", "POS", "END"]].copy()
        df_h2 = df.loc[df["DEPTH"] >= 7][["#CHROM", "POS", "END"]].copy()
        df_h1 = (
            BedTool.from_dataframe(df_h1)
            .merge(d=500)
            .to_dataframe(names=["#CHROM", "POS", "END"])
        )
        df_h2 = (
            BedTool.from_dataframe(df_h2)
            .merge(d=500)
            .to_dataframe(names=["#CHROM", "POS", "END"])
        )
        df_h1.to_csv("~{output_h1}",sep="\t",index=False)
        df_h2.to_csv("~{output_h2}",sep="\t",index=False)
            
        CODE
    >>>

    runtime {
        memory: samtools_mem_gb + " GB"
        cpu: samtools_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File output_depth_file = output_depth
        File output_h1_file = output_h1
        File output_h2_file = output_h2
    }
}


task create_data_table {
    input {
        File python_script
        File truvari_vcf
        File sample_file

        String sample_id
        String docker
        File callable_bed
        String data_table_mem_gb
        String data_table_n_cpu
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 50)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_depth = base_filename + ".depth.txt" 
    String output_h1 = base_filename + "_callable_regions_h1_500.bed.gz"
    String output_h2 = base_filename + "_callable_regions_h2_500.bed.gz"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h
        
        python ~{python_script} 

    >>>

    runtime {
        memory: data_table_mem_gb + " GB"
        cpu: data_table_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File output_depth_file = output_depth
        File output_h1_file = output_h1
        File output_h2_file = output_h2
    }
}
