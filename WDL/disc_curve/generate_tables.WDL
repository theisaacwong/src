version 1.0

workflow discovery_curve {

    input {
        String docker = "iqwong/disccurve:0.1"

        Array[File] input_h1_files
        Array[File] input_h2_files
        File sv_vcf
        File sample_file
        File sample_order


        String n_cpu = 4
        String mem_gb = 8
    }


    call create_tables {
        input:
            bam_file = input_bams[scatter_index],
            sample_id = input_samples[scatter_index],
            callable_bed = callable_bed,
            docker = docker,
            samtools_mem_gb = samtools_mem_gb,
            samtools_n_cpu = samtools_n_cpu
    }
    
    output {
        Array[File] output_depth_files = create_callable.output_depth_file
        Array[File] output_h1_files = create_callable.output_h1_file
        Array[File] output_h2_files = create_callable.output_h2_file
    }
}

task create_tables {
    input {
        Array[File] input_h1_files
        Array[File] input_h2_files
        File sv_vcf
        File sample_file
        File sample_order
        File reference_index
        File python_script
        String name

        String docker
        String mem_gb
        String n_cpu
        Int disk_size
    }

    String out_bed = name + ".out.bed"
    String out_tsv = name + ".out.tsv"
    String out_hap = name + ".out.hap"
    String out_gts = name + ".out.gts"
    String out_counts = name + ".out.counts"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h

        mkdir -p /callable_files/

        find / | grep "_callable_regions_h._500.bed.gz" | xargs -I {} mv {} /callable_files/

        python ~{python_script} \
            --vcf ~{sv_vcf} \
            --samplefile ~{sample_file} \
            --bed ~{out_bed} \
            --tsv ~{out_tsv} \
            --haps ~{out_hap} \
            --dir /callable_files \
            --gt ~{out_gts} \
            --fai ~{reference_index} \
            --counts ~{out_counts} \
            --order ~{sample_order}
    >>>

    runtime {
        memory: mem_gb + " GB"
        cpu: n_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File output_bed = out_bed
        File output_tsv = out_tsv
        File output_hap = out_hap
        File output_gts = out_gts
        File output_counts = out_counts
    }
}

