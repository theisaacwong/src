version 1.0

workflow discovery_curve {

    input {
        String docker

        File callable_bed
        Array[File] input_bams
        Array[String] input_samples


        File truvari_vcf
        File sample_file
    }


    scatter (scatter_index in range(length(input_bams))) {
        call create_callable {
            input:
                bam_file = input_bams[scatter_index],
                sample_id = input_samples[scatter_index],
                callable_bed = callable_bed
                docker = docker
        }
    }

    call create_data_table {


    }

    call create_plot {


    }
}

task create_callable {
    input {
        File bam_file
        String sample_id
        String docker
        File callable_bed
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 50)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_depth = base_filename + ".depth.txt" 

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h
        
        samtools depth -@ ${nthreads} -b ~{callable_bed} ~{bam_file} > ~{output_depth}

        python <<CODE
            


        CODE
        >>>
    >>>

    runtime {
        memory: samtools_mem_gb + " GB"
        cpu: samtools_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File reads_1 = output_fq1
        File reads_2 = output_fq2
    }
}