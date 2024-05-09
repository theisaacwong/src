version 1.0

workflow samtoolsCollate {
    input {
        Array[File] bam_files
        Array[File] bam_index_files
        String docker

        String samtools_n_cpu = 4
        String samtools_mem_gb = 8
    }
   
    scatter (scatter_index in range(length(bam_files))) {
        call samtools_collate {
            input: 
                bam_file = bam_files[scatter_index],
                bam_index = bam_index_files[scatter_index],
                docker = docker,
                samtools_n_cpu = samtools_n_cpu,
                samtools_mem_gb = samtools_mem_gb
        }
    }

    output {
        Array[File] collated_fq_R1 = samtools_collate.reads_1,
        Array[File] collated_fq_R2 = samtools_collate.reads_2
    }
}


task samtools_collate {
    input {
        File bam_file
        File bam_index
        String docker
        String samtools_n_cpu
        String samtools_mem_gb
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 100)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_fq1 = base_filename + "_1.fastq.gz"
    String output_fq2 = base_filename + "_2.fastq.gz"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h
        samtools collate -@ ${nthreads} -u -O ~{bam_file} | samtools fastq -@ ${nthreads} -1 ~{output_fq1} -2 ~{output_fq2} -0 /dev/null -s /dev/null -n
        # samtools collate -@ ${nthreads} -o TEMP_BAM.bam ~{bam_file}
        # samtools fastq -@ ${nthreads} -1 ~{output_fq1} -2 ~{output_fq2} -0 /dev/null -s /dev/null -n TEMP_BAM.bam
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

