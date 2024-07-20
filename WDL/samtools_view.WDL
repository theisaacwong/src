version 1.0

workflow samtoolsView {
    input {
        Array[File] bam_files
        Array[File] bam_index_files
        String docker
        File reference_fasta
        File reference_fai
        File reference_dict

        String samtools_n_cpu = 4
        String samtools_mem_gb = 8
    }
   
    scatter (scatter_index in range(length(bam_files))) {
        call samtools_view {
            input: 
                bam_file = bam_files[scatter_index],
                bam_index = bam_index_files[scatter_index],
                docker = docker,
                samtools_n_cpu = samtools_n_cpu,
                samtools_mem_gb = samtools_mem_gb,
                reference_fasta = reference_fasta,
                reference_fai = reference_fai,
                reference_dict = reference_dict
        }
    }

    output {
        Array[File] chr1_crams = samtools_view.cram_out
        Array[File] chr1_crais = samtools_view.crai_out
    }
}


task samtools_view {
    input {
        File bam_file
        File bam_index
        String docker
        String samtools_n_cpu
        String samtools_mem_gb

        File reference_fasta
        File reference_fai
        File reference_dict
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 100)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_chr1 = base_filename + ".chr1.cram"
    String output_chr1_crai = base_filename + ".chr1.cram.crai"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        free -h

        samtools view --cram -h -o ~{output_chr1} -@ ${nthreads} -T ~{reference_fasta} ~{bam_file} chr1
        samtools index ~{output_chr1}
    >>>

    runtime {
        memory: samtools_mem_gb + " GB"
        cpu: samtools_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File cram_out = output_chr1
        File crai_out = output_chr1_crai
    }
}

