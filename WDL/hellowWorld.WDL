version 1.0

workflow testWorkflow {
    input {
        Array[File] bams
        Array[File] bais
        Array[String] sample_ids
        String docker

    }
    scatter (scatter_index in range(length(sample_ids))) {
        call getHeader as headers {
            input: 
                bam = bams[scatter_index],
                bai = bais[scatter_index],
                sample_id = sample_ids[scatter_index],
                docker = docker
            }
    }
    output {
        Array[String] sample_all = headers.sample_id
        Array[file] header = headers.header_file 
    }
}

task getHeader {
    input {
        File bam
        File bai
        String sample_id
        String docker
        
    }

    Int disk_size = ceil(size(bam, "GiB") + 1)
    String base_filename = basename(bam, ".bam")
    String output_filename = base_filename + ".header"

    command <<<
        set -euo pipefail
        samtools view -H ~{bam} > ~{output_filename}
    >>>

    runtime {
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: ~{docker}
    }

    output {
        sample_id = sample_id,
        header_file = output_filename
    }
}