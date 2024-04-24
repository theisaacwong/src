version 1.0

workflow runLocityper {
    input {
        Array[File] bam_files
        Array[File] bam_index_files
        Array[String] sample_ids
        String docker
        File reference
        File reference_index
        File counts_jf
        File db_path
        String locus_name
        String locus_coordinates
        File alleles_fa
    }
   
    call generate_counts {
    input:
        reference = reference,
        docker = docker
    }

    call generate_db {
    input:
        reference = reference,
        reference_index = reference_index,
        counts_jf = generate_counts.counts_file,
        locus_name = locus_name,
        locus_coordinates = locus_coordinates,
        alleles_fa = alleles_fa,
        docker = docker
    }
    
    scatter (scatter_index in range(length(bam_files))) {
        call samtools_collate {
            input: 
                bam_file = bam_files[scatter_index],
                bam_index = bam_index_files[scatter_index],
                docker = docker
        }
        call locityper_preprocess {
            input:
                input_fq1 = samtools_collate.reads_1,
                input_fq2 = samtools_collate.reads_2,
                counts_file = generate_counts.counts_file,
                reference = reference,
                reference_index = reference_index,
                docker = docker
        }
        call locityper_genotype {
            input: 
                sample_id = sample_ids[scatter_index],
                input_fq1 = samtools_collate.reads_1,
                input_fq2 = samtools_collate.reads_2,
                db_targz = generate_db.db_tar,
                preprocess_tar = locityper_preprocess.preprocess_tar,
                docker = docker
        }
    }

    output {
        Array[File] results = locityper_genotype.genotype_tar 
    }
}

task samtools_collate {
    input {
        File bam_file
        File bam_index
        String docker
    }

    Int disk_size = ceil(size(bam_file, "GiB") + 100)
    String base_filename = basename(basename(bam_file, ".bam"), ".cram")
    String output_fq1 = base_filename + "_1.fastq.gz"
    String output_fq2 = base_filename + "_2.fastq.gz"

    command <<<
        set -euo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        samtools collate -@ ${nthreads} -u -O ~{bam_file} | samtools fastq -@ ${nthreads} -1 ~{output_fq1} -2 ~{output_fq2} -0 /dev/null -s /dev/null -n
    >>>

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File reads_1 = output_fq1
        File reads_2 = output_fq2
    }
}

task locityper_preprocess {
    input {
        File input_fq1
        File input_fq2
        File counts_file
        File reference
        File reference_index
        String docker
    }

    Int disk_size = ceil(size(input_fq1, "GiB") + size(input_fq2, "GiB")) * 10 + 50
    String output_tar = "locityper_prepoc.tar.gz"

    command <<<
        set -euo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        mkdir -p locityper_prepoc
        locityper preproc -i ~{input_fq1} ~{input_fq2} \
            -j ~{counts_file} \
            -@ ${nthreads} \
            --technology illumina \
            -r ~{reference} \
            -o locityper_prepoc
        tar -czf ~{output_tar} locityper_prepoc
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File preprocess_tar = output_tar
    }
}

task locityper_genotype {
    input {
        String sample_id
        File input_fq1
        File input_fq2
        File preprocess_tar
        File db_targz
        String docker
    }

    Int disk_size = ceil(size(input_fq1, "GiB") + size(input_fq2, "GiB")) * 10 + 50
    String output_tar = sample_id + ".tar.gz"

    command <<<
        set -euo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        
        mkdir -p locityper_prepoc
        mkdir -p db

        tar -xvzf ~{preprocess_tar} -C locityper_prepoc 
        tar -xvzf ~{db_targz} -C db 
        mkdir -p ~{sample_id}

        locityper genotype -i ~{input_fq1} ~{input_fq2} \
            -d db \
            -p locityper_prepoc
            -@ ${nthreads} \
            --debug 2 \
            -o ~{sample_id}
        tar -czf ~{output_tar} ~{sample_id}
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File genotype_tar = output_tar
    }
}

task generate_db {
    input {
        File reference
        File reference_index
        File counts_jf
        String locus_name
        String locus_coordinates
        File alleles_fa
        String docker
    }

    Int disk_size = 10
    String output_tar = locus_name + ".db.tar.gz"

    command <<<
        set -euo pipefail

        locityper add -d ~{locus_name}.db \
            -r ~{reference} \
            -j ~{counts_jf} \
            -l ~{locus_name} ~{locus_coordinates} ~{alleles_fa}

        tar -czf output_tar ~{locus_name}.db
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File db_tar = output_tar
    }
}


task generate_counts {
    input {
        File reference
        String docker
    }

    Int disk_size = 10
    String output_file = "counts.jf"

    command <<<
        set -euo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        jellyfish count --canonical \
            --lower-count 2 \
            --out-counter-len 2 \
            --mer-len 25 \
            --threads ${nthreads} \
            --size 3G \
            --output ~{output_file} \
            ${reference}
    >>>

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File counts_file = output_file
    }
}

