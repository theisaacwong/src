version 1.0

workflow runLocityper {
    input {
        Array[String] sample_ids
        Array[File] fastq_r1
        Array[File] fastq_r2
        String docker
        File reference
        File reference_index
        File counts_jf
        File db_path
        String locus_name
        String locus_coordinates
        File alleles_fa
        File weights_file


        String locityper_n_cpu = 4
        String locityper_mem_gb = 8
    }
   

    call generate_db {
    input:
        reference = reference,
        reference_index = reference_index,
        counts_jf = counts_jf,
        locus_name = locus_name,
        locus_coordinates = locus_coordinates,
        alleles_fa = alleles_fa,
        docker = docker
    }
    
    scatter (scatter_index in range(length(sample_ids))) {
        call locityper_preprocess_and_genotype {
            input:
                sample_id = sample_ids[scatter_index],
                input_fq1 = fastq_r1[scatter_index],
                input_fq2 = fastq_r2[scatter_index],
                db_targz = generate_db.db_tar,
                counts_file = counts_jf,
                reference = reference,
                weights_file = weights_file,
                locus_name = locus_name,
                reference_index = reference_index,
                docker = docker,
                locityper_n_cpu = locityper_n_cpu,
                locityper_mem_gb = locityper_mem_gb

        }
    }

    output {
        Array[File] results = locityper_preprocess_and_genotype.genotype_tar
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
    String output_tar = locus_name + "db.tar.gz"

    command <<<
        set -euxo pipefail

        locityper add -d ~{locus_name}.db \
            -r ~{reference} \
            -j ~{counts_jf} \
            -l ~{locus_name} ~{locus_coordinates} ~{alleles_fa}
        
        echo "compressing DB"
        tar -czf ~{output_tar} ~{locus_name}.db
        echo "done compressing DB"

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

task locityper_preprocess_and_genotype {
    input {
        File input_fq1
        File input_fq2
        File counts_file
        File reference
        File reference_index
        File weights_file
        String docker
        String locityper_n_cpu
        String locityper_mem_gb
        File db_targz
        String sample_id
        String locus_name
    }

    Int disk_size = ceil(size(input_fq1, "GiB") + size(input_fq2, "GiB")) + 80
    String output_tar = sample_id + "." + locus_name + ".tar.gz"

    command <<<
        set -euxo pipefail
        nthreads=$(nproc)
        echo "using ${nthreads} threads"
        mkdir -p locityper_prepoc
        locityper preproc -i ~{input_fq1} ~{input_fq2} \
            -j ~{counts_file} \
            -@ ${nthreads} \
            --technology illumina \
            -r ~{reference} \
            -o locityper_prepoc

        mkdir -p db
        tar --strip-components 1 -C db -xvzf ~{db_targz}
        mkdir -p out_dir

        echo -e {locus_name}"\t"{weights_file} > weights.tsv

        locityper genotype -i ~{input_fq1} ~{input_fq2} \
            -d db \
            -p locityper_prepoc \
            -@ ${nthreads} \
            --reg-weights weights.tsv \
            --debug 2 \
            -o out_dir
        tar -czf ~{output_tar} out_dir

    >>>

    runtime {
        memory: locityper_mem_gb + " GB"
        cpu: locityper_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        docker: docker
    }

    output {
        File genotype_tar = output_tar
    }
}


