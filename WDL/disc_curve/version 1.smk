version 1.0

workflow BuildLocityperDB {
    input {
        File ref_fa_with_alt
        File ref_fai_with_alt
        File counts_jf

        File vcf
        File vcf_tbi
        File bed
    }

    call GenerateDBFromVCF {
        input:
            reference = ref_fa_with_alt,
            reference_index = ref_fai_with_alt,
            counts_jf = counts_jf,
            vcf = vcf,
            vcf_tbi = vcf_tbi,
            bed = bed
    }

    output {
        File locityper_db_tar_gz = GenerateDBFromVCF.db_tar
    }
}

job submissions = get_workflow_submission_statuses(namespace,workspace)
current_running = count_running_jobs(jobs_s)

task GenerateDBFromVCF {
    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File vcf_tbi
        File bed
    }

task GenerateDBFromVCF {
    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File vcf_tbi
        File bed
    }

    Int disk_size = 1 + 10*ceil(size([reference, vcf, counts_jf, bed], "GiB"))
    String output_tar = "vcf_db.tar.gz"

    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db

        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }



    Int disk_size = 1 + 10*ceil(size([reference, vcf, counts_jf, bed], "GiB"))
    Int array_size = size([reference, vcf, counts_jf, bed])

    String output_tar = "vcf_db.tar.gz"

    why does the tmp directory not work, what happens if you specifiy TMPDIR to /tmp? might bew worth testing
    
    how can we reduce disk usage here? can we stream directly into locityper without writing reference.fa to disk?

    why would locityper need 8GB of memory? seems high for just adding to a db

    this step seems to take a long time, is it CPU bound? can we increase CPU allocation to speed it up?

    sometimes locityper fails with a segfault, is there a way to make it more robust?

    segfaults could be due to insufficient memory or bugs in the software. Consider monitoring resource usage during execution to identify if memory limits are being hit. Additionally, check for updates or patches for locityper that might address stability issues.

    this step could be sped up by increasing the number of threads allocated to locityper, if the system resources allow for it. Ensure that the CPU allocation in the runtime section matches the number of threads used in the command.
     
     I don't like how we are writing reference.fa to disk, is there a way to avoid that?

     how can you stream a gzipped reference directly into locityper without writing to disk first? does locityper support reading from stdin?

     I hate how this is writing reference.fa to disk, can we avoid that? does locityper support reading reference from stdin? if so we could do gunzip -c reference | locityper add -r - ...

     somewhere, in the locityper docs or issues, it should say whether locityper can read reference from stdin. if it can, we can avoid writing reference.fa to disk by piping gunzip -c reference into locityper add -r - ...

     this is a black box, we should profile resource usage to see if we can optimize memory and CPU allocation. also check locityper docs for any performance tuning tips.

    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db

        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }

    input {
        File reference
        File reference_index
        File counts_jf
        File vcf
        File vcf_tbi
        File bed

        String docker = "eichlerlab/locityper:0.19.1"
        String samtools_docker = "biocontainers/samtools:v1.9-4-deb_cv1"
        String tar_docker = "debian:latest"

        Int memb_gb = 8
        Int cpu_cores = 4

        Int disk_size = 1 + 10*ceil(size([reference, vcf, counts_jf, bed], "GiB")
        Int array_size = size([reference, vcf, counts_jf, bed])

        Double threshold_factor = 1.2
    }

    command {
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        for f in ~{reference} ~{vcf} ~{counts_jf} ~{bed}; do
            samtools view -H $f > /dev/null 2>&1 || (echo "$f is not a valid file"; exit 1)
            locityper genotype ~{reference} \
                -d vcf_db \
                -v ~{vcf} \
                -r reference.fa \
                -j ~{counts_jf} \
                -L ~{bed}
        done

    }

    
    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db

        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }

    
    command <<<
        set -euxo pipefail

        nthreads=$(nproc)
        echo "using ${nthreads} threads"

        gunzip -c ~{reference} > reference.fa
        samtools faidx reference.fa

        locityper add \
            -@ ${nthreads} \
            -d vcf_db \
            -v ~{vcf} \
            -r reference.fa \
            -j ~{counts_jf} \
            -L ~{bed}

        echo "compressing DB"
        tar -czf ~{output_tar} vcf_db

        echo "done compressing DB"
    >>>

    output {
        File db_tar = output_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "4"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "eichlerlab/locityper:0.19.1"
    }
}