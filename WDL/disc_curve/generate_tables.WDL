version 1.0

workflow discovery_curve {

    input {
        String docker
        Array[File] input_h1_files
        Array[File] input_h2_files
        File sv_vcf
        File sv_vcf_index
        File sample_file
        File sample_order
        File reference_index
        File python_script
        

        String name

        String n_cpu = 4
        String mem_gb = 8
        Int disk_size
    }


    call create_tables {
        input:
            input_h1_files = input_h1_files,
            input_h2_files = input_h2_files,
            sv_vcf = sv_vcf,
            sv_vcf_index = sv_vcf_index,
            sample_file = sample_file,
            sample_order = sample_order,
            reference_index = reference_index,
            python_script = python_script,
            name = name,
            docker = docker,
            mem_gb = mem_gb,
            n_cpu = n_cpu,
            disk_size = disk_size
    }
    
    output {
        File output_bed = create_tables.output_bed
        File output_tsv = create_tables.output_tsv
        File output_hap = create_tables.output_hap
        File output_gts = create_tables.output_gts
        File output_counts = create_tables.output_counts
    }
}

task create_tables {
    input {
        Array[File] input_h1_files
        Array[File] input_h2_files
        File sv_vcf
        File sv_vcf_index
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

        find /cromwell_root/ | grep "_callable_regions_h._500.bed.gz" | xargs -I {} mv {} /callable_files/

        /opt/conda/envs/audano-curve/bin/python ~{python_script} \
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

