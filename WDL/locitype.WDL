version 1.0

workflow runLocityper {
    input {
        String docker
        File reference
        File reference_index
        File counts_jf
        String locus_name
        String locus_coordinates
        File alleles_fa
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
    
    output {
        File results = generate_db.db_tar 
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
    String output_tar = "db.tar.gz"

    command <<<
        set -euxo pipefail

        mkdir -p db

        locityper add -d db \
            -r ~{reference} \
            -j ~{counts_jf} \
            -l ~{locus_name} ~{locus_coordinates} ~{alleles_fa}
        
        echo "compressing DB"
        tar -czf output_tar db
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
