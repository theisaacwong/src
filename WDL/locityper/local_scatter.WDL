version 1.0

import "src/WDL/locityper/ValidateVariants.WDL" as local_locityper

workflow local_scatter {

    input {
        Array[File] crams
        Array[File] crais
        
        String sample_id

        File ref_fa_uncompressed
        File ref_fai_uncompressed

        File counts_jf
        File bed

        File? sample_map

        File locityper_db_tar_gz

        Int locityper_n_cpu = 32
        Int locityper_max_gts = 100000
    }

    scatter (scatter_index in range(length(crams))) {
        call local_locityper {
            input:
                cram = crams[scatter_index],
                crai = crais[scatter_index],
                sample_id = sample_id,
                ref_fa_uncompressed = ref_fa_uncompressed,
                ref_fai_uncompressed = ref_fai_uncompressed,
                counts_jf = counts_jf,
                bed = bed,
                sample_map = sample_map,
                locityper_db_tar_gz = locityper_db_tar_gz,
                locityper_n_cpu = locityper_n_cpu,
                locityper_max_gts = locityper_max_gts
        }
    }

    output {
        Array[File] results_tar_gzs = local_locityper.results_tar_gz
        Array[File] summary_csvs = local_locityper.summary_csv
    }

}