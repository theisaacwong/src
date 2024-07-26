#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(plyr)

INPUT_PAF_FOFN <- args[1] # will need to be merged
INDEX_FAI <- "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta.fai"

OUTPUT_FILE <- args[2]

df_index <- read.table(INDEX_FAI, header=FALSE, comment.char="")
colnames(df_index) <- c("chr", "length", "V3", "V4", "V5")

list_by_bed <- lapply(INPUT_PAF_FOFN %>% readLines, function(INPUT_PAF){


  if(Sys.info()["sysname"]=="Windows"){
    df_bed <- read.table(text = system2("wsl", c("bedtools", "merge", "-i", "TEMP.TSV"), stdout = TRUE), header=FALSE, comment.char = "")
  }else{
    print(INPUT_PAF)
    system2("cut", c("-f6,8,9", INPUT_PAF, ">", "TEMP_0.TSV"))
    system2("bedtools", c("sort",  "-i", "TEMP_0.TSV", ">", "TEMP_1.tsv"))
    system2("bedtools", c("merge", "-i", "TEMP_1.tsv", ">", "TEMP_2.tsv"))
    df_bed <- read.table("TEMP_2.tsv", header=FALSE, comment.char = "")
  }
  colnames(df_bed)[1:3] <- c("chr", "start", "end")

  df_bed$length <- df_bed$end - df_bed$start

  MY_CHROMOSOMES <- sort(unique(df_bed$chr))
  coverage_by_chrom <- lapply(MY_CHROMOSOMES, function(x){
    curr_sum <- sum(df_bed$length[df_bed$chr == x])
    return(curr_sum)
  }) %>% as.data.frame()
  names(coverage_by_chrom) <- MY_CHROMOSOMES
  coverage_by_chrom$SAMPLE <- INPUT_PAF  # INPUT_BED %>% str_extract("....._(fa|mo|p1|s1|s2)")
  #  rownames(coverage_by_chrom) <- # INPUT_BED %>% str_extract("....._(fa|mo|p1|s1|s2)")
  return(coverage_by_chrom[, c("SAMPLE", MY_CHROMOSOMES)])
})


df_coverage <- do.call(rbind.fill, list_by_bed)
df_coverage$chrM <- NULL

df_percent_coverage <- colnames(df_coverage)[-c(1)] %>% lapply(function(x){
  curr_chrom <- x
  curr_len <- df_index$length[df_index$chr==curr_chrom]
  return(df_coverage[x] / curr_len)
}) %>% as.data.frame()
df_percent_coverage$SAMPLE <- df_coverage$SAMPLE
df_percent_coverage <- df_percent_coverage[colnames(df_coverage)]

write.table(df_percent_coverage, OUTPUT_FILE, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



