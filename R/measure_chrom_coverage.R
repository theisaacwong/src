#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


INPUT_BED_FOFN <- args[1] # must be merged
INDEX_FAI <- args[2]

OUTPUT_FILE <- args[3]

INDEX_FAI <- "C:/docs/autism/coverage/T2T-CHM13v2.fasta.fai"
INPUT_BED_FOFN <- "C:/docs/autism/coverage/my_beds.txt"


df_index <- read.table(INDEX_FAI, header=FALSE)
colnames(df_index) <- c("chr", "length", "V3", "V4", "V5")



list_by_bed <- lapply(INPUT_BED_FOFN %>% readLines, function(INPUT_BED){
  print(INPUT_BED)

  if(Sys.info()["sysname"]=="Windows"){
    df_bed <- read.table(text = system2("wsl", c("bedtools", "merge", "-i", INPUT_BED), stdout = TRUE), header=FALSE, comment.char = "")
  }else{
    df_bed <- read.table(text = system2("bedtools", c("merge", "-i", INPUT_BED), stdout = TRUE), header=FALSE, comment.char = "")
  }
  colnames(df_bed)[1:3] <- c("chr", "start", "end")

  df_bed$length <- df_bed$end - df_bed$start

  MY_CHROMOSOMES <- sort(unique(df_bed$chr))
  coverage_by_chrom <- lapply(MY_CHROMOSOMES, function(x){
    curr_sum <- sum(df_bed$length[df_bed$chr == x])
    return(curr_sum)
  }) %>% as.data.frame()
  names(coverage_by_chrom) <- MY_CHROMOSOMES
  coverage_by_chrom$SAMPLE <- INPUT_BED %>% str_extract("....._(fa|mo|p1|s1|s2)")
  rownames(coverage_by_chrom) <- INPUT_BED %>% str_extract("....._(fa|mo|p1|s1|s2)")
  return(coverage_by_chrom[, c("SAMPLE", MY_CHROMOSOMES)])
})
df_coverage <- do.call(rbind, lapply(list_by_bed, function(x) x[match(names(list_by_bed[[1]]), names(x))]))

df_percent_coverage <- colnames(df_coverage)[-c(1)] %>% lapply(function(x){
  curr_chrom <- x
  curr_len <- df_index$length[df_index$chr==curr_chrom]
  return(df_coverage[x] / curr_len)
}) %>% as.data.frame()
df_percent_coverage$SAMPLE <- df_coverage$SAMPLE
df_percent_coverage <- df_percent_coverage[colnames(df_coverage)]

write.table(df_percent_coverage, OUTPUT_FILE, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



