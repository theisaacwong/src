#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(plyr)
library(GenomicRanges)
library(IRanges)
library(stringr)

INPUT_PAF_FOFN <- args[1]
INDEX_FAI <- "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta.fai"

OUTPUT_FILE <- args[2]

df_index <- read.table(INDEX_FAI, header=FALSE, comment.char="")
colnames(df_index) <- c("chr", "length", "V3", "V4", "V5")

df_coverage <- lapply(INPUT_PAF_FOFN %>% readLines, function(INPUT_PAF){
  print(INPUT_PAF)
  df0 <- fread(INPUT_PAF, header=FALSE, fill=TRUE) %>% as.data.frame()
  df0$V23 <- NULL
  df0$V24 <- NULL
  df0$V25 <- NULL

  df1 <- lapply(df0$V1 %>% u, function(x){
    temp1 <- df0[df0$V1==x, ]
    AS_SCORE <- temp1$V15 %>% str_remove("AS:i:") %>% as.numeric()
    return(temp1[which(AS_SCORE==max(AS_SCORE)), ])
  }) %>% do.call(rbind, .) %>% as.data.frame()

  gr1 <- GRanges(seqnames = df1$V6,
                 ranges=IRanges(start=df1$V8,
                                end=df1$V9))
  gr1 <- sortSeqlevels(gr1)
  gr1 <- sort(gr1)

  df2 <- reduce(gr1) %>% as.data.frame(stringsAsFactors=FALSE)
  df2$seqnames <- df2$seqnames %>% as.character()
  df2$length <- df2$end - df2$start

  MY_CHROMOSOMES <- sort(unique(df2$seqnames))
  coverage_by_chrom <- lapply(MY_CHROMOSOMES, function(x){
    curr_sum <- sum(df2$length[df2$seqnames == x])
    return(curr_sum)
  }) %>% as.data.frame()
  names(coverage_by_chrom) <- MY_CHROMOSOMES
  coverage_by_chrom$SAMPLE <- INPUT_PAF # %>% str_extract("....._(fa|mo|p1|s1|s2)")
  return(coverage_by_chrom[, c("SAMPLE", MY_CHROMOSOMES)])

}) %>% do.call(rbind.fill, .)

df_coverage$chrM <- NULL

df_percent_coverage <- colnames(df_coverage)[-c(1)] %>% lapply(function(x){
  curr_chrom <- x
  curr_len <- df_index$length[df_index$chr==curr_chrom]
  return(df_coverage[x] / curr_len)
}) %>% as.data.frame()
df_percent_coverage$SAMPLE <- df_coverage$SAMPLE
df_percent_coverage <- df_percent_coverage[colnames(df_coverage)]

write.table(df_percent_coverage, OUTPUT_FILE, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



