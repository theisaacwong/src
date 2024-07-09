#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


SAMPLE_ID <- args[1]



if(Sys.info()["sysname"]=="Windows"){
  vcf_header <- system2("wsl", args = c("zcat", file, ' | grep -m 1 "^#CHROM"'), stdout = TRUE) %>% str_split("\t") %>% unlist # read in just the column header
}else{
  vcf_header <- system2("zcat", args = c(file, ' | grep -m 1 "^#CHROM"'), stdout = TRUE) %>% str_split("\t") %>% unlist # read in just the column header
}
