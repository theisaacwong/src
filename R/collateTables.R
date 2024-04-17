#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if(any(c("-h", "--help") %in% args)){
  print("Rscript collateTables.R [directory] [file_pattern] [output]")
  quit("no")
}

my_dir <- args[1]
my_pattern <- args[2]
my_output <- args[3]


my_files <- list.files(path = my_dir,
                       pattern = my_pattern,
                       full.names = TRUE,
                       recursive = TRUE)

df0 <- lapply(my_files, function(x){

  df_temp <- fread(x) %>% as.data.frame()
  df_temp$SOURCE_FILE <- x
  return(df_temp)

}) %>% do.call(rbind, .)

fwrite(df0, file = my_output, quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
