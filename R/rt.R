#!/usr/bin/env Rscript
library(data.table)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
df0 <- fread(input) %>% as.data.frame() %>% t %>% as.data.frame()
colnames(df0) <- c("")
print(df0)
