#!/usr/bin/env Rscript

source("/net/eichler/vol28/home/iwong1/nobackups/src/R/my_libraries.R")

library(parallel)

MY_CORES <- 120

print(Sys.time())
print("reading in GTs")

df00 <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/batch2/aou_1kg_tier1.table.txt", sep="\t") %>% as.data.frame()
df00$ID <- lapply(df00$ID, function(x){
  unlist(str_split(x, ";"))[1]
}) %>% unlist
variantToSamples <- HashMap$new()
temp1 <- lapply(seq_down(df00), function(x){
  curr_var <- df00$ID[x]
  variantToSamples$put(curr_var, HashMap$new())
  temp2 <- lapply(4:ncol(df00), function(y){
    curr_sample <- colnames(df00)[y]
    curr_genotype <- df00[x,y]
    variantToSamples$get(curr_var)$put(curr_sample, curr_genotype)
  })
})

df00 <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/batch2/temp/svs_for_locityper_batch2_julie_yang_050126.table.txt", sep="\t") %>% as.data.frame()
df00$ID <- lapply(df00$ID, function(x){
  unlist(str_split(x, ";"))[1]
}) %>% unlist
temp1 <- lapply(seq_down(df00), function(x){
  curr_var <- df00$ID[x]
  variantToSamples$put(curr_var, HashMap$new())
  temp2 <- lapply(4:ncol(df00), function(y){
    curr_sample <- colnames(df00)[y]
    curr_genotype <- df00[x,y]
    variantToSamples$get(curr_var)$put(curr_sample, curr_genotype)
  })
})

print(Sys.time())
print("reading in files")
files <- list.files("/net/eichler/vol28/home/iwong1/nobackups/aou/batch2/terra_outputs/results", full.names = TRUE, recursive = TRUE, pattern = "gts.filtered.csv")
df01 <- mclapply(files, function(x){
  df_temp <- fread(x, fill=TRUE) %>% as.data.frame()
  df_temp <- df_temp[df_temp$genotype!="*", ]
  return(df_temp)
}, mc.cores=MY_CORES)
save(df01, file="df01_premerge.rda")

to_keep <- lapply(df01, ncol) %>% unlist
save(df01[[to_keep != 8]], "errors.rda")

df01 <- df01[to_keep==8] %>% do.call(rbind, .)

print(Sys.time())
print("converting GTs")

df01$ltgt_to_gt <- mclapply(seq_down(df01), function(x){
  rvals <- lapply(unlist(df01$genotype[x] %>% str_split(",")), function(gt){
    rval <- "ooblets"
    if(gt=="GRCh38"){
      rval <- "0"
    } else {
      which_hap  <- gt %>% str_extract("(?<=\\.).$") %>% as.numeric()
      which_samp <- gt %>% str_extract(".*(?=..$)")
      # print(paste0(gt,  ": ", which_samp, " ::: ", which_hap))
      rval <- unlist(str_split(variantToSamples$get(df01$locus[x])$get(which_samp), "\\|"))[which_hap]
    }
    return(rval)
  }) %>% unlist %>% paste0(collapse = "|")

  return(rvals)
}, mc.cores=MY_CORES) %>% unlist

df01$ltgt_to_gt[df01$ltgt_to_gt=="0|1"] <- "1|0"

print(Sys.time())
print("saving")

save(df01, file="df01.rda")

print(Sys.time())
print("writing file")

write.table(df01, "locityper_summary_20260515.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
