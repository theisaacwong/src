#!/usr/bin/env Rscript

source("/net/eichler/vol28/home/iwong1/nobackups/src/R/my_libraries.R")

library(parallel)

MY_CORES <- 120

print(Sys.time())
print("reading in GTs")

df00 <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/batch4/terra_outputs/table_short.tsv", sep="\t") %>% as.data.frame()
df00$ID <- lapply(df00$ID, function(x){
  unlist(str_split(x, ";"))[1]
}) %>% unlist
variantToSamples <- HashMap$new()
temp1 <- lapply(seq_down(df00), function(x){
  curr_var <- df00$ID[x]
  variantToSamples$put(curr_var, HashMap$new())
  temp2 <- lapply(4:ncol(df00), function(y){
    curr_sample <- colnames(df00)[y] %>% as.character()
    curr_genotype <- df00[x,y]
    variantToSamples$get(curr_var)$put(curr_sample, curr_genotype)
  })
})

print(Sys.time())
print("reading in files")
#files <- list.files("/net/eichler/vol28/home/iwong1/nobackups/aou/batch4/terra_outputs/csv_outputs", full.names = TRUE, recursive = TRUE, pattern = "gts.filtered.csv")
#df01 <- mclapply(files, function(x){
#  df_temp <- fread(x, fill=TRUE) %>% as.data.frame()
#  df_temp <- df_temp[df_temp$genotype!="*", ]
#  return(df_temp)
#}, mc.cores=MY_CORES)
#save(df01, file="/net/eichler/vol28/home/iwong1/nobackups/aou/batch4/terra_outputs/df01_premerge.rda")
#load("/net/eichler/vol28/home/iwong1/nobackups/aou/batch3/terra_outputs/df01_premerge.rda")
#to_keep <- lapply(df01, ncol) %>% unlist
#df01 <- df01[to_keep==8] %>% do.call(rbind, .)
#fwrite(df01, "df01_post_merge_batch4.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#save(df01, file="df01_postmerge.rda")
#load("df01_postmerge.rda")

df01 <- fread("df01_post_merge_batch4.tsv") %>% as.data.frame()

print(Sys.time())
print("converting GTs")

multiallelics <- c("MUC1", "MUC2", "MUC3A", "MUC4", "MUC5AC", "MUC5B", "MUC6", "MUC7", "MUC12", "MUC17", "MUC20", "MUC21", "MUC22", "DMBT1", "FCGR", "LPA", "HPR", "TPSAB1")

df01$ltgt_to_gt <- mclapply(seq_down(df01), function(x){
  if(df01$locus[x] %in% multiallelics){ return(df01$genotype[x]) }
  rvals <- lapply(unlist(df01$genotype[x] %>% str_split(",")), function(gt){
    rval <- gt
    if(gt=="GRCh38"){
      rval <- "0"
    } else {
      which_hap  <- gt %>% str_extract("(?<=\\.).$") %>% as.numeric()
      which_samp <- gt %>% str_extract(".*(?=..$)") %>% as.character()
      # print(paste0(gt,  ": ", which_samp, " ::: ", which_hap))
      if(variantToSamples$containsKey(df01$locus[x])){
        if(variantToSamples$get(df01$locus[x])$containsKey(which_samp))
        {
          rval <- unlist(str_split(variantToSamples$get(df01$locus[x])$get(which_samp), "\\|"))[which_hap]
        }
      }
    }
    return(rval)
  }) %>% unlist %>% paste0(collapse = "|")

  return(rvals)
}, mc.cores=MY_CORES) %>% unlist

df01$ltgt_to_gt[df01$ltgt_to_gt=="0|1"] <- "1|0"

#print(Sys.time())
#print("saving")
#save(df01, file="df01.rda")

print(Sys.time())
print("writing file")

fwrite(df01, "locityper_summary_20260604.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
