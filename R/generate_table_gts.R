#!/usr/bin/env Rscript

source("/net/eichler/vol28/home/iwong1/nobackups/src/R/my_libraries.R")

library(parallel)

MY_CORES <- 120

print(Sys.time())
print("reading in GTs")

df00 <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/ehr_gts_mans/table.txt", sep="\t") %>% as.data.frame()
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

print(Sys.time())
print("reading in files")

files <- list.files("/net/eichler/vol28/home/iwong1/nobackups/aou/dl_parallel", full.names = TRUE, recursive = TRUE, pattern = ".csv")
df01 <- mclapply(files, function(x){
  df_temp <- fread(x) %>% as.data.frame()
  return(df_temp)
}, mc.cores=MY_CORES) %>% do.call(rbind, .)

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

# df01$ltgt_to_gt[df01$ltgt_to_gt=="0|1"] <- "1|0"

print(Sys.time())
print("converting IDs")

df_map <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/ehr_gts_mans/004-isaac-first-batch_SVs_of_interest_013025_CDS_name_map.tsv") %>% as.data.frame()
codeToName <- HashMap$new()
codeToName$populate(df_map, "new_name", "V1")

df01$converted_names <- mclapply(seq_down(df01), function(x){
  rvals <- lapply(unlist(df01$genotype[x] %>% str_split(",")), function(gt){
    rval <- "GRCh38"
    if(gt!="GRCh38"){
      # which_hap  <- gt %>% str_extract("(?<=\\.).$") %>% as.numeric()
      which_samp <- gt %>% str_extract(".*(?=..$)")
      if(str_detect(which_samp, ":")){
        which_samp <- which_samp %>% str_extract("(?<=:).*?$")
      }
      # print(paste0(gt,  ": ", which_samp, " ::: ", which_hap))
      rval <- codeToName$get(which_samp)
      if(str_detect(rval, ":")){
        temp1 <- unlist(str_split(rval, ":"))
        rval <- temp1[length(temp1)]
      }
    }
    return(rval)
  }) %>% unlist %>% paste0(collapse = ",")
  
  return(rvals)
}, mc.cores=MY_CORES) %>% unlist

df01$genotype <- df01$converted_names
df01$converted_names <- NULL

df_panel <- fread("/net/eichler/vol28/home/iwong1/nobackups/aou/ehr_gts_mans/Locityper_SV_Panel_1KG_TestRun.tsv") %>% as.data.frame()

df01$gene <- mclapply(df01$locus, function(x){
  temp1 <- df_panel$GENE[df_panel$SVID==x]
  if(length(temp1)==0){temp1 <- "None"}
  paste0(temp1, collapse=",")
}, mc.cores=MY_CORES) %>% unlist

print(Sys.time())
print("saving")

save(df01, file="df01.rda")

print(Sys.time())
print("writing file")

write.table(df01, "locityper_ehr_gt_summary.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
