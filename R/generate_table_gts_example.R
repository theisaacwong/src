#!/usr/bin/env Rscript

source("/net/eichler/vol28/home/iwong1/nobackups/src/R/my_libraries.R")

library(parallel)

MY_CORES <- 120

print(Sys.time())
print("reading in GTs")

# this is a table format of VCF where the rows are variants, columns are samples, and fields are the genotype values,
# bcftools query -l example.vcf.gz| tr "\n" "\t" > table.txt ; bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' example.vcf.gz >> table.txt
# then do some manual editing for the column names should get you a similar file
df00 <- fread("table.txt", sep="\t") %>% as.data.frame()
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
# this is just a list of files, something like:
#       find -type f | grep .gts.filtered.csv > gts.filtered.csv.files
files <- readLines("gts.filtered.csv.files")
#files <- list.files("/path/to/files", full.names = TRUE, recursive = TRUE, pattern = "gts.filtered.csv") # this is slower than the above example
df01 <- mclapply(files, function(x){
  df_temp <- fread(x, fill=TRUE) %>% as.data.frame()
  df_temp <- df_temp[df_temp$genotype!="*", ]
  return(df_temp)
}, mc.cores=MY_CORES)
save(df01, file="df01_premerge.rda")
#load("df01_premerge.rda")
to_keep <- lapply(df01, ncol) %>% unlist
df01 <- df01[to_keep==8] %>% do.call(rbind, .)
fwrite(df01, "df01_post_merge_batch5.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#df01 <- fread("df01_post_merge_batch5.tsv") %>% as.data.frame()

print(Sys.time())
print("converting GTs")


df01$ltgt_to_gt <- mclapply(seq_down(df01), function(x){
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


print(Sys.time())
print("writing file")

fwrite(df01, "locityper_summary.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
