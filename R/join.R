#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output <- args[length(args)]

df_temp <- 0
IS_HEADER=FALSE
df_final <- read.table(args[1], sep="\t", header=FALSE)

if(df_final$V1[1]=="SAMPLE"){
	IS_HEADER <- TRUE
	df_final <- read.table(args[1], sep="\t", header=IS_HEADER)
}

n_last <- length(args) - 1

for(i in 2:n_last){
	df_temp <- read.table(args[i], sep="\t", header=IS_HEADER)

	if(IS_HEADER){
		df_final <- merge(df_final, df_temp, by = "SAMPLE")
	} else {
		df_final <- merge(df_final, df_temp, by = "V1")
	}
}

write.table(df_final, output, sep="\t", col.names=IS_HEADER, row.names=FALSE, quote=FALSE)


