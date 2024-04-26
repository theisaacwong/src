#!/usr/env/bin R
## Extract private variants ##
## load packages ##
require(readr)
require(tidyr)
require(dplyr)
require(parallel)
require(data.table)
library(stringr)

# command line arguments
args = commandArgs(trailingOnly = T)
ped_file = args[1]
counts_file = args[2]
file = args[3]
het_ac = as.numeric(args[4])
hom_ac = as.numeric(args[5])

seq_down <- function(x){ 1:nrow(x)}
MIN_DP_THRESHOLD <- 20

## define functions ##
# extract and format rare variants from a family
extract_rare = function(file){
  famid = gsub("inheritance/", "", gsub(".inheritance.vcf.gz", "", file))
  vcf_header <- system2("zcat", args = c(file, ' | grep -m 1 "^#CHROM"'), stdout = TRUE) %>% str_split("\t") %>% unlist # read in just the column header
  col_index_fa <- grep("\\.fa$", vcf_header) # determine which column is for fa/mo
  col_index_mo <- grep("\\.mo$", vcf_header)
  if(any(col_index_fa %in% col_index_mo) | length(col_index_fa)!=1 | length(col_index_mo)!=1){ # if regex encountered a conflict, stop
    stop(paste0("Error in column names. \n\tFA: ", vcf_header[col_index_fa], "\n\tMO: ", vcf_header[col_index_mo]))}
  COL_NAME_TEMPLATE <- c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  col_names_vcf <- vcf_header
  col_names_vcf[col_index_fa] <- "FA" # relabel original column names
  col_names_vcf[col_index_mo] <- "MO"
  col_names_vcf[1:9] <- COL_NAME_TEMPLATE
  col_classes <- list(integer=2, # fread is faster if you tell it what data types are in each column
                      numeric=6,
                      character=c(1,3:5,6:length(col_names_vcf)))

  skip_lines <- system2("zcat", args = c(file, " | sed '/^#CHROM/q' | wc -l"), stdout = TRUE) %>% as.numeric() # skip header lines
  raw = suppressWarnings(
    fread(file,
          sep = "\t",
          header=FALSE,
          skip = skip_lines,
          colClasses = col_classes,
          col.names = col_names_vcf) %>% as.data.frame())

	info = c(famid, nrow(raw))
	info_cols = names(raw)[1:8]

	pro_sex = paste(unique(ped$V6[which(grepl("p", ped$V2) & ped$V1 == famid)]), collapse = ",")
	if(any(grepl("s", ped$V2) & ped$V1 == famid)){
		sib_sex = paste(unique(ped$V6[which(grepl("s", ped$V2) & ped$V1 == famid)]), collapse =",")
	}else{
		sib_sex = "NA"
	}


	# for each row in the FORMAT column, get a list of which index corresponds to each field
	my_FORMAT <- raw$FORMAT %>% str_split(":")
	index_GT <- lapply(my_FORMAT, function(x){which(x=="GT")}) %>% unlist
	index_AD <- lapply(my_FORMAT, function(x){which(x=="AD")}) %>% unlist
	index_DP <- lapply(my_FORMAT, function(x){which(x=="DP")}) %>% unlist


	# get a vector of indexes for columns corresponding to fa, mo, p1, p2, etc. written more verbose to scale dynamically with family size
	which_columns_to_split <- which(colnames(raw) %in% c("FA", "MO") | grepl("\\.p[0-9]$", colnames(raw)))
	family_columns_to_cbind <- lapply(which_columns_to_split, function(COL_INDEX){
	  split_member <- raw[, COL_INDEX] %>% str_split(":")
	  # identify which columns have more FORMAT fields than the FA column has fields, we won't be able to get good AD or DP for for these
	  bool_useable <- (my_FORMAT %>% lapply(length) %>% unlist)==(split_member %>% lapply(length) %>% unlist)
	  columns_member <- lapply(seq_down(raw), function(x){
	    if(bool_useable[x]==FALSE){
	      if(split_member[[x]][1] %>% grepl("\\/|\\|")){ # check to see if the first field contains GT information
	        return(c(split_member[[x]][1], 0,0))
	      } else {
	        return(c("./.", 0, 0))
	      }
	    }

	    # return the corresponding row from current family member, and the corresponding column for each corresponding field
	    return(c(split_member[[x]][index_GT[x]],
	             split_member[[x]][index_AD[x]],
	             split_member[[x]][index_DP[x]]))
	  }) %>% do.call(rbind, .) %>% as.data.frame()
	  colnames(columns_member) <- paste0(colnames(raw)[COL_INDEX], "_", c("GT", "AD", "DP"))
	  return(columns_member)
	}) %>% do.call(cbind, .) %>% as.data.frame()
	raw <- cbind(raw, family_columns_to_cbind)


	# Get which rows all have DP values greater than threshold for variable number of columns
  which_raw_DP_columns <-  colnames(raw) %>% grep("_DP", .)

  temp_DP_bool <- lapply(which_raw_DP_columns, function(x){
    temp1 <- raw[, x]
    temp1[temp1=="."] <- 0
    temp1 <- as.integer(temp1)
    return(temp1 >= MIN_DP_THRESHOLD)
  }) %>% do.call(cbind, .) %>% as.data.frame()
  BOOL_DP_PASS <- lapply(seq_down(temp_DP_bool), function(x){
    all(temp_DP_bool[x, ])
  }) %>% unlist

  BOOL_QUAL_PASS      <- as.numeric(raw$QUAL) > 50
  BOOL_MENDEL_PASS    <- grepl("MENDEL=True", raw$INFO)
  BOOL_INTERSECT_PASS <- grepl("set=Intersect", raw$INFO)
  BOOL_REF_PASS       <- nchar(raw$REF) > 1 & nchar(raw$REF) == nchar(raw$ALT)

  BOOL_INHERIT_FILTER <- (BOOL_DP_PASS & BOOL_QUAL_PASS & BOOL_MENDEL_PASS & BOOL_INTERSECT_PASS) | BOOL_REF_PASS

	# extract SNVs called by GATK and FreeBayes # ys_corrected_extract_rare.R
	inherited = raw[BOOL_INHERIT_FILTER, ]

	info = c(info, nrow(inherited))
	rm(raw)
	# extract rare variants
	snv_id = inherited %>% unite(SNVID, c(X.CHROM, POS, REF, ALT), sep = ":")
	rm(inherited)
	family = semi_join(snv_id, rare, by = "SNVID") %>% separate(SNVID, c("X.CHROM", "POS", "REF", "ALT"), sep = ":")
	info = c(info, nrow(family))
	rm(snv_id)
	print(info)
	if(nrow(family) > 0){
		family$CARRIER = ifelse(grepl("^1/0|^0/1", family$FA) & !family$X.CHROM %in% c("chrX", "chrY"),
		                        "fa",
		                        ifelse(family$X.CHROM %in% c("chrX", "chrY") & grepl("^1/1", family$FA),
		                                  "fa",
		                                  ifelse(family$X.CHROM != "chrY" & grepl("^1/0|^0/1", family$MO),
		                                              "mo",
		                                              NA)))
		family$FAMILY = famid
		family$PROBAND_SEX = pro_sex
		family$SIBLING_SEX = sib_sex
		family$INFO = paste(family$INFO, ";CARRIER=", family$CARRIER, ";PROBAND_SEX=", family$PROBAND_SEX, ";SIBLING_SEX=", family$SIBLING_SEX, ";FAMID=", family$FAMILY, sep = "")
		## add more columns for column number == family members, want to track: GT and AD from input file for each individual, need to check if each individual has enough read support,

		info_cols <- colnames(family)[c(1:8, grep("_(GT|AD|DP)$", colnames(family)))]
		write.table(family[which(is.na(family$CARRIER) == F), info_cols], file = paste("results/rare/families/", famid, ".txt", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F)

		info_cols <- colnames(family)[grep("_(GT|AD|DP)$", colnames(family))]
		famid = gsub("inheritance/", "", gsub(".inheritance.vcf.gz", "", file))
		cat(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t", paste0(info_cols, collapse = "\t")), sep = "\n", file = paste("results/rare/", famid, ".header", sep = ""), append = T)
		rm(family)
		gc()
	} else {
	  cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", sep = "\n", file = paste("results/rare/", famid, ".header", sep = ""), append = T)
	}
	return("done")
}

## initalize data ##
# read pedigree file
ped = read.table(ped_file, sep =" ", stringsAsFactors = F)

# calculate minor allele count
families = scan("family_ids.txt", "c")
n_parents = length(families) * 2

print(c(het_ac, hom_ac, n_parents))
# read variant counts VCF
print("reading variant counts")
counts = suppressWarnings(read_delim(counts_file, delim = "\t", comment = "#", progress = F, col_names = F))
rare = counts %>%
  separate(X8, c("POP_GT_COUNT", NA), sep = ";") %>%
  dplyr::mutate(POP_GT_COUNT = gsub("POP_GT_COUNT=", "", POP_GT_COUNT)) %>%
  separate(POP_GT_COUNT, c("hom_ref", "het", "hom_alt"), sep = ",", convert = T) %>%
  dplyr::filter(het <= het_ac & hom_alt <= hom_ac) %>%
  dplyr::transmute(SNVID = paste(X1, X2, X4, X5, sep = ":"))

rm(counts)
## load file names ##
files = dir("inheritance", ".inheritance.vcf.gz$", full.names = T)
files = dir("./", ".inheritance.vcf.gz$", full.names = T)

# load header
n = as.numeric(system(paste("zcat ", files[1], " | head -n 7000 | grep -c '#'", sep = ""), intern = T)) - 1
header = read.table(files[1], sep = "\n", comment.char = "@", nrow = n)

info = as.character(header[grep("INFO=|FILTER=", header$V1),])
new_info = gsub(">", '">', gsub("Description=", 'Description="', info))
new_info = c(new_info, '##INFO=<ID=FAMID,Number=1,Type=String,Description="Family ID">', '##INFO=<ID=CARRIER,Number=1,Type=String,Description="Carrier parent">', '##INFO=<ID=PROBAND_SEX,Number=1,Type=String,Description="Proband sex">', '##INFO=<ID=SIBLING_SEX,Number=1,Type=String,Description="Sibling sex">')

tmp = as.character(header[grep("FORMAT|INFO", header$V1, invert = T),])
new_header = c(tmp[1:2], sort(new_info), tmp[3:length(tmp)])

# output new header
print("outputting VCF header")
famid = gsub("inheritance/", "", gsub(".inheritance.vcf.gz", "", file))
write.table(new_header, file = paste("results/rare/", famid, ".header", sep = ""), sep = "\n", col.names = F, quote = F, row.names = F)


# extract rare variants from each family (runs 10 families in parallel)
print("extracting rare variants")
extract_rare(file)
