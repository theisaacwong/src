#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


HELP_OPTS <- c("--help", "-h", "--h", "-help", "--version", "-v")
if(any(HELP_OPTS %in% args)){
  cat("usage: Rscript pinp -1 [file1] -2 [file2]\n")
  cat("  --not: for reverse %in%\n")
  cat("  -o [output.file]: specifiy outfile instead of stdout\n")
  quit(save="no")

}

BOOL_NOT <- !("--not" %in% args)

MANDATORY_ARGS <- c("-1", "-2")
if(all(MANDATORY_ARGS %in% args) == FALSE){
  stop("-1 and/or -2 are mandatory")
  quit(save="no")
}

INPUT_1 <- args[which(args=="-1") + 1]
INPUT_2 <- args[which(args=="-2") + 1]

if(INPUT_1=="stdin" & INPUT_2=="stdin"){
  stop("only one input can be stdin")
  quit(save="no")
}


readIn <- function(STARTING_SIZE = 1000, SCALING_FACTOR = 2){
  mydata <- rep("", STARTING_SIZE)
  i <- 1
  f <- file("stdin")
  open(f)
  while(length(line <- readLines(f,n=1)) > 0) {
    mydata[i] <- line
    i <- i+1
    if(i >= length(mydata)){
      mydata <- c(mydata, rep("", length(mydata)*SCALING_FACTOR))
    }
  }
  return(mydata[mydata != ""])
}

if(INPUT_1 == "stdin"){
  DATA_1 <- readIn()
} else {
  DATA_1 <- read.table(INPUT_1, header = FALSE)$V1
}
if(INPUT_2 == "stdin"){
  DATA_2 <- readIn()
} else {
  DATA_2 <- read.table(INPUT_2, header = FALSE)$V1
}


OUTPUT <- DATA_1[(DATA_1 %in% DATA_2) == BOOL_NOT]

if("-o" %in% args){
  write.table(OUTPUT, args[which("-o"==args)+1], col.names = FALSE, row.names = FALSE, quote=FALSE)

} else {
  for(i in seq_along(OUTPUT)){
    cat(paste0(OUTPUT[i], "\n"))
  }
}


