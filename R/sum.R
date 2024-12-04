#!/usr/bin/env Rscript

STARTING_SIZE <- 1000
SCALING_FACTOR <- 2
mydata <- rep("", STARTING_SIZE)
i <- 1

f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
#  write(line, stderr())
	mydata[i] <- line
	i <- i+1
	if(i >= length(mydata)){
		mydata <- c(mydata, rep("", length(mydata)*SCALING_FACTOR))
	}	


  # process line
}

mydata <- mydata[mydata != ""]

mydata <- as.numeric(mydata)
cat(paste0(sum(mydata), "\n"))


