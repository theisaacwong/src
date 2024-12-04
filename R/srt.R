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

counts <- as.data.frame(sort(table(mydata), decreasing = TRUE))
colnames(counts) <- c("String", "Freq")
counts <- counts[, c(2,1)]


temp <- lapply(1:nrow(counts), function(x){
	cat(paste0(counts$Freq[x], "\t", counts$String[x], "\n"))
})


