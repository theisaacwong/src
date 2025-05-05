#!/usr/bin/env Rscript
library(stringr)

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
}

mydata <- mydata[mydata != ""]
mynodes <- str_extract(mydata, "(?<=q@)e[:digit:]{3}(?=.grid)")


my_order <- c(paste0("e", c(217:226, 227:245, 246:247, 248:260, paste0("00", 1:9), paste0("0", 10:16), paste0("0", 17:22))))
my_order[grepl("e[a-z].", my_order)] <- my_order[grepl("e[a-z].", my_order)] %>% str_remove("^e")
my_age <- c(rep("2012", length(217:226)),
            rep("2015", length(227:245)),
            rep("2016", length(246:247)),
            rep("2018", length(248:260)),
            rep("2023", length(paste0("00", 1:9))),
            rep("2023", length(paste0("0", 10:16))),
            rep("2025", length(paste0("0", 17:22))))
df_block <- data.frame(order = my_order,
                       age = my_age,
                       stringsAsFactors = FALSE)



counts <- as.data.frame(sort(table(mynodes), decreasing = TRUE))
colnames(counts) <- c("String", "Freq")
counts <- counts[, c(2,1)]


temp <- lapply(1:nrow(counts), function(x){
  cat(paste0(counts$Freq[x], "\t", counts$String[x], "\n"))
})
