#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(googlesheets4)
options(gargle_oauth_email = TRUE)

if(length(args)==1){
  tea_type <- args[1]
} else if(length(args)==0){
  tea_type <- "costco_green"
} else {
  tea_type <- paste0(args, collapse = "_")
}

temp1 <- data.frame(year = format(Sys.Date(), "%Y") %>% as.numeric(),
                    month = format(Sys.Date(), "%m") %>% as.numeric(),
                    day = format(Sys.Date(), "%d") %>% as.numeric(),
                    time = format(Sys.time(), "%H:%M"),
                    tea = tea_type)
sheet_append("https://docs.google.com/spreadsheets/d/1zZP9tZZqLydz7tnPKlQTDWkiuuVRIwGO6m63kkV_8wc/edit#gid=0",
             data=temp1)

