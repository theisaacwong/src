library("stringr")
library("ggplot2")
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tidyverse"))
library("fortunes")
library("cowsay")
library("parallel")
library("googlesheets4")

# COBALT
# On average, any data manipulation that can be described in a sentence or two of English can be
# programmed in one line in R. If you find yourself writing a long 'for' loop to do something that
# sounds simple, take a step back and research if an existing combination of functions can easily
# handle your request.
# -- Erik Iverson
# R-help (June 2010)

srt <- function(x){sort(table(x), decreasing = TRUE)}

u <- function(x) {unique(x)}
lu <- function(x) {length(unique(x))}

heda <- function(x) {head(x)}
heaD <- function(x) {head(x)}
h <- function(x) {head(x)}

read <- function(x){fread(x) %>% as.data.frame()}
write <- function(data, file) {write.table(x = data, file = file, quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")}

seq_down <- function(x){ 1:nrow(x)}
lo <- function(x){sapply(x, length)==1}
ld <- function(x){sapply(x, length) %>% table}



wrapper <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}


# initialize with
# myMap <- HashMap$new()
HashMap <- setRefClass("HashMap",
                       fields = list(hash = "environment"),
                       methods = list(
                         initialize = function(...) {
                           hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 1000L)
                           callSuper(...)
                         },
                         # get = function(key) { unname(base::get(key, hash)) },
                         get = function(key) { base::get(key, hash) },
                         put = function(key, value) { base::assign(key, value, hash) },
                         append = function(key, appendant){
                           if(base::exists(key, hash)){
                             base::assign(key, base::append(unname(base::get(key, hash)) ,appendant), hash)
                           } else {
                             base::assign(key, appendant, hash)
                           }
                         },
                         containsKey = function(key) { base::exists(key, hash) },
                         keySet = function() { base::ls(hash) },
                         toString = function() {
                           keys <- ls(hash)
                           rval <- vector("list", length(keys))
                           names(rval) <- keys
                           for(k in keys){
                             # rval[[k]] <- unname(base::get(k, hash))
                             rval[[k]] <- base::get(k, hash)
                           }
                           return(rval)
                         },
                         populate = function(dataframe, keyCol, valCol){
                           if(is.character(keyCol) | is.character(valCol)){
                             keyCol = colnames(dataframe)[which(colnames(dataframe) == keyCol) ]
                             valCol = colnames(dataframe)[which(colnames(dataframe) == valCol) ]
                           }
                           for(i in 1:nrow(dataframe)){
                             base::assign(dataframe[i, keyCol], dataframe[i, valCol], hash)
                           }
                         },
                         populate_append = function(dataframe, keyCol, valCol){
                           if(is.character(keyCol) | is.character(valCol)){
                             keyCol = colnames(dataframe)[which(colnames(dataframe) == keyCol) ]
                             valCol = colnames(dataframe)[which(colnames(dataframe) == valCol) ]
                           }
                           for(i in 1:nrow(dataframe)){
                             key <- dataframe[i, keyCol]
                             appendant <- dataframe[i, valCol]
                             if(base::exists(key, hash)){
                               base::assign(key, base::append(unname(base::get(key, hash)) ,appendant), hash)
                             } else {
                               base::assign(key, appendant, hash)
                             }
                           }
                         },
                         head = function(){
                           keys <- utils::head(ls(hash))
                           rval <- vector("list", length(keys))
                           names(rval) <- keys
                           for(k in keys){
                             rval[[k]] <- unname(base::get(k, hash))
                           }
                           return(rval)
                         },
                         random = function(n=10){
                           keys <- sample(ls(hash), n)
                           rval <- vector("list", length(keys))
                           names(rval) <- keys
                           for(k in keys){
                             rval[[k]] <- unname(base::get(k, hash))
                           }
                           return(rval)
                         },
                         size = function() { length(hash) },
                         isEmpty = function() { length(hash)==0 },
                         remove = function(key) {
                           .self$hash[[key]] <- NULL
                           rm(key, envir = .self$hash)},
                         clear = function() { hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 1000L) }
                       )
)


# get
"%G%" <- function(hashMap, key) {hashMap$get(key)}


bname <- function(file_path){
  lapply(file_path, function(x){
    temp1 <- str_split(x, "/") %>% unlist
    temp1[length(temp1)]  %>% return
  }) %>% unlist %>% return
}

# default(write.table) <- list(row.names <- FALSE, col.names <- TRUE, sep <- "\t", quote <- FALSE)
# formals(utils::write.table, envir = environment(utils::write.table))$row.names <- FALSE
# formals(utils:::write.table)$col.names <- TRUE
# formals(utils:::write.table)$sep <- "\t"
# formals(utils:::write.table)$quote <- FALSE

if(interactive()){
  fortune() %>% paste(collapse = "\n\t-") %>% say(by="egret")
}





