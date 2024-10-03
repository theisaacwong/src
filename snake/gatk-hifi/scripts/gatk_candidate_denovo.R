library("optparse")

option_list <- list(
  make_option(c('-i', '--inputfile'), action='store', type='character', default='datafile', help='Input file which is a file that results from running the naive inheritance script'),
  make_option(c('-f', '--family'), action='store', type='character', help='Family ID'),
  make_option(c('-s', '--sampleid'), action='store', type='character', help='Sample ID')
)

opt <- parse_args(OptionParser(option_list = option_list))

inputfile = opt$inputfile

infile <- read.delim(inputfile, header=F)

type1 <- infile[which(infile$V9 %in% c("GT:AD:DP:GQ:PL", "GT:AD:DP:GQ:PGT:PID:PL:PS")),]
type1$fatherAD <- apply(type1, 1, function(x) strsplit(x[10], ":")[[1]][2])
type1$motherAD <- apply(type1, 1, function(x) strsplit(x[11], ":")[[1]][2])
type1$fatherDP <- apply(type1, 1, function(x) strsplit(x[10], ":")[[1]][3])
type1$motherDP <- apply(type1, 1, function(x) strsplit(x[11], ":")[[1]][3])
type1$probandAB <- apply(type1, 1, function(x) as.numeric(as.character(strsplit(strsplit(x[12], ":")[[1]][2], ",")[[1]][2]))/as.numeric(as.character(strsplit(x[12], ":")[[1]][3])))
type1$probandDP <- apply(type1, 1, function(x) strsplit(x[12], ":")[[1]][3])
type1$probandGQ <- apply(type1, 1, function(x) strsplit(x[12], ":")[[1]][4])

type1 <- type1[grep(',0', type1$fatherAD),]
type1 <- type1[grep(',0', type1$motherAD),]

pro <- type1[grep("denovo_kid", type1$V8),]
protype1 <- pro[which(as.numeric(as.character(pro$probandAB)) > 0.3 & as.numeric(as.character(pro$fatherDP)) > 10 & as.numeric(as.character(pro$motherDP)) > 10 & as.numeric(as.character(pro$probandGQ)) > 20 & as.numeric(as.character(pro$probandDP)) > 10),]

protype1$fatherAD <- NULL
protype1$motherAD <- NULL
protype1$fatherDP <- NULL                                                                                                                               
protype1$motherDP <- NULL 

probandGATK <- unique(protype1)
name <- paste(opt$family, opt$sampleid, sep=".")
write.table(probandGATK, file=paste(name, ".denovo.gatk.tsv", sep=""), sep="\t", quote=F, row.names=F)
