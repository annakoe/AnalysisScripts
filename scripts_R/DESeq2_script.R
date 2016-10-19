#!/usr/bin/env Rscript

### Example Usage:
### ./DESeq2_script.R -c path_to/DeSeq2/Inputfiles/Batch8_p300_counts.csv -e path_to/additional_files/ColData_forDESeq2/Batch8_DeSeq2_ColData

library("optparse")


## Set user options
option_list = list(
  make_option(c("-c", "--counttable"), type="character", default=NULL, 
              help="path to table of counts", metavar="character"),
  make_option(c("-e", "--experiment"), type="character", default=NULL, 
              help="path to experiment info (col data)", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default='DESeq2_results.csv', 
              help="outputfilename", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
### catch user input
opt = parse_args(opt_parser)

### break if no user input
if (is.null(opt$counttable) | is.null(opt$experiment) ){
  print_help(opt_parser)
  stop("At least two arguments must be supplied.\n", call.=FALSE)}
  
library ("DESeq2")
#read inputfile into dataframe
counts<-read.table(opt$counttable, sep=",", header=TRUE, na.strings="NaN", check.names=FALSE, row.names=1)
experiment <- read.table(opt$experiment, sep=",", header=TRUE, blank.lines.skip = TRUE, row.names = "Sample.name")
 
counts[is.na(counts)] <- 0  #replace NAs with 0

dds<- DESeqDataSetFromMatrix(countData=counts, colData = experiment, design = ~treatment)

#converting counts to integer mode
dds<-DESeq(dds)
res<-results(dds)

write.csv(res, opt$outputfile)

