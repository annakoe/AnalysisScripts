#!/usr/bin/env Rscript

### Example Usage:
### ./Print_counts_missing_as_0.R -f ../DeSeq2/Inputfiles/BatchSS2608_Set7_counts.csv

library("optparse")
suppressMessages(library("genefilter"))

##Functions

find_last_element <- function(my_string, character_to_split_at) {
  split.string=unlist(strsplit(as.character(my_string), character_to_split_at))
  last_string_element=split.string[length(split.string)]
  return(last_string_element)
}

find_nth_element <- function(my_string, character_to_split_at, position) {
  split.string=unlist(strsplit(as.character(my_string), character_to_split_at))
  nth_string_element=split.string[1]
  return(nth_string_element)
}

## Set user options
option_list = list(
  make_option(c("-f", "--inputfile"), type="character", default=NULL, 
              help="path to inputfile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
### catch user input
opt = parse_args(opt_parser)

### break if no user input
if (is.null(opt$inputfile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied.\n", call.=FALSE)
}

#read the inputfile into a dataframe
dataframe<-read.table(opt$inputfile, sep=",", header=TRUE, na.strings="", check.names=FALSE, row.names=1)
#replace NAs with 0
dataframe[is.na(dataframe)] <- 0
#keep a copy of the original df
original_df<-dataframe

#extract the column names from the dataframe and make lowercase
samplenames=toupper(colnames(dataframe))

#if column name contains P300, SET7, SV, or PV change the column name to "sample"
colnames(dataframe)[grepl("P300",samplenames) | grepl("SET7",samplenames) | grepl("SV",samplenames) | grepl("PV",samplenames)] <- "sample"
#if column name contains AV or EMT change the column name to "control"
colnames(dataframe)[grepl("EMT",samplenames) | grepl("AV",samplenames)] <- "control"

groups<-colnames(dataframe)
results_by_row <- rowttests(as.matrix(dataframe), factor(groups), tstatOnly = FALSE)

original_df["p.value.ttest.unadj"]<-results_by_row$p.value
original_df$sample.mean<-apply(dataframe[colnames(dataframe)=="sample"],1 ,mean)
original_df$control.mean<-apply(dataframe[colnames(dataframe)=="control"],1 ,mean)
original_df$FC<-original_df$sample.mean/original_df$control.mean
original_df$L2FC<-log2(original_df$sample.mean/original_df$control.mean)

#split the inputfilename and take last element
inputfilename=find_last_element(opt$inputfile, "/")
outputfilename=find_nth_element(inputfilename, "[.]")
outputfilename=paste0(outputfilename, "_missing_as_0.csv")

write.csv(original_df, outputfilename)