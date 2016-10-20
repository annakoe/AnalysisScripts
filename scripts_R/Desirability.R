#!/usr/bin/env Rscript

### Example Usage:
### ./Desirability.R -f /Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Enrichment_analysis/DESeq2/Candidates_using_desireability_function/p300s_DESeq_L2FC_NApvalue2zero.csv -w 4 -e 1

library("optparse")
library("desiR")
library("beeswarm")
library("ggplot2")

## Set user options
option_list = list(
  make_option(c("-f", "--inputfile"), type="character", default=NULL, 
              help="path to inputfile", metavar="character"),
  make_option(c("-w", "--weight_fold_change"), type="character", default=NULL, 
              help="weight to be assigned to Log2FC values (optional, default=2)", metavar="character"),
  make_option(c("-e", "--weight_error"), type="character", default=NULL, 
              help="weight to be assigned to SE values (optional, default=1)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
### catch user input
opt = parse_args(opt_parser)

### break if no user input
if (is.null(opt$inputfile)) 
  {
  print_help(opt_parser)
  stop("At least one argument must be supplied.\n", call.=FALSE)}

###Functions

plot_histogram_Log2FC<-function(df.column, experiment, batch){ 
  png(paste("Hist_L2FC_" ,experiment, " _ ", batch, ".png", sep=""))
  hist(df.column, breaks=30, col="midnightblue", border="white", main=paste(experiment, " - ", batch, sep=""), xlab="Log2FC")
  des.line(df.column, "d.high", des.args=c(cut1 = 0.05, cut2 = 2, des.min = 0.1))
  dev.off()
}

plot_histogram_LfcSE<-function(df.column, experiment, batch){
  png(paste("Hist_LfcSE_" ,experiment, " _ ", batch, ".png", sep=""))
  hist(df.column, breaks=30, col="grey", border="white", main=paste(experiment, " - ", batch, sep=""), xlab="LfcSE")
  des.line(df.column, "d.low", des.args=c(cut1 = 1.04, cut2 =1.3, des.min = 0.1))
  dev.off()
}

get_variables_from_desirability<-function(col.list){
  list_of_variables <- vector(mode="list", length=length(col.list))
  index=1
  for (column in col.list){
    name=paste("d.", column, sep = "")
    list_of_variables[index]<-name
    index=index+1
  }
  return(unlist(list_of_variables))
}

#read inputfile into dataframe
my_df<-read.table(opt$inputfile, sep=",", header=TRUE, na.strings="", check.names=FALSE, row.names=1)

filename=toupper(opt$inputfile)
if (grepl("SET7",filename)){
  experiment = "Set7"
} else if (grepl("P300",filename)){
  experiment = "p300"
} else {experiment = "unknown"
}

message("The experiment was for:")
print(experiment)
message("Fitting data to desirability functions and plotting histograms...")

Log2FC.cols<-colnames(my_df)[grepl("log2FoldChange", colnames(my_df))]
LfcSE.cols<-colnames(my_df)[grepl("LfcSE", colnames(my_df))]

my_df[grepl("log2FoldChange", colnames(my_df))][is.na(my_df[grepl("log2FoldChange", colnames(my_df))])]<-0
#wherever the column is log2FoldChange and the entry is NA, change that NA to 0
# my_df[grepl("log2FoldChange", colnames(my_df))] subsets columns that have log2Foldchange in the colname
# [is.na(my_df[grepl("log2FoldChange", colnames(my_df))])] checks if these columns are NA

my_df[grepl("LfcSE", colnames(my_df))][is.na(my_df[grepl("LfcSE", colnames(my_df))])]<-3
#wherever the column is LfcSE and the entry is NA, change that NA to 3

for (column in Log2FC.cols){
  this_column<-my_df[,grepl(column, colnames(my_df))] #subset all rows and column that has the right colname
  #here,cutoffs and min.desirability are hardcoded but could be made optional at user input:
  assign(paste("d.", column, sep = ""), d.high(this_column, cut1 = 0.05, cut2 = 2, des.min = 0.1)) 
  # assign(variable_name, content), allows us to make a variable name in the loop and assign content to it
  plot_histogram_Log2FC(df.column = this_column, experiment=experiment, batch=column)
}
  
for (column in LfcSE.cols){
  this_column<-my_df[,grepl(column, colnames(my_df))]
  #here,cutoffs and min.desirability are hardcoded but could be made optional at user input:
  assign(paste("d.", column, sep = ""), d.low(this_column, cut1 = 1.04, cut2 =1.3, des.min = 0.1))
  plot_histogram_LfcSE(df.column=this_column, experiment=experiment, batch=column)
}

message("Calculating overall desirability score for...")
vector.Log2FC<-get_variables_from_desirability(Log2FC.cols)
vector.LfcSE<-get_variables_from_desirability(LfcSE.cols)
vector_of_names<-c(vector.Log2FC, vector.LfcSE)

print(vector.Log2FC)
message(paste("...with a weighting of ", opt$weight_fold_change, sep=""))
message("and for:")
print(vector.LfcSE)
message(paste("...with a weighting of: ", opt$weight_error, sep=""))

d_weights=as.numeric(c(rep(opt$weight_fold_change, length(vector.Log2FC)), rep(opt$weight_error,length(vector.LfcSE))))

my_df$Desirabilty<-d.overall(sapply(vector_of_names, function(x) eval(parse(text=x))), weights=d_weights)
#sapply applies eval(parse()) to all elements of the vector
#otherwise eval() only evalutes the last element, therefore need to pass it over vector

message("plotting scores...")
png(paste("Desirability_" ,experiment, "_" , opt$weight_fold_change, "x_L2FC_", opt$weight_error, "x_LfcSE_test.png", sep=""))
plot(rev(sort(my_df$Desirabilty)), type="l", xlab="rank", ylab="Desirability")
dev.off()

ordered <- my_df[order(my_df$Desirabilty, decreasing=TRUE), ]
write.csv(ordered, paste(experiment,"_candidates_", opt$weight_fold_change, "x_L2FC_", opt$weight_error, "x_LfcSE_test.csv", sep=""))
message("DONE!")