#!/usr/bin/env Rscript

### Example Usage:
### ./Plot_candidates.R -f /Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Enrichment_analysis/DESeq2/Candidates_using_desireability_function/p300_candidates_4xL2FC_1xLfcSE.csv -c ~/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Enrichment_analysis/Naive_Approach_gRNA_enrichment/Bayesian_counts/missing_as_0/All_p300_counts_ttest.csv
library("optparse")
library("desiR")
library("beeswarm")
library("ggplot2")

### Set user options
option_list = list(
  make_option(c("-f", "--inputfile"), type="character", default=NULL, 
              help="path to inputfile", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="path to counts", metavar="character"),
  make_option(c("-o", "--outputfilename"), type="character", default=NULL, 
              help="path to counts", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
### catch user input
opt = parse_args(opt_parser)


###Functions:

check_user_input<-function(user_arg){
  if (is.null(user_arg)) {
    print_help(opt_parser)
    stop("At least three arguments must be supplied.\n", call.=FALSE)}
}

split_samplesheet<-function(samplesheet){
  #takes a vector of samplenames and splits them into samples and controls
  samplenames=toupper(samplesheet)
  samples<-samplesheet[grepl("P300",samplenames) | grepl("SET7-",samplenames) | grepl("SV",samplenames) | grepl("PV",samplenames)]
  controls<-samplesheet[grepl("EMT",samplenames) | grepl("AV",samplenames)]
  return.list<-list("controls" = controls, "samples" = samples)
  return(return.list)
}

plot_counts<-function(dataframe, controls, samples, gRNA, rank, experiment){
  plot.controls<-unlist(dataframe[gRNA, controls, drop=TRUE])
  plot.samples<-unlist(dataframe[gRNA, samples, drop=TRUE])
  
  png(paste(opt$outputfilename,"_Rank_", rank, "_", gRNA, "_", experiment, ".png", sep=""))
  data.plt<-list(plot.controls, plot.samples)
  boxplot(data.plt, main=paste('Candidate gRNA - ', gRNA), cex.main=1.4, ylab="gRNA count (UMI)", cex.lab=1.2, las=1, cex.axis=1.2, col = c(rgb(120, 120, 120, maxColorValue=255), rgb(224, 108, 16, maxColorValue=255)), outcol=NA, names=c('controls', 'samples'), outlier.colour = NA)
  beeswarm(data.plt, pch=16, cex=1.2, add=TRUE)
  dev.off()
}

plot_counts_for_each_batch<-function(batch.names){
  for (batch in batch.names){
    split.batch <- split_samplesheet(get(batch)) #get() takes string and looks for variable with that string as name
    samples = split.batch$samples
    controls = split.batch$controls
    
    message(paste("For experiment ", batch, " found the following samples (", length(samples), " samples in total): ", sep=""))
    print(samples)
    message("...and the following controls (", length(controls), " in total):", sep="")
    print(controls)
    
    message("Plotting counts...")
    rank=1
    for (gRNA in top20) {
      plot_counts(dataframe=df.counts, controls=controls, samples=samples, gRNA=gRNA, rank=rank, experiment=batch)
      rank = rank+1
    }
  }
}
  


###SCRIPT
# check if user input, break if no user input
check_user_input(opt$inputfile)
check_user_input(opt$counts)
check_user_input(opt$outputfilename)

#read inputfiles into dataframes
df.D<-read.table(opt$inputfile, sep=",", header=TRUE, na.strings="", check.names=FALSE, row.names=1)
df.counts<-read.table(opt$counts, sep=',', header=TRUE, check.names=FALSE, row.names=1)

all<-split_samplesheet(colnames(df.counts))

top20<-rownames(df.D)[1:20]
rank=1

for (gRNA in top20) {
  plot_counts(dataframe=df.counts, controls=all$controls, samples=all$samples, gRNA=gRNA, rank=rank, experiment='all')
  rank = rank+1
}

###plot counts for individual experiments

message("Looking for samplesheets...")

filename=toupper(opt$inputfile)

if (grepl("SET7",filename)) {
  
  message("Loading samplesheets for Set7 experiments")
  Set7.Batch3<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch3.txt")[,1])
  Set7.Batch5<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch5.txt")[,1])
  Set7.Batch0209<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS0209_Set7.txt")[,1])
  Set7.Batch2608<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS2608_Set7.txt")[,1])
  
  batch.names<-c("Set7.Batch3", "Set7.Batch5", "Set7.Batch0209", "Set7.Batch2608")
  plot_counts_for_each_batch(batch.names)
  message("DONE!")
  
} else if (grepl("P300",filename)) {
  
  message("Loading samplesheets for p300 experiments")
  p300.Batch4<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch4.txt")[,1])
  p300.Batch8<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch8.txt")[,1])
  p300.Batch0209<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS0209_p300.txt")[,1])
  p300.Batch2608<-as.vector(read.table("/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS2608_p300.txt")[,1])

  batch.names<-c("p300.Batch4", "p300.Batch8", "p300.Batch0209", "p300.Batch2608")
  plot_counts_for_each_batch(batch.names)
  message("DONE!")

} else stop("Neither p300 nor Set7. cannot find samplefiles.\n", call.=FALSE)

