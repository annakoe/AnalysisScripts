# Need this library
library(reshape2)


#########################################
# Load and prepare a data file
#########################################

# Length of barcode
D <- 7
# Load up one of the data files (needs to be in the current directory)
data <- read.csv("1-Set7-C4-R1_gRNA_Q20_aligned_2mismatches1gap_uniquelymapped_readID_gRNA_5bc_3bc_length14_no_error_correction.csv",header=FALSE)
# vector of all the unique gRNA names
gRNA <- unique(data$V1)
# Total number of unique gRNAs
G <- length(gRNA)


#########################################
# Generate datasets of barcodes
#########################################

# Preallocate a list structure to hold the barcode datasets
Y <- vector('list',G)

# This loop goes through each gRNA, pulls out all the associated barcodes and puts them in a character matrix
for(mu in 1:G){
   ind <- which(data[[1]]==gRNA[mu])
   N <- length(ind)
   # Converts into character matrix (not the most elegent way...)
   Y[[mu]] <- matrix(as.vector(melt(lapply(as.character(data[[2]][ind]),strsplit,split=""))$value),nrow=N,ncol=2*D,byrow=TRUE)
}


#########################################
# Fit model for each gRNA
#########################################

# Preallocate a list of model resuls
res <- vector('list',G)

# Loop through gRNAs, for each one fit a model and get the corrected read count
# This can be parallelised for speed
for(mu in 1:G){
   # Begin tryCatch (catches any errors instead of stopping the loop)
   tryCatch({
      # Indices for barcodes matched that the current gRNA
      ind <- which(data[[1]]==gRNA[mu])
      # Vector of read counts
      counts <- data[[3]][ind]
      # Fit the model
      res[[mu]] <- fit_model(Y[[mu]], counts)
   }, error=function(e) NULL) #End tryCatch
} # End loop over gRNAs