hamming <- function(Y){

   # Computes the pairwise Hamming distance between a matrix of genomic DNA sequences
   #
   # Inputs:
   #    Y:       N by D character matrix where each row is a barcode pair and 
   #             each element must equal T,C,G,A.
   #
   # Outpouts:
   #   H:       An N by N matrix of pairwise Hamming distances
   #
   #
   # Copyright James Barrett 2015
   # Version: 1.0.0
   # Date: 9 Nov 2015
   # Contact: regmjeb@ucl.ac.uk
   
   D <- 7
   N <- nrow(Y)
   
   YT <- Y=="T"
   HT <- YT %*% t(YT)
   
   YC <- Y=="C"
   HC <- YC %*% t(YC)
   
   YG <- Y=="G"
   HG <- YG %*% t(YG)
   
   YA <- Y=="A"
   HA <- YA %*% t(YA)
   
   H <- 1- (HT + HC + HG + HA)/D
   
   return(H)   
}
