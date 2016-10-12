LL <- function(Y, counts, Q){
   
   # Fits a model to observed count data. Returns parameter estimates and BIC score.
   #
   # Inputs:
   #    Y:       N by 2*D character matrix where each row is a barcode pair and 
   #             each element must equal T,C,G,A. D=7 is the barcode length.
   #    counts:  numeric vector of length N containing the read count for each barcode pair
   #    Q:       scalar, number of latent components
   #    
   #
   # Outputs:
   #    W1:     numeric vector of length N where component i indicates which latent component
   #            the first barcode i was generated from. A value of zero means the second
   #            barcode i was associated with a latent barcode instead of the first.
   #    W2:     as above but for the second barcodes
   #    X1:     Q by D character matrix of latent first barcodes. Each row is a barcode. 
   #    X2:     as above but for the second barcodes
   #    matches: total number of matches between barcodes and associated observed barcodes 
   #    mismatches: total number of mismatches
   #    beta:   Noise hyperparameter
   #    BIC:    Bayes information criterion (used for model selection)
   #
   # Example:
   # Y <-  matrix(c("T","C","C","C","C","C","C","G","T","G","A","T","C","T",
   # "T","G","G","T","G","G","T","G","A","A","A","G","C","A",
   # "T","G","T","T","G","G","T","C","T","C","C","C","G","T"), 3, 14, byrow=T) 
   # counts <- c(1,3,1)
   # Q <- 1
   # result <- LL(Y, counts, Q)
   #
   #
   # Copyright James Barrett 2015
   # Version: 1.0.0
   # Date: 9 Nov 2015
   # Contact: regmjeb@ucl.ac.uk
   
   D <- 7         # Length of barcode
   N <- nrow(Y)   # Total number of pairs
   
   # Split matrix Y into two
   Y1 <- Y[,1:7]    # First 7bp barcodes
   Y2 <- Y[,8:14]   # Second barcodes
   
   
   #---------------------------------------------------------------#
   # Initial clustering (used to define initial guesses for the latent barcodes)
   #---------------------------------------------------------------#
   
   # Generate N by N matrix of pariwise Hamming distances for Y1
   H <- hamming(Y1)
   # Cluster hierarchically (use plot(fit) to visualise)
   fit <- hclust(as.dist(H), method = 'ward.D')
   # Use the clustering to split into Q clusters
   W1.new <- cutree(fit, k=Q)
   
   H <- hamming(Y2)
   fit <- hclust(as.dist(H), method = 'ward.D')
   W2.new <- cutree(fit, k=Q)
   
   # For speed (only marginal gain)
   count.matrix1 <- vector('list',D)
   count.matrix2 <- vector('list',D)
   r <- c("T","C","G","A")
   for(b in 1:D){
      count.matrix1[[b]] <- matrix(0,N,4)
      ind <- Y1[,b]=="T"
      count.matrix1[[b]][ind,1] <- counts[ind]
      ind <- Y1[,b]=="C"
      count.matrix1[[b]][ind,2] <- counts[ind]
      ind <- Y1[,b]=="G"
      count.matrix1[[b]][ind,3] <- counts[ind]
      ind <- Y1[,b]=="A"
      count.matrix1[[b]][ind,4] <- counts[ind]
      
      count.matrix2[[b]] <- matrix(0,N,4)
      ind <- Y2[,b]=="T"
      count.matrix2[[b]][ind,1] <- counts[ind]
      ind <- Y2[,b]=="C"
      count.matrix2[[b]][ind,2] <- counts[ind]
      ind <- Y2[,b]=="G"
      count.matrix2[[b]][ind,3] <- counts[ind]
      ind <- Y2[,b]=="A"
      count.matrix2[[b]][ind,4] <- counts[ind]
   }
   
   
   #---------------------------------------------------------------#
   # MAP barcode and W solutions
   #---------------------------------------------------------------#
   
   MAX <- 10      # Maximum number of iterations to find X and W
   counter <- 0   # Count how many iterations have been done so far
   
   # This while loop will find X, then W and interatively update their values until they
   # converge to stable values.
   
   while (counter < MAX){
      
      # Allocate Q by D matrix of latent barcodes for first barcode
      X1 <- matrix(0,Q,D)
      for (b in 1:Q){
         for(d in 1:D){
            # This line computes the most common nucleotide in all barcodes that were 
            # generated from latent barcode b (as specified by the value of W1. Read 
            # counts are taken into account.
            X1[b,d] <- r[which.max(colSums(count.matrix1[[d]][W1.new==b,,drop=FALSE]))]
         }
      }
      
      X2 <- matrix(0,Q,D)
      for (b in 1:Q){
         for(d in 1:D){
            X2[b,d] <- r[which.max(colSums(count.matrix2[[d]][W2.new==b,,drop=FALSE]))]
         }
      }
      
      # Next update the values of W
      W1.old <- W1.new   # Keep track of the old values so we can test if a stable solution has been reached
      W2.old <- W2.new
      for (i in 1:N){
         
         # H1 and H2 are vectors of edit distances between barcode i and the latent barcodes
         # Equivalent to Hamming distance in this case.
         H1 <- colSums(Y1[i,]!=t(X1))
         H2 <- colSums(Y2[i,]!=t(X2))
         
         # Test whether the first or second barcode has a smaller minimum edit distance for observation i
         if(min(H1)<=min(H2)){
            # If the first barcode is closer then W1.new[i] tells us which latent barcode i comes from
            W1.new[i] <- which.min(H1)   
            # The zero in W2.new[i] means i belongs to the one of the first barcodes
            W2.new[i] <- 0
         } else {
            W2.new[i] <- which.min(H2)   
            W1.new[i] <- 0
         }
      }
      
      # Stop if the solutions aren't changing anymore
      if (identical(W1.old,W1.new) & identical(W2.old,W2.new)) break
      counter <- counter + 1
   }
   W1 <- W1.new
   W2 <- W2.new
   
   
   #---------------------------------------------------------------#
   # MAP beta and BIC score
   #---------------------------------------------------------------#
   
   # Here we count the total number of nucleotide matches and mismatches
   
   N.match <- numeric(N)
   for (i in seq(1,N)){
      
      # Only cound matches between the observed and latent barcodes that have been associated with each other
      if(min(H1)<=min(H2)){
         N.match[i] <- sum(Y1[i,]==X1[W1[i],])   
      } else {
         N.match[i] <- sum(Y2[i,]==X2[W2[i],])
      }
   }
   
   matches <- sum(N.match*counts)
   mismatches <- D*sum(counts)-matches
   beta <- mismatches/(matches+mismatches) # Noise hyperparameter
   
   # Log likelihood
   LL <- matches*log(1-beta) + mismatches*log(beta)
   # Bayes information criterion score
   BIC <- -2*LL + (2*D)*Q*log(N)
   
   # Return a list with any relevant quantities
   return(list(W1=W1,W2=W2,X1=X1,X2=X2,matches=matches,mismatches=mismatches,beta=beta,BIC=BIC))
}
