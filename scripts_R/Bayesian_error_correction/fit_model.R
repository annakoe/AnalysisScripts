fit_model <- function(Y, counts){
   
   # Wrapper code for fitting a model to barcode count data. Returns a fitted model 
   # (see LL.R) and the corrected read count.
   #
   # Inputs:
   #    Y:       N by 2*D character matrix where each row is a barcode pair and 
   #             each element must equal T,C,G,A. D=7 is the barcode length.
   #    counts:  numeric vector of length N containing the read count for each barcode pair
   #    
   #
   # Outputs:
   #    model:  list where element q is the output from the LL function
   #    rho:    corrected read count
   #
   #
   # Copyright James Barrett 2015
   # Version: 1.0.0
   # Date: 9 Nov 2015
   # Contact: regmjeb@ucl.ac.uk
   
   
   # Total number of observed barcode pairs
   N <- length(counts)   
   # Preallocate a list for the output
   results <- list(model=vector('list',N), BIC=rep(NA,N))
   
   # ---------- Begin loop over N -------- #
   for (q in 1:N){
      
      # BREAK if there's only one barcode pair observed  
      if (N==1){
         results$rho <- 1
         break
      }
      
      # Fit model using LL function
      results$model[[q]] <- LL(Y, counts, q)
      results$BIC[q] <- results$model[[q]]$BIC
      results$rho <- which.min(results$BIC) # Current best estimate for corrected read count
      
      # BREAK if there are no mismatches and q=1. This means the model has achieved a 
      # perfect fit so we can stop.
      if ((results$model[[q]]$mismatches==0) & (q==1)){
         results$rho <- 1
         break
      }
      
      # BREAK if there are no mismatches. 
      if (results$model[[q]]$mismatches==0){
         results$rho <- which.min(results$BIC[1:(q-1)])
         break
      }
      
      # Search at least until q=10, then check to see if we've hit the minimum BIC score.
      # The BIC is non monotonic so don't search in the last five values of q.
      if ((q>=10) & (which.min(results$BIC[1:q])<(q-5))) break
   }
   # ---------- End loop over N -------- #
   
   return(results)
}