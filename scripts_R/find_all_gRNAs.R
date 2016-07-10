#!/usr/bin/env Rscript

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19) 
#loads BSgenome and human genome

#define a dictionnary for which sequences to run in a DNASTringSet object,
#useful if wanting to find more than one sequence, see later sections
dict0 <- DNAStringSet(c("GNNNNNNNNNNNNNNNNNNNNGG"))
names(dict0)<-c("GN20GG")

#define the function 'writeHits' that writes output
writeHits <- function(seqname, matches, strand, file="", append=FALSE) {
  if (file.exists(file) && !append)
    warning("existing file ", file, " will be overwritten with 'append=FALSE'")
  
  if (!file.exists(file) && append)
    warning("existing file ", file, " will be overwritten with 'append=FALSE'")
  
  #define output format as chr, start, stop, DNA strand and name 
  hits <- data.frame(seqname=rep.int(seqname, length(matches)), start=start(matches), end=end(matches), strand=rep.int(strand, length(matches)), patternID=names(matches), check.names=FALSE)
  
  #write output to file
  write.table(hits, file=file, append=append, quote=FALSE, sep="\t", row.names=FALSE, col.names=!append)}


#define the function GenomeSearch_masked, that searches for all occurrences of a degenerate sequence in the repeat masked human genome
GenomeSearch_masked <- function(dict0, outfile=""){
  library(BSgenome.Hsapiens.UCSC.hg19) 
  seqnames<-seqnames(Hsapiens) 
  seqnames_in1string <- paste(seqnames, collapse=", ") 
  cat("Target:", providerVersion(Hsapiens), "chromosomes", seqnames_in1string, "\n") 
  append <- FALSE

#loop through chromosomes
for (seqname in seqnames) 
{
  subject <- Hsapiens[[seqname]]
  
  #repeat-mask the chromosome
  active(masks(subject))["RM"] <- TRUE
  
  cat(">>> Finding all hits in chromosome", seqname, "...\n")
  
  #within the loop through the chomosomes, loop through the sequences 
  #defined in dict0 and find matches
  for (i in seq_len(length(dict0))) {
    patternID <- names(dict0)[i]
    pattern <- dict0[[i]]
    plus_matches <- matchPattern(pattern, subject, fixed=c(pattern=FALSE, subject=TRUE))
    names(plus_matches) <- rep.int(patternID, length(plus_matches))
    #invoke the writeHits function defined above
    writeHits(seqname, plus_matches, "+", file=outfile, append=append)
    append <- TRUE
    
    rcpattern <- reverseComplement(pattern)
    minus_matches <- matchPattern(rcpattern, subject, fixed=c(pattern=FALSE, subject=TRUE))
    names(minus_matches) <- rep.int(patternID, length(minus_matches))
    
    #invoke the writeHits function defined above
    writeHits(seqname, minus_matches, "-", file=outfile, append=append)
  }
  cat(">>> DONE\n") 
}
}

GenomeSearch_masked(dict0, outfile="GN20GG_masked_Hsapiens_hg19.txt")

