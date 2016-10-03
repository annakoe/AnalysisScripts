#!/usr/bin/env Rscript

#Using Bioconductor 3.2 (BiocInstaller 1.20.3), R 3.2.2
#loads BSgenome and human genome
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19) 
library(optparse)


## Set user options (defaults to finding GN19-NGG)
option_list = list(
  make_option(c("-p", "--pattern"), type="character", default="GNNNNNNNNNNNNNNNNNNNNGG", 
              help="supply pattern to search for in genome (degenerate bases are allowed)", metavar="character"),
  make_option(c("-n", "--pattern.name"), type="character", default=NULL, 
              help="supply a column identifier describing the pattern", metavar="GN20GG"),
  make_option(c("-o", "--outputfile"), type="character", default="GN20GG_masked_allregions.txt", 
              help="supply an outputfilename", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
### catch user input
opt = parse_args(opt_parser)

#define a dictionnary for which sequences to run in a DNASTringSet object,
#useful if wanting to find more than one sequence, see later sections
dict0 <- DNAStringSet(c(opt$pattern))
names(dict0)<-c(opt$pattern.name)

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
  #load a chromosome into memory
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

#call the Genome search function
GenomeSearch_masked(dict0, outfile=opt$outputfile)