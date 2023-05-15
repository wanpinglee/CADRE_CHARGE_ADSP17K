#!/usr/bin/env RScript

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
	cat("\n")
	cat("Usage: Rscript --vanilla generate_gds.R <VCF_file> <Out_Prefix>\n")
	q(save="no")
}

library(SeqArray)
library(SeqVarTools)

##### Inputs #####
vcffile <- args[1]

##### Output #####
outprefix <- args[2]
gdsfile <- paste0(outprefix, ".gds")

# Load VCF
seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA", parallel=12, verbose=TRUE)
