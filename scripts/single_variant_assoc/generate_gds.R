#!/usr/bin/env RScript

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
	cat("\n")
	cat("Usage: Rscript --vanilla generate_gds.R <VCFs_Path> <Out_Prefix>\n")
	cat("\tVCF filenames should be with the prefix \"chr\" and the suffix \"vcf.gz\".\n")
	q(save="no")
}

library(SeqArray)
library(SeqVarTools)

##### Inputs #####
vcffile <- paste0(args[1], "/chr", 1:22, ".vcf.gz")

##### Output #####
outprefix <- args[2]
gdsfile <- paste0(outprefix, ".gds")

# Load VCF
seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA", parallel=12, verbose=TRUE)
