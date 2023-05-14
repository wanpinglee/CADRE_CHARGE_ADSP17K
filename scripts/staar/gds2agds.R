###########################################################
# The code was modified by the one composed
#       by Xihao Li, Zilin Li.
# https://github.com/xihaoli/STAARpipeline-Tutorial
###########################################################

rm(list=ls())

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
        cat("\n")
        cat("Usage: Rscript --vanilla gds2agds.R <CHR_GDS> <CHR> <FAVOR_Anno_folder>\n")
	cat("\tNote: <CHR_GDS> will be modified!!!\n")
	cat("\tInput:\n")
	cat("\t\t<CHR_GDS> is a GDS file for a chromosome..\n")
	cat("\t\t<CHR> is a number.\n")
	cat("\t\tThe folder contains Anno_chr_chr<CHR>_*.csv.\n")
        q(save="no")
}

gc()

##########################################################################
#           Input
##########################################################################

### gds file
gds <- args[1]
chr <- as.numeric(args[2])
### annotation file (output of Annotate.R)
FAVOR_Anno_file <- args[3]

###########################################################################
#           Main Function 
###########################################################################

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)

### read annotation data
FunctionalAnnotation <- read_csv(FAVOR_Anno_file,
col_types=list(col_character(),col_double(),col_double(),col_double(),col_double(),
col_double(),col_double(),col_double(),col_double(),col_double(),
col_character(),col_character(),col_character(),col_double(),col_character(),
col_character(),col_character(),col_character(),col_character(),col_double(),
col_double(),col_character()))

dim(FunctionalAnnotation)

## rename colnames
colnames(FunctionalAnnotation)[2] <- "apc_conservation"
colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"
colnames(FunctionalAnnotation)[9] <- "apc_protein_function"

## open GDS
genofile <- seqOpen(gds, readonly = FALSE)
seqSetFilterChrom(genofile, include=chr)

Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation, compress="LZMA_RA")

seqClose(genofile)
