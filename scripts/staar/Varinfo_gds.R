###########################################################
# The code was modified by the one composed
#       by Xihao Li, Zilin Li.
# https://github.com/xihaoli/STAARpipeline-Tutorial
###########################################################

rm(list=ls())

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
        cat("\n")
        cat("Usage: Rscript --vanilla Varinfo_gds.R <GDS> <CHR> <FAVORdatabase_chrsplit.csv> <Out_Path>\n")
	cat("\tInput:\n")
	cat("\t\t<CHR_GDS> is a GDS file for a chromosome.\n")
	cat("\t\t<CHR> is a number.\n")
        cat("\t\t<FAVORdatabase_chrsplit.csv> could be found on https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/\n")
	cat("\tOutput:\n")
	cat("\t\tA folder <Out_Path>/chr<CHR> will be created.\n")
	cat("\t\tMultiple files, VarInfo_chr<CHR>_*.csv, will be created in the folder.\n")
        q(save="no")
}

gc()

##########################################################################
#           Input
##########################################################################

### Targeted GDS
genoGDS <- args[1]

### Targeted chromosome
chr <- as.numeric(args[2])

### DB split information
file_DBsplit <- args[3] # The file could be found on https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv

### output
output_path <- args[4]

###########################################################################
#           Main Function 
###########################################################################

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
genofile <- seqOpen(genoGDS)
seqSetFilterChrom(genofile, include=chr)

## Generate VarInfo
CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)
seqResetFilter(genofile)

for(kk in 1:dim(DB_info)[1])
{
	print(kk)

	VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
	VarInfo <- data.frame(VarInfo)

	### make directory
	system(paste0("mkdir -p ", paste0(output_path, "/chr", chr)))

	write.csv(VarInfo,paste0(output_path,"/chr",chr,"/VarInfo_chr",chr,"_",kk,".csv"), quote=FALSE, row.names = FALSE)
}

seqClose(genofile)
