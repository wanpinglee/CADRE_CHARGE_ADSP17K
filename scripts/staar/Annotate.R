###########################################################
# The code was modified by the one composed
#       by Xihao Li, Zilin Li.
# https://github.com/xihaoli/STAARpipeline-Tutorial
###########################################################

rm(list=ls())

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
        cat("\n")
        cat("Usage: Rscript --vanilla Annotate.R <xsv_command> <FAVOR_DB> <CHR> <FAVORdatabase_chrsplit.csv> <Out_Path>\n")
	cat("\tNote: VarInfo_chr<CHR>_*.csv should be ready on <Out_Path>/chr<CHR>.\n")
	cat("\tInput:\n")
	cat("\t\t<xsv_command>: See the installation instruction on https://github.com/BurntSushi/xsv#installation.\n")
	cat("\t\t<FAVOR_DB>: The CSV files can be downloaded on https://favor.genohub.org/ (FAVORannotator tab).\n")
	cat("\t\t<CHR> is a number.\n")
        cat("\t\t<FAVORdatabase_chrsplit.csv> could be found on https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/\n")
	cat("\tOutput:\n")
	cat("\t\tMultiple files, Anno_chr_chr<CHR>_*.csv, will be created on <Out_Path>/chr<CHR>.\n")
        q(save="no")
}

#gc()
#memory.limit(size=10000) # Size in MB

##########################################################################
#           Input
##########################################################################

### xsv directory
xsv <- args[1]

### DB file
DB_path <- args[2]

chr <- as.numeric(args[3])

### DB split information 
file_DBsplit <- args[4] # The file could be found on https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv

### output
output_path <- args[5]

### anno channel (subset)
anno_colnum <- c(1,8:12,15,16,19,23,25:36)

###########################################################################
#           Main Function 
###########################################################################

### chromosome number
## annotate (seperate)
DB_info <- read.csv(file_DBsplit,header=TRUE)
chr_splitnum <- sum(DB_info$Chr==chr)

for(kk in 1:chr_splitnum)
{
	print(kk)
	system(paste0("mkdir -p ", paste0(output_path, "/chr", chr)))
	system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Anno_chr",chr,"_",kk,".csv"))
}

## merge info
Anno <- paste0(output_path,"chr",chr,"/Anno_chr",chr,"_",seq(1:chr_splitnum),".csv ")
merge_command <- paste0(xsv," cat rows ",Anno[1])

for(kk in 2:chr_splitnum)
{
	merge_command <- paste0(merge_command,Anno[kk])
}

merge_command <- paste0(merge_command,"> ",output_path,"chr",chr,"/Anno_chr",chr,".csv")

system(merge_command)

## subset
anno_colnum_xsv <- c()
for(kk in 1:(length(anno_colnum)-1))
{
	anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
}
anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])

system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,"chr",chr,"/Anno_chr",chr,".csv > ",output_path,"chr",chr,"/Anno_chr",chr,"_STAARpipeline.csv"))

