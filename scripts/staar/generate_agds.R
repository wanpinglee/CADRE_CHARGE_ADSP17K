
##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
	cat("\n")
	cat("Usage: Rscript --vanilla generate_gds.R <VCFs_Path> <Out_Path>\n")
	cat("\tVCF filenames should be with the prefix \"chr\" and the suffix \"vcf.gz\".\n")
	q(save="no")
}

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

##### Inputs #####
vcf_dir <- args[1]

##### Output #####
## Below is copied from https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/Association_Analysis_PreStep.r ##
## file directory for the output files
output_path <- args[2]
dir.geno <- output_path

## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/filter"

## annotation name. The first eight names are used to define masks in gene-centric analysis, do not change them!! 
## The others are the annotation you want to use in the STAAR procedure, and they are flexible to change.
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## channel name of the annotations. Make sure they are matched with the name, especially for the first eight one!! 
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info","/genecode_comprehensive_exonic_category","/metasvm_pred","/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

###########################################################
#           Main Function 
###########################################################

#### aGDS directory
agds_dir <- paste0(dir.geno, "/chr", seq(1,4), ".gds") 
save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

#### Annotation name catalog (alternatively, can skip this part by providing Annotation_name_catalog.csv with the same information)
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))

jobs_num <- matrix(rep(0,12), nrow=4)

for(chr in 1:4) {
	#vcffile <- paste0(vcf_dir, "/chr", chr, ".vcf.gz")
	gdsfile <- paste0(dir.geno, "/chr", chr, ".gds")
	#seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA", parallel=12, verbose=TRUE)

	genofile <- seqOpen(gdsfile)

	filter <- seqGetData(genofile, QC_label)
	position <- as.numeric(seqGetData(genofile, "position"))
	jobs_num[chr,1] <- chr
	jobs_num[chr,2] <- min(position)
	jobs_num[chr,3] <- max(position)

	seqClose(genofile)
}

# Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))
# Sliding Window
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))
# SCANG
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")
jobs_num <- as.data.frame(jobs_num)

save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))


