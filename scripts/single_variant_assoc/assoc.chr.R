

args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)!=5) {
	cat("\n")
	cat("Usage: Rscript --vanilla assoc.R <gds> <pheno> <null.model> <out_prefix> <chr>\n")
	cat("\tThe script will search <Data_Prefix>.gds and <Data_Prefix>.nullModel.rds as inputs.\n")
	cat("\tThe script will generate <Data_Prefix>.assoc.rds and <Data_Prefix>.assoc.csv.\n")
	q(save="no")
}

library(SeqArray)
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(vcfR)
library(GWASTools)
library(dplyr)


## Read in GDS file ##
gdsfile <- args[1]
phenofile <- args[2]
nullmodelfile <- args[3]

# Output
out_prefix <- args[4]
chr <- args[5]
assocfile <- paste0(out_prefix, ".chr", chr, ".assoc.rds")
assoccsv <- paste0(out_prefix, ".chr", chr, ".assoc.csv")

print(gdsfile)
gds <- seqOpen(gdsfile)

#Read in phenotype dataa
#File format (50 columns)
#sample.id	Sex	Age	Study DSS	Sample Set	Phenosource	AD	Sequencing Center
#Sample Name	Study	LSAC	WGS WES runID	Read length BP	Sequencer	PCR free
#Pipeline type and version	Production date Genome build	GATK gVCF GWAS concordance	Sex check	
#Freemix Chipmix Sequencing Platform 4	Platform	PCR	Baylor_HiSeq20002500	Broad_HiSeq20002500     
#Genentech_HiSeq20002500	Illumina_HiSeq20002500	WashU_HiSeq20002500	Baylor_HiSeqXTen	Broad_HiSeqXTen
#NYGC_HiSeqXTen	USUHS_HiSeqXTen	WashU_HiSeqXTen	USUHS_NovaSeq	Center_Platform length	Baylor_100	Baylor_150
#Broad_100	Broad_150	Genentech_100	Illumina_100	NYGC_150	USUHS_150	WashU_100	WashU_150	Center_Length
pheno<-read.table(phenofile, head=T, sep='\t')


#Create list of IDs from GDS file
sampleIds <- data.frame( sample.id=seqGetData(gds, 'sample.id') )
nrow(sampleIds)

#Left join phenotype data with GDS data
#Matches IDs where possible (some GDS IDs not in phenotype file)
#Does not re-order GDS IDs
pheno5 <- dplyr::left_join(sampleIds, pheno)
nrow(pheno5)


#Create metadata for AnnotatedDataFrame
#names(pheno5)
metadata <- data.frame(labelDescription=c("Sample ID","Sex","Age","Study DSS","Sample Set","Phenosource","AD","Sequencing Center","Sample Name","Study","LSAC","WGS WES","runID","Read length BP","Sequencer","PCR free","Pipeline type and version","Production date","Genome build","GATK gVCF GWAS concordance","Sex check","Freemix","Chipmix","Sequencing Platform 4","Platform","PCR","Baylor_HiSeq20002500","Broad_HiSeq20002500","Genentech_HiSeq20002500","Illumina_HiSeq20002500","WashU_HiSeq20002500","Baylor_HiSeqXTen","Broad_HiSeqXTen","NYGC_HiSeqXTen","USUHS_HiSeqXTen","WashU_HiSeqXTen","USUHS_NovaSeq","Center_Platform","length","Baylor_100","Baylor_150","Broad_100","Broad_150","Genentech_100","Illumina_100","NYGC_150","USUHS_150","WashU_100","WashU_150","Center_Length", "APOE_e2", "APOE_e4"), row.names=names(pheno5))


#Create annotated data frame
annot <- AnnotatedDataFrame(pheno5, metadata)

#Merge GDS and annotated data frame to create SeqVarData
seqData <- SeqVarData(gds, sampleData=annot)
seqSetFilterChrom(seqData, include=chr)


#Left join PC data with GDS IDs
#Does not re-order GDS IDs

Nullmod_GRM<-readRDS(nullmodelfile)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, Nullmod_GRM, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)

saveRDS(object = assoc, file = assocfile)
write.csv(assoc, file = assoccsv)
