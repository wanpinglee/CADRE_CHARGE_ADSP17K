
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
	cat("\n")
	cat("Usage: Rscript --vanilla nullModel.R <GDS> <pcair.rds> <pcrel.rds> <Pheno_file> <out_prefix>\n")
	cat("\tThe script will generate <Data_Prefix>.nullModel.rds.\n")
	q(save="no")
}


library(SeqArray)
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(vcfR)
library(GWASTools)
library(dplyr)


# Input files:
# 1. Combined GDS file: "/gds_comb.gds"
# 2. Cleaned phenotype file: "gcad_17k_phenotype_v3_cl_20220307.csv"
# 3. PCs: "pcair.RDS"
# 4. GRM: "pcair.RDS"

## Read in GDS file ##
gdsfile <- args[1]
pcairfile <- args[2]
pcrelfile <- args[3]

#2.1 Create a SeqVarData object
#Read in phenotype dataa
#File format (50 columns)
#sample.id	Sex	Age	Study DSS	Sample Set	Phenosource	AD	Sequencing Center
#Sample Name	Study	LSAC	WGS WES runID	Read length BP	Sequencer	PCR free
#Pipeline type and version	Production date Genome build	GATK gVCF GWAS concordance	Sex check	
#Freemix Chipmix Sequencing Platform 4	Platform	PCR	Baylor_HiSeq20002500	Broad_HiSeq20002500     
#Genentech_HiSeq20002500	Illumina_HiSeq20002500	WashU_HiSeq20002500	Baylor_HiSeqXTen	Broad_HiSeqXTen
#NYGC_HiSeqXTen	USUHS_HiSeqXTen	WashU_HiSeqXTen	USUHS_NovaSeq	Center_Platform length	Baylor_100	Baylor_150
#Broad_100	Broad_150	Genentech_100	Illumina_100	NYGC_150	USUHS_150	WashU_100	WashU_150	Center_Length
#pheno<-read.table("/mnt/adsp/users/leew/17k_gatk/PCA/info/adsp_self_report/aa.pheno", head=T, sep='\t')
pheno<-read.table(args[4], head=T, sep='\t')

# Output
nullmodelfile <- paste0(args[5], ".nullModel.rds")

#Load GDS
gds <- seqOpen(gdsfile)

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

#3.2 Replacement: Using existing PCs
#Import PCA output
pcs <- readRDS(pcairfile)

#Convert to data frame
pc.df <- as.data.frame(pcs$vectors)
#nrow(pc.df)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pc.df, pData(annot), by="sample.id")

#Left join PC data with GDS IDs
#Does not re-order GDS IDs
pc.df <- dplyr::left_join(sampleIds, pc.df)
nrow(pc.df)

#4.1 Null Model
#Create annotated data frame
annot <- AnnotatedDataFrame(pc.df)

#Add annotated data frame to seqData
sampleData(seqData) <- annot

# Import GRM (covariance matrix from pcrelate output)
pcrel <- readRDS(pcrelfile)
grm <- pcrelateToMatrix(pcrel, scaleKin=2)

#Subset GRM to sample IDs in phenotype file
#grm_sub<-grm[rownames(grm) %in% pheno$sample.id,colnames(grm) %in% pheno$sample.id]
#nrow(grm_sub)

#Output grm subset
#saveRDS(object = grm_sub, file = "17k_PCs_GRM_PhenoSubset_20220317.RDS")

#Merge phenotype data with PC data
idsToKeep_df<-as.data.frame(rownames(grm))
names(idsToKeep_df)<-"sample.id"
pc.df_sub <- left_join(idsToKeep_df,pc.df)
nrow(pc.df_sub)

rownames(pc.df_sub)<-pc.df_sub$sample.id

# Fit the null model
#Reference for Center-Length dummy variable: USUHS_150
#Final Null Model (2022.03.21)
#Nullmod_CenterLengthPCR <- fitNullModel(seqData, outcome="AD", covars=c("Sex", "Baylor_100","Baylor_150","Broad_100","Broad_150","Genentech_100","Illumina_100","NYGC_150","WashU_100","WashU_150","PCR.free","PC1","PC2","PC3","PC4","PC9","PC11"), cov.mat=grm_sub, verbose=FALSE, family="binomial")
#Nullmod_CenterLengthPCR <- fitNullModel(seqData, outcome="AD", covars=c("Sex", "Baylor_100","Baylor_150","Broad_100","Broad_150","Genentech_100","Illumina_100","NYGC_150","WashU_100","WashU_150","PCR.free","PC1","PC2","PC3","PC4","PC9","PC11"), cov.mat=grm, family="binomial")
#Nullmod_CenterLengthPCR <- fitNullModel(seqData, outcome="AD", covars=c("Sex", "Baylor_100","Baylor_150","Broad_100","Broad_150","Genentech_100","Illumina_100","NYGC_150","WashU_100","WashU_150","PCR.free","PC1","PC2","PC3","PC4","PC9","PC11"), cov.mat=grm)

Nullmod_CenterLengthPCR <- fitNullModel(seqData, outcome="AD", covars=c("Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "Age", "Study.DSS", "Phenosource", "PCR.free", "Baylor_100","Baylor_150","Broad_100","Broad_150","Genentech_100","WashU_150", "APOE_e2", "APOE_e4"), cov.mat=grm, family="binomial")

saveRDS(object = Nullmod_CenterLengthPCR, file = nullmodelfile)

