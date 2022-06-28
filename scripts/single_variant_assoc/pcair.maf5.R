
##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
	cat("\n")
	cat("Usage: Rscript --vanilla pcair.r <GDS> <Sample_Annot File> <Out_Prefix>\n")
	cat("\tVCF filenames should be with the prefix \"chr\" and the suffix \"vcf.gz\".\n")
	cat("\tSample Annot file is tab-delimited with the header: sample.id Population SuperPopulation Project\n\n")
	q(save="no")
}

library(SeqArray)
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(SNPRelate)
library(dplyr)
library(GGally)
library(RColorBrewer)

##### Inputs #####
#vcffile <- paste0(args[1], "/chr", 1:22, ".vcf.gz")
gdsfile <- args[1]

# Annot file format is
# sample.id  Population  SuperPopulation Project
# HG00096  GBR EUR 1KG
# Note: The order of samples should match to VCF.
#annotfile <- "/mnt/adsp/users/leew/17k_gatk/PCA/info/adsp_self_report/aa.sample.txt"
annotfile <- args[2]

# Ref ID file format
#HG00096
#HG00097
#HG00099
#refIDs <- read.table("/mnt/data3/old-master/leew/17k_gatk/PCA/sampleIdFinalRef.txt")

##### Output #####
outprefix <- args[3]
kingfile <- paste0(outprefix, ".king.rds")
pcairfile <- paste0(outprefix,".pcair.rds")
pcrelfile <- paste0(outprefix, ".pcrel.rds")
pccsv <- paste0(outprefix, ".pc.csv")

# Load VCF
#seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA", parallel=12, verbose=TRUE)

gds <- seqOpen(gdsfile)
annot <- read.table(annotfile, head=T, sep='\t')
annot$outcome <- rnorm(nrow(annot))
metadata <- data.frame(labelDescription=c("sample id", "population", "super population", "project", "simulated phenotype"))
annot <- AnnotatedDataFrame(annot, metadata)
#all.equal(annot$sample.id, seqGetData(gds, "sample.id"))

seqData <- SeqVarData(gds, sampleData=annot)

snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, num.thread = 12, maf = 0.05,
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)

saveRDS (object = snpset, file = paste0(outprefix, ".snpset.rds"))
saveRDS (object = pruned, file = paste0(outprefix, ".pruned.rds"))

king <- snpgdsIBDKING(gds, num.thread = 30, verbose=FALSE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)

saveRDS (object = king, file = kingfile)

pcs <- pcair(seqData, 
             kinobj=kingMat,
             divobj=kingMat,
             snp.include=pruned,
	     num.cores = 30)
#, unrel.set = refIDs)

saveRDS (object = pcs, file = pcairfile)

#PC-Relate
seqSetFilter(seqData, variant.id=pruned)
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], training.set=pcs$unrels, BPPARAM=BiocParallel::SerialParam())
seqResetFilter(seqData, verbose=FALSE)

saveRDS(object = pcrel, file = pcrelfile)

# Output PCS
#annot1 <- pData(annot) %>%
#  +         select(sample.id, SuperPopulation, Project)
pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
#pc.df <- left_join(pc.df, annot1, by="sample.id")

write.csv(pc.df, pccsv)

#ggplot(pc.df, aes(PC1, PC2, color=Population)) + geom_point() +
#  scale_color_manual(values=pop.cols)
