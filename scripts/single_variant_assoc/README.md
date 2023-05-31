---

## Files Description
### Input file
- Per chromosome VCF

### Output file
- gds, pcari.rds, pcrel.gds, king.rds, pc.csv, .nullModel.rds, assoc.rds and assoc.csv

---

### Execution

1. Generate GDS  [scripts/single_variant_assoc/generate_gds.R](generate_gds.R)
	- VCFs should be named by chr#.vcf.gz
2. Perform a Principal Components Analysis [scripts/single_variant_assoc/pcair.maf5.R](pcair.maf5.R)
	- PC-AiR accounts for sample relatedness (known or cryptic) to provide accurate ancestry inference.
	- Sample Annot file is tab-delimited with a header, consisting of sample.id, Population, SuperPopulation, and Project.
	- The analysis will only use sample.id in the annot file.
3. Generate null model [scripts/single_variant_assoc/nullModel.maf5.R](nullModel.maf5.R)
	- fitNullModel fits a model with random effects specified by covariance structures.
	- fitNullModel allows for the inclusion of a polygenic random effect using a kinship matrix or genetic relationship matrix (GRM).
	- Pheno file is tab-delimited consisting of covariant with a header.
4. Perform association tests [scripts/single_variant_assoc/assoc.chr.R](assoc.chr.R)
	- assocTestSingle performs genotype association tests with the outcome using the null model obtained from Step 3.

---

## Dependencies
### Tools
- Biobase
- dplyr
- GENESIS
- GGally
- SeqArray
- SeqVarTools
