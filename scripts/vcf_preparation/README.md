---

## Files Description
### Input file
- Bi-allelic **QC** VCFs (one VCF for a autosomal) with 16,905 samples released in October 2021
- Multi-allelic VCFs (one VCF for a autosomal) with 16,906 samples released in March 2021

### Output file
- VCFs (one VCF for a autosomal) with sample and variant filtering done

---

### Execution

1. For a bi-allelic variant, keep it if (check [bi_allelic.sh](bi_allelic.sh))
	- VFLAGS_one_subgroup=0
	- (ABHet_one_subgroup > 0.25 && ABHet_one_subgroup < 0.75) || ABHet_one_subgroup = '.'
	- AN > 1690 (GT missing rate < 5%)

2. For a multi-allelic variant, (check [multi_allelic.sh](multi_allelic.sh))
	- Coordinate the sample list with bi-allelic VCFs
	- Removes the field PGT, PID, PL, and PS to reduce filesizes
	- Normalize and left-align snps and indels: bcftools norm -f $REF -m -any
	- Change low-quality GT to missing: bcftools filter -S .  -i "FMT/GQ >= 20 & FMT/DP >= 10"
	- Keep GT only
	- Update AC, AF, and AN
	- Note: ALT could be '*'

---

## Dependencies
### Tools
- HTS library v1.12 (has embedded in [libs](libs))
- bcftools v1.12 (has embedded in [libs](libs))
