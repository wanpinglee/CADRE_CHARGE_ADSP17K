**Note:** This repository collects scripts for CADRE and CHARGE joint project that aims to conduct association analyses of Alzheimer’s Disease (AD) by using ADSP 17K whole genome sequence (WGS) data which is NIAGADS R3 release.

---

## Project Files Description
### Input file
- Bi-allelic **QC** VCFs (one VCF for a autosomal) with 16,905 samples released in October 2021
- Multi-allelic VCFs (one VCF for a autosomal) with 16,906 samples released in March 2021
- Phenotype file

### Output file
- VCFs (one VCF for a autosomal) with sample and variant filtering done
- Summary csv files of association analyses

---

## Steps
1. VCF preparation [README](scripts/vcf_preparation/README.md)
2. Single-variant association analysis for variant with MAF > 0.5% [README](scripts/single_variant_assoc/README.md)
3. Gene-based testing for coding variant with MAF < 1% using STAAR [README](scripts/staar/README.md)
4. Rare noncoding variants


## Dependencies
### Tools
- HTS library v1.12 (has embedded in [libs](libs))
- bcftools v1.12 (has embedded in [libs](libs))

---


## License
The implementation is available for academic and nonprofit use for free [LICENSE.md](LICENSE.md).

---

## Contacts

Wan-Ping Lee - [Email](wan-ping.lee@pennmedicine.upenn.edu), Research Assistant Professor | Associate Director for IT of NIAGADS, Department of Pathology and Laboratory Medicine, University of Pennsylvania

Gina Marie Peloso - [Email](gpeloso@bu.edu), Associate Professor, Department of Biostatistics, Boston University

---
## Acknowledgments


