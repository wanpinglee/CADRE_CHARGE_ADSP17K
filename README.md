**Note:** This repository collects scripts for the callset preparation of the CHARGE and CADRE joint project that aims to conduct association analyses of Alzheimer’s Disease (AD) by using ADSP 17K whole genome sequence (WGS) data which is NIAGADS R3 release.

---

## Project Files Description
### Input file
- Bi-allelic **QC** VCFs (one VCF for a autosomal) with 16,905 samples released in October 2021
- Multi-allelic VCFs (one VCF for a autosomal) with 16,906 samples released in March 2021

### Output file
- VCFs (one VCF for a autosomal) with sample and variant filtering done

---

## Dependencies
### Tools
- NIAGADS/GCAD [compact_vcf](https://bitbucket.org/NIAGADS/compact_vcf/src/master/)
- HTS library v1.12 (has embedded in [libs](libs))
- bcftools v1.12 (has embedded in [libs](libs))

---

## Getting started
### Installation
```
git clone --recursive https://wanpinglee_penn@bitbucket.org/wanpinglee_penn/callset_reparation.git
cd callset_reparation
make
```

### Execution

1. Keep variants in bi-allelic VCFs if
	- VFLAGS_one_subgroup=0
	- (ABHet_one_subgroup > 0.25 && ABHet_one_subgroup < 0.75) || ABHet_one_subgroup = '.'
	- AN > 1690 (GT missing rate > 5%)

2. Multi-allelic variants
	-gentoype level QC:
 bcftools filter -S .  -i "FMT/GQ >= 20 & FMT/DP >= 10"  vcf_filename
	-variant level QC:

## License
The implementation is available for academic and nonprofit use for free [LICENSE.md](LICENSE.md).

---

## Contacts
Gina Marie Peloso - [Email](gpeloso@bu.edu)

Wan-Ping Lee - [Email](wan-ping.lee@pennmedicine.upenn.edu)

---
## Acknowledgments


