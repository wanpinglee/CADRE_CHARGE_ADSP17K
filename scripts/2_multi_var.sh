#!/bin/bash

if [ "$#" -ne 3 ]
then
      echo "usage: $0 <VCF> <REF> <OUT_DIR>"
      exit 1
fi

VCF=$1
REF=$2
OUT_DIR=$3

mkdir -p $OUT_DIR

# AN > 1690 since we like to filter variants if their GT missing > 0.05.
# bin/compact_vcf reads from stdin and generates to stdout

bin/bcftools filter -S . -i "FMT/GQ >= 20 & FMT/DP >= 10" $VCF \
  | bin/compact_vcf \
  | bin/bcftools norm -f $REF -m- \
  | bin/bcftools view -v snps -f 'PASS' \
  | bin/bcftools +fill-tags - -- -t AC,AN,AF \
  | bin/bcftools view -i "AN > 1690" -Oz -o $OUT_DIR/$(basename $VCF)
