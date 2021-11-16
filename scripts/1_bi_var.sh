#!/bin/bash

if [ "$#" -ne 2 ]
then
      echo "usage: $0 <VCF> <OUT_DIR>"
      exit 1
fi

VCF=$1
OUT_DIR=$2

mkdir -p $OUT_DIR

# AN > 1690 since we like to filter variants if their GT missing > 0.05.
# bin/compact_vcf reads from stdin and generates to stdout

bin/bcftools view $VCF \
  | bin/compact_vcf \
  | bin/bcftools view -i "VFLAGS_one_subgroup=0 && ((ABHet_one_subgroup > 0.25 && ABHet_one_subgroup < 0.75) || ABHet_one_subgroup = '.') && AN > 1690" -Oz -o $OUT_DIR/$(basename $VCF)
