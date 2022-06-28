#!/bin/bash

if [ "$#" -ne 3 ]
then
      echo "usage: $0 <VCF> <Sample_List> <OUT_DIR>"
      exit 1
fi

VCF=$1
SAMPLES=$2
OUT_DIR=$3
THREADS=2

SCRIPT_DIR=$(dirname "$0")
export BCFTOOLS_PLUGINS=$SCRIPT_DIR/../libs/bcftools-1.14/plugins

mkdir -p $OUT_DIR

# AN > 1690 since we like to filter variants if their GT missing > 0.05.
# bin/compact_vcf reads from stdin and generates to stdout

PARAM="VFLAGS_one_subgroup=0 && ((ABHet_one_subgroup > 0.25 && ABHet_one_subgroup < 0.75) || ABHet_one_subgroup = '.') && AN > 1690"

$SCRIPT_DIR/../bin/bcftools view --threads $THREADS --force-samples -c 1 -S $SAMPLES -Ou -f PASS $VCF \
  | $SCRIPT_DIR/../bin/bcftools annotate --threads $THREADS -x "FORMAT" \
  | $SCRIPT_DIR/../bin/bcftools view --threads $THREADS -i $PARAM \
  | $SCRIPT_DIR/../bin/bcftools +fill-tags -Oz -o $OUT_DIR/$(basename $VCF) --threads $THREADS - -- -t AN,AC,AF,hwe
