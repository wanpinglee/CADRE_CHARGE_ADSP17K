#!/bin/bash

if [ "$#" -ne 4 ]
then
      echo "usage: $0 <VCF> <REF> <Sample_List> <OUT_DIR>"
      echo "<Sample_List> indicates all samples in the bi-alleic VCFs."
      exit 1
fi

VCF=$1
REF=$2
SAMPLES=$3
OUT_DIR=$4
THREADS=2

SCRIPT_DIR=$(dirname "$0")
export BCFTOOLS_PLUGINS=$SCRIPT_DIR/../libs/bcftools-1.14/plugins

mkdir -p $OUT_DIR

$SCRIPT_DIR/../bin/bcftools view --threads $THREADS --force-samples -c 1 -S $SAMPLES -Ou -f PASS $VCF \
  | $SCRIPT_DIR/../bin/bcftools annotate --threads $THREADS -Ou -x "FORMAT/PGT,FORMAT/PID,FORMAT/PL,FORMAT/PS" \
  | $SCRIPT_DIR/../bin/bcftools norm --threads $THREADS -m -any -f $REF -O u \
  | $SCRIPT_DIR/../bin/bcftools filter --threads $THREADS -S . -Ou -i "FMT/GQ >= 20 & FMT/DP >= 10" \
  | $SCRIPT_DIR/../bin/bcftools annotate --threads $THREADS -x "FORMAT" \
  | $SCRIPT_DIR/../bin/bcftools sort --threads $THREADS \
  | $SCRIPT_DIR/../bin/bcftools +fill-tags -Oz -o $OUT_DIR/$(basename $VCF) --threads $THREADS - -- -t AN,AC,AF,hwe
