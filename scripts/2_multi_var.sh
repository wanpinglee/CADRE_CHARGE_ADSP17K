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

mkdir -p $OUT_DIR
SCRIPT_DIR=$(dirname "$0")

$SCRIPT_DIR/../bin/bcftools view --threads 2 -S $SAMPLES -Ou -f PASS $VCF \
  | $SCRIPT_DIR/../bin/bcftools annotate --threads 2 -Ou -x "FORMAT/PGT,FORMAT/PID,FORMAT/PL,FORMAT/PS" \
  | $SCRIPT_DIR/../bin/bcftools norm --threads 2 -m -any -f $REF -O u \
  | $SCRIPT_DIR/../bin/bcftools filter --threads 2 -S . -Ou -i "FMT/GQ >= 20 & FMT/DP >= 10" \
  | $SCRIPT_DIR/../bin/bcftools annotate --threads 2 -x "FORMAT" -Oz -o $OUT_DIR/$(basename $VCF)
