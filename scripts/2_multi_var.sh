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

cd /mnt/data3/old-master/leew/tools/callset_preparation/scripts

../bin/bcftools view --threads 2 -S $SAMPLES -Ou -f PASS $VCF \
  | ../bin/bcftools annotate --threads 2 -Ou -x "FORMAT/PGT,FORMAT/PID,FORMAT/PL,FORMAT/PS" \
  | ../bin/bcftools norm --threads 2 -m -any -f $REF -O u \
  | ../bin/bcftools filter --threads 2 -S . -Ou -i "FMT/GQ >= 20 & FMT/DP >= 10" \
  | ../bin/bcftools annotate --threads 2 -x "FORMAT" -Oz -o $OUT_DIR/$(basename $VCF)
