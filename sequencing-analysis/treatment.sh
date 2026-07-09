#!/usr/bin/env bash

cp */*mutpos mutposreportsA/

cd mutposreportsA

python3 edit_mutpos_reports.py


awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} /^Chrom\t/{next} {print}' *mutposreports > compilatedclean.tab
