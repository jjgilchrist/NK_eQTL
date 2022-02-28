#!/bin/bash

for ID in {1..5}
do

#look up GWAS catalog ID and trait from list of 5
gwas_cat=`awk \
'NR=="'"$ID"'" {print $1}' \
gwas.cat`
trait=`awk \
'NR=="'"$ID"'" {print $2}' \
gwas.cat`

# process summary statistics file for FUSION. munge_sumstats.py available from: https://github.com/bulik/ldsc
./munge_sumstats.py \
--sumstats "$gwas_cat".tsv.gz \
--out $trait

done
