#!/bin/bash

#Inputs are sorted, tabix-indexed vcf and bed files (see: https://qtltools.github.io/qtltools/)
#Calculate PCs of expression matrix.
#Excluding phenotypes with evidence of cis eQTL at baseline
QTLtools \
pca \
--bed total.bed.gz \
--exclude-phenotypes baseline_probes_with_cis_eqtl.txt \
--scale \
--center \
--out total.QTLtools

#Generate 51 covariate files for 0-50 PCs from covariate matrix of 1st 50 PCs.
#SGE_TASK_ID=1-1000
PC=$((SGE_TASK_ID - 2))
awk -v line=$SGE_TASK_ID \
'BEGIN {print}; NR < line {print $0}' \
total.QTLtools.pca50 \
 | tail -n +2 \
> total.QTLtools.pca$PC

#Run approximate pass of analysis in cis to optimise PC inclusion.
#Permutation passes in cis including 0-50 PCs of expression matrix.
#100 permutations.
#SGE_TASK_ID=1-1000

for PC in {0..50}
do
QTLtools \
cis \
--vcf total.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca$PC \
--permute 100 \
--chunk ${SGE_TASK_ID} 1000 \
--out perm_1000.${SGE_TASK_ID}
done


#Run full permutation pass in cis (for 32 PCs).
#Permutation pass in cis. 
#10,000 permutations. 
#Grouped phenotypes (>1 probe mapping to 1 gene) with --grp-best.
#SGE_TASK_ID=1-1000

PC=32
QTLtools \
cis \
--vcf totalvcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca$PC \
--permute 10000 \
--grp-best \
--chunk ${SGE_TASK_ID} 1000 \
--out perm_10000.${SGE_TASK_ID}

#Run conditional pass in cis (for 32 PCs).
#permutations_all.thresholds.txt is output from runFDR_cis.R script (https://qtltools.github.io/qtltools/) from cis permutation run above.

PC=32
QTLtools \
cis \
--vcf total.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca$PC \
--grp-best \
--mapping permutations_all.thresholds.txt \
--chunk ${SGE_TASK_ID} 1000 \
--out cond.${SGE_TASK_ID}

#extract peak independent hits from conditional output 

cat cond.* > nk.conditional_full.txt
cat nk.conditional_full.txt \
| awk '{ if ($19 == 1) print $0}' \
> nk.conditional_top_variants.txt

#Run RTC for conditional output
#hotspots_b37_hg19.bed from: https://qtltools.github.io/qtltools/

QTLtools \
rtc \
--vcf  total.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca32 \
--hotspot hotspots_b37_hg19.bed \
--gwas-cis gwas_loci_for_rtc.txt \
nk.conditional_top_variants.txt \
--normal \
--conditional \
--grp-best \
--out rtc_results.nk_eqtl_cis_conditional.txt





