#!/bin/bash

#Inputs are sorted, tabix-indexed vcf and bed files (see: https://qtltools.github.io/qtltools/)
#Randomly sample 1000 phenotypes to build null distribution for trans analysis
QTLtools \
trans \
--vcf total.vcf.gz \
--bed total.bed.gz \
--sample 1000 \
--normal \
--out nk.trans.sample

#Run approximate pass of analysis in trans to optimise PC inclusion.
#Approximate passes in trans including 0-50 PCs of expression matrix.
#SGE_TASK_ID=1-50
PC=$((SGE_TASK_ID - 1))
QTLtools \
trans \
--vcf total.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca$PC \
--adjust nk.trans.sample.best.txt.gz \
--normal \
--threshold 0.1 \
--out ./pc$PC/trans.adjust

#Determine numbers of significant eSNPs in trans with runFDR_atrans.R script (https://qtltools.github.io/qtltools/)
#SGE_TASK_ID=1-50
PC=$((SGE_TASK_ID - 1))
Rscript \
runFDR_atrans.R \
./pc$PC/trans.adjust.best.txt.gz \
./pc$PC/trans.adjust.hits.txt.gz \
0.05 \
trans.pc$PC

#nominal pass in trans inlcuding 12 PCs, outputting snp:phenotype pairs p<1e-5
QTLtools \
trans \
--vcf total.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca12 \
--nominal \
--threshold 1e-5 \
--out trans.nominal

#full permutation pass with 1000 permutations
#SGE_TASK_ID=1-1000
QTLtools \
trans \
--vcf total.recode.vcf.gz \
--bed total.bed.gz \
--cov total.QTLtools.pca12 \
--threshold 1e-5 \
--permute \
--exclude-phenotypes probe.excl \
--out trans.perm${SGE_TASK_ID} \
--seed ${SGE_TASK_ID}

#calculate FDR using for each permutation using runFDR_ftrans.R (https://qtltools.github.io/qtltools/)
Rscript \
runFDR_ftrans.R \
trans.nominal.hits.txt.gz \
trans.perm"$SGE_TASK_ID".hits.txt.gz \
perm."$SGE_TASK_ID"

#join outputs of FDR script together
sort \
-k 1b,1 \
<( awk \
'{print $1 "." $4, $10}' \
perm.1 ) \
> perm.total

for NO in {2..1000}
do
join \
perm.total \
<( sort \
-k 1b,1 \
<( awk \
'{print $1 "." $4, $10}' \
perm.$NO ) ) \
> perm.tmp
mv \
perm.tmp \
perm.total
done

#calculate the median of FDR estimates for each permutation in Rscript
d <- read.table("perm.total" , header = F, row.names = 1)
d[] <- lapply(d, function(x) {
    as.numeric(as.character(x)) })
med <- c()
med <- apply(d, 1, function(x) {
    median(x) })
out <- cbind(rownames(d), med)
write.table(out, "perm.medians", col.names = F, row.names = F, quote = F)

#join medians back to SNP, phenotype info
join \
<( sort \
-k 1b,1 \
<( awk \
'{print $1 "." $4, $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
perm.1 ) ) \
perm.medians \
> perm.total.medians

#output trans eSNP:phenotype list at FDR 0.05
awk \
'$11<0.05 {print $0}' \ 
perm.total.medians \
> perm100.fdr005.trans

