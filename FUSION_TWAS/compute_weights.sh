#!/bin/bash
#compute functional weights using cis eQTL mapping data in NK cells
#Requires Plink & R/3.4.0


#define probe ID, chromosome, BP location from list: probe_loc_fusion.txt
probe=`awk \
'NR=="'"$SGE_TASK_ID"'" {print $1}' \
probe_loc_fusion.txt`

chr=`awk \
'NR=="'"$SGE_TASK_ID"'" {print $2}' \
probe_loc_fusion.txt`

bp=`awk \
'NR=="'"$SGE_TASK_ID"'" {print $3}' \
probe_loc_fusion.txt`

#set 1Mb window
LOW="$(($bp-1000000))"
HIGH="$(($bp+1000000))"

gene="$(($SGE_TASK_ID+1))"

#extract expression data from normalised expression matrix
awk \
-v col="$gene" \
'{print 0, $1, $col}' \
expn_for_fusion.txt \
> $probe.pheno

# Get the locus genotypes for all samples and set current gene expression as the phenotype
plink2 \
--vcf /well/hill/James/nk.eqtl/data/nk_245.impute2.for_fusion.vcf dosage=DS \
--make-bed \
--pheno $probe.pheno \
--out $probe \
--chr $chr \
--from-bp $LOW \
--to-bp $HIGH \
--extract ./fusion/software/LDREF/1000G.EUR.$chr.bim \
--force-intersect

#run FUSION.compute_weights.R to compute weights, including covariate matrix of 32 expression PCs and 10 genetic PCs
Rscript \
./fusion/software/fusion_twas-master/FUSION.compute_weights.R \
--bfile $probe \
--tmp $probe.tmp \
--out $probe \
--covar covar_pc32_genopc10.txt \
--verbose 1 \
--save_hsq \
--PATH_gcta ./fusion/software/fusion_twas-master/gcta_nr_robust \
--PATH_plink ./plink \
--PATH_gemma ./fusion/software/gemma-0.98.1-linux-static \
--models blup,bslmm,lasso,top1,enet \
--rn TRUE \
--hsq_p 0.05 \
--probe $probe


# Append heritability output to hsq file
cat $probe.hsq >> total.hsq

# Clean-up
rm -f $probe.hsq $probe.*

