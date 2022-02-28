#!/bin/bash
#requires R/3.4.0

#FUSION run for chromosome 22, performing colocalization tests
Rscript \
./fusion/software/fusion_twas-master/FUSION.assoc_test.R \
--sumstats PASS_Ulcerative_Colitis.sumstats \
--weights weight001_unique_ens.pos \
--weights_dir ./weights/ \
--ref_ld_chr ./fusion/software/LDREF/1000G.EUR. \
--chr 22 \
--coloc_P 0.001 \
--GWASN  27432 \
--PANELN panels.txt \
--out coloc.Ulcerative_Colitis.chr22.dat

#FUSION post-processing
Rscript \
./fusion/software/fusion_twas-master/FUSION.post_process.R \
--sumstats PASS_Ulcerative_Colitis.sumstats \
--input nk_wt_001.Ulcerative_Colitis.chr22.dat \
--out nk_wt_001.Ulcerative_Colitis.chr22.analysis \
--ref_ld_chr ./fusion/software/LDREF/1000G.EUR. \
--chr 22
