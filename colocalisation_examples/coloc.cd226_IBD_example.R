library(coloc)

#read in matrix of summary stats for test GWAS trait and NK cell eQTL mapping data for phenotypes (with evidence of a cis eQTL) within 250kb of a suggestive (p<1e-6) GWAS locus
total.coloc <- read.table("total.coloc.GCST003043", header = F)
colnames(total.coloc) <- c("rsid", "beta", "se", "maf", "eqtl_n", "probe", "chr", "bp", "beta_gwas", "se_gwas", "p_gwas")

#as example subset for probe_id=3940754 (CD226)
d1 <- na.omit(subset(total.coloc, probe==3940754))
#read in LD matrix for SNPs in vicinity of CD226
ld <- as.matrix(read.table("3940754.ld", header = T, row.names = 1))
#prepare eQTL list for coloc
eqtl.list <- list(beta = d1$beta, varbeta = (d1$se)^2, type = "quant", N = 245, MAF = d1$maf, snp = as.character(d1$rsid),LD=ld)
#prepare GWAS trait (as example inflammatory bowel disease) list for coloc
gwas.list <- list(pvalues = d1$p_gwas, type = "cc", N = (12882+21770), s = (12882/(12882+21770)), MAF = d1$maf, snp = as.character(d1$rsid),LD=ld)
out <- coloc.signals(eqtl.list, gwas.list,p12=1e-5, method="mask", mode="iterative") 
pp4 <- max(out$summary[,8])

