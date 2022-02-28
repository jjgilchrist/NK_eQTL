library(moloc)

#read in cis eQTL mapping aummary stats for CD226 across 5 cell types
nk.moloc <- read.table("nk.cd226.txt", header = T)
cd4.moloc <- read.table("cd4.cd226.txt", header = T)
cd8.moloc <- read.table("cd8.cd226.txt", header = T)
neut.moloc <- read.table("neut.cd226.txt", header = T)
cd14.moloc <- read.table("cd14.cd226.txt", header = T)

#create list of data frames for moloc
moloc.list <- list(nk.moloc, neut.moloc, cd14.moloc, cd4.moloc, cd8.moloc)

#run moloc
out <- moloc_test(moloc.list)
