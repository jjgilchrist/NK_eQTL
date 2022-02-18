library(qvalue)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(powerEQTL)
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

###Figure S1: The effects of principal component incorporation for cis and trans eQTL mapping.###
#Read in cis data
cis.data <- read.table("FigureS1.cis_data.txt", header = T)

#Figure S1 cis.
cols <- brewer.pal(8,"Set2")
title <- expression(paste(italic("cis"), " eQTL"))
cis_plot <- ggplot(data=cis.data, aes(x=pcs, y=sig.genes)) +
  geom_line(color=cols[3])+
  geom_point(color=cols[3])+
  labs(x="PCs", y = "#eQTLs - FDR<0.05")+
  ylim(0,5000)+
  ggtitle(title) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title=element_text(size=25, face = "bold"))

#Read in trans data
trans.data <- read.table("FigureS1.trans_data.txt", header = T)

#Figure S1 trans.
title <- expression(paste(italic("trans"), " eQTL"))
trans_plot <- ggplot(data=trans.data, aes(x=pcs, y=sig.out)) +
  geom_line(color=cols[3])+
  geom_point(color=cols[3])+
  labs(x="PCs", y = "#eQTLs - FDR<0.05")+
  ylim(0,100)+
  ggtitle(title) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title=element_text(size=25, face = "bold"))

#Construct Figure S1
FigureS1 <- ((cis_plot)/(trans_plot))

###Figure S2: Gene expression in NK cells and sharing of eQTL between NK cells and other primary immune cells.###
#read in median expression data and whether each eQTL is shared/unique
expn.uniq <- read.table("FigureS2.data.txt", header = T)

#Figure S2
cols <- brewer.pal(8,"Set2")
p_labels = data.frame(expt = c("base"), label = c("italic(P)==0.152"))
p1 = ggplot(expn.uniq, aes(x=factor(unique), y=exp))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.05, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1), labels = c("Shared", "Unique")) + aes(fill = factor(unique), col=factor(unique)) + scale_fill_manual(values = cols[c(2,1)]) + scale_colour_manual(values = cols[c(2,1)])
p5 = p4 + ylab("Median expression")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text.y=element_text(size=12), axis.text.x=element_text(size=15),
                             axis.title=element_text(size=15), axis.title.x=element_blank())
FigureS2 = p6 + geom_text(x=1.5, y=13.5, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

###Figure S3: Functional annotation of NK cell eQTL.###
#Read in data for frequency plots representing the density of eQTL loci surrounding functional genomic tracks in primary NK cell eQTL.
total.density <- read.table("FigureS3A.data.txt", header = F)
colnames(total.density) <- c("start", "stop", "count", "feature")
total.density$pos <- (total.density$start+total.density$stop)/2
state.labs <- c("CTCF", "Enhancer", "Promoter flank", "Repressed", "Transcribed", "Promoter/TSS", "Weak Enhancer")
names(state.labs) <- c("CTCF", "E", "PF", "R", "T", "TSS", "WE")

#Figure S3A
p_labels = data.frame(feature = c("CTCF", "E", "PF", "R", "T", "TSS", "WE"), 
                      label = c("NS", "italic(P)==0.0015", "NS", "italic(P)<2.00%*%10^-8", "italic(P)<2.00%*%10^-8", "italic(P)<2.00%*%10^-8", "italic(P)==0.0087"))
cols <- brewer.pal(8,"Set2")
FigureS3A<- ggplot(total.density, aes(x=pos/1000000, y=count)) +
  geom_line() +
  theme_bw() +
  xlab("Distance to QTLs (Mb)") +
  ylab("#annotations/kb") +
  ggtitle("Primary") +
  facet_grid(. ~ feature, labeller = labeller(feature = state.labs)) +
  theme(panel.spacing.x = unit(4, "mm")) +
  aes(fill = as.factor(feature), col=as.factor(feature)) + scale_fill_manual(values = cols[c(8,4,8,3,4,4,4)]) + scale_colour_manual(values = cols[c(8,4,8,3,4,4,4)]) +
  theme(
    strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
      size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
    axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold")) +
  geom_text(x=0, y=(max(total.density$count)+100), aes(label=label), data=p_labels, parse=T, inherit.aes=F, size = 4) +
  ylim(NA, (max(total.density$count)+100))

#Read in data for frequency plots representing the density of eQTL loci surrounding functional genomic tracks in conditional NK cell eQTL.
conditional.density <- read.table("FigureS3B.data.txt", header = F)
colnames(conditional.density) <- c("start", "stop", "count", "feature")
conditional.density$pos <- (conditional.density$start+conditional.density$stop)/2
state.labs <- c("CTCF", "Enhancer", "Promoter flank", "Repressed", "Transcribed", "Promoter/TSS", "Weak Enhancer")
names(state.labs) <- c("CTCF", "E", "PF", "R", "T", "TSS", "WE")
p_labels = data.frame(feature = c("CTCF", "E", "PF", "R", "T", "TSS", "WE"), 
                      label = c("NS", "italic(P)==0.0018", "NS", "italic(P)==0.0013", "NS", "NS", "NS"))

#Figure S3B
cols <- brewer.pal(8,"Set2")
FigureS3B<- ggplot(conditional.density, aes(x=pos/1000000, y=count)) +
  geom_line() +
  theme_bw() +
  xlab("Distance to QTLs (Mb)") +
  ylab("#annotations/kb") +
  ggtitle("Conditional") +
  facet_grid(. ~ feature, labeller = labeller(feature = state.labs)) +
  theme(panel.spacing.x = unit(4, "mm")) +
  aes(fill = as.factor(feature), col=as.factor(feature)) + scale_fill_manual(values = cols[c(8,4,8,3,8,8,8)]) + scale_colour_manual(values = cols[c(8,4,8,3,8,8,8)]) +
  theme(
    strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
      size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
    axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold")) +
  geom_text(x=0, y=(max(conditional.density$count)+100), aes(label=label), data=p_labels, parse=T, inherit.aes=F, size = 4) +
  ylim(NA, (max(conditional.density$count)+100))

#Read in data for frequency plots representing the density of eQTL loci surrounding functional genomic tracks in unique NK cell eQTL.
unique.density <- read.table("FigureS3C.data.txt", header = F)
colnames(unique.density) <- c("start", "stop", "count", "feature")
unique.density$pos <- (unique.density$start+unique.density$stop)/2
state.labs <- c("CTCF", "Enhancer", "Promoter flank", "Repressed", "Transcribed", "Promoter/TSS", "Weak Enhancer")
names(state.labs) <- c("CTCF", "E", "PF", "R", "T", "TSS", "WE")
p_labels = data.frame(feature = c("CTCF", "E", "PF", "R", "T", "TSS", "WE"), 
                      label = c("NS", "NS", "NS", "italic(P)==0.00033", "italic(P)==0.0022", "italic(P)==2.20%*%10^-7", "italic(P)==0.00053"))

FigureS3C<- ggplot(unique.density, aes(x=pos/1000000, y=count)) +
  geom_line() +
  theme_bw() +
  xlab("Distance to QTLs (Mb)") +
  ylab("#annotations/kb") +
  ggtitle("Unique") +
  facet_grid(. ~ feature, labeller = labeller(feature = state.labs)) +
  theme(panel.spacing.x = unit(4, "mm")) +
  aes(fill = as.factor(feature), col=as.factor(feature)) + scale_fill_manual(values = cols[c(8,8,8,3,4,4,4)]) + scale_colour_manual(values = cols[c(8,8,8,3,4,4,4)]) +
  theme(
    strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
      size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
    axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold")) +
  geom_text(x=0, y=(max(unique.density$count)+100), aes(label=label), data=p_labels, parse=T, inherit.aes=F, size = 4) +
  ylim(NA, (max(unique.density$count)+100))

#read in data describing under/over-represenation TFBS in total primary NK cell eQTL
fenrich <- read.table("FigureS3D.data.txt", header = T)

#Figure S3D
FigureS3D <- ggplot(fenrich, aes(x = -log(p,10), y = or, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Primary") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,4,3)]) +
  scale_colour_manual(values = cols[c(8,4,3)])  +
  ylab("Odds ratio") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(fenrich, is_annotate==1), aes(label=name), size=4, col = c("black"))

#read in data describing under/over-represenation TFBS in total conditional NK cell eQTL
fenrich <- read.table("FigureS3E.data.txt", header = T)

#Figure S3E
FigureS3E <- ggplot(fenrich, aes(x = -log(p,10), y = or, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Conditional") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,4,3)]) +
  scale_colour_manual(values = cols[c(8,4,3)])  +
  ylab("Odds ratio") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(fenrich, is_annotate==1), aes(label=name), size=4, col = c("black"))

#read in data describing under/over-represenation TFBS in total conditional NK cell eQTL
fenrich <- read.table("FigureS3E.data.txt", header = T)

#Figure S3E
Figure3E <- ggplot(fenrich, aes(x = -log(p,10), y = or, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Conditional") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,4,3)]) +
  scale_colour_manual(values = cols[c(8,4,3)])  +
  ylab("Odds ratio") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(fenrich, is_annotate==1), aes(label=name), size=4, col = c("black"))

#read in data describing under/over-represenation TFBS in total unique NK cell eQTL
fenrich <- read.table("FigureS3F.data.txt", header = T)

#Figure S3F
FigureS3F <- ggplot(fenrich, aes(x = -log(p,10), y = or, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Unique") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,4,3)]) +
  scale_colour_manual(values = cols[c(8,4,3)])  +
  ylab("Odds ratio") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(fenrich, is_annotate==1), aes(label=name), size=4, col = c("black"))

#read in data comparing feature enrichment in conditional vs primary NK cell eQTL
c.vs.p <- read.table("FigureS3G.data.txt", header = T)

#Figure S3G
FigureS3G <- ggplot(c.vs.p, aes(x = -log(p,10), y = fc, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Conditional vs Primary") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,3,4)]) +
  scale_colour_manual(values = cols[c(8,3,4)])  +
  ylab("Fold change") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(c.vs.p, is_annotate==1), aes(label=feature), size=4, col = c("black"))

#read in data comparing feature enrichment in shared vs unique NK cell eQTL
s.vs.u <- read.table("FigureS3H.data.txt", header = T)

#Figure S3H
FigureS3H <- ggplot(s.vs.u, aes(x = -log(p,10), y = fc, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1) +
  geom_point(aes(size = 1, alpha = 0.72), pch = 21, show.legend = FALSE) +
  ggtitle("Unique vs Shared") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), plot.title=element_text(size=25, face = "bold"))+
  scale_fill_manual(values = cols[c(8,4,3)]) +
  scale_colour_manual(values = cols[c(8,4,3)])  +
  ylab("Fold change") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(s.vs.u, is_annotate==1), aes(label=feature), size=4, col = c("black"))

##Construct Figure S3
FigureS3 <- ((FigureS3A)/(FigureS3B)/(FigureS3C)/((FigureS3D|FigureS3E|FigureS3F))/((FigureS3G|FigureS3H))) + plot_annotation(tag_levels = 'A')


###Figure S4: GOBP enrichment of NK cell-specific eQTL.###
#read in list of genes with eQTL specific to NK cells
test <- as.character(read.table("FigureS4.test_data.txt", header = F)$V1)
#read in list of genes tested
background <- as.character(read.table("FigureS4.background_data.txt", header = F)$V1)

#test for enrichment using XGR
eTerm_GOBP <- xEnricherGenes(data=test, background=background, ontology="GOBP", ontology.algorithm="none", test = "hypergeo", RData.location=RData.location, min.overlap=7, p.adjust.method = "BH")
out.GOBP <- subset(xEnrichViewer(eTerm_GOBP, top_num=40, details = F), adjp<0.05)
out.GOBP$GOBP <- rownames(out.GOBP)
out.GOBP$GOBP <- factor(out.GOBP$GOBP, levels = c(row.names(out.GOBP[order(out.GOBP$fc), ])))

#Figure S4
FigureS4 <- ggplot(out.GOBP, aes(x = GOBP, y = fc)) +
  geom_point(aes(fill = -log10(pvalue), size = nOverlap), pch = 21) +
  scale_fill_distiller(palette = "Spectral") +
  expand_limits(y=c(0,15)) +
  scale_y_continuous(breaks=c(-10, -5,0,5,10)) +
  ylab("Fold change") +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16)) +
  coord_flip() +
  geom_text(aes(label = name),
            size = 4.3, hjust = 0, nudge_y=0.25) +
  scale_color_manual(values = c("dodgerblue", "firebrick2"))

###Figure S5: Effect of background GWAS trait on NK cell eQTL GWAS enrichment.###
#read in data comparing enrichment of NK cell eQTL for GWAS traits according to background (height vs. total body impedence)
enrich.background <- read.table("FigureS5.data.txt", header = T)
#correlation between fold changes observed using either background
cor.test(enrich.background$height_fc, enrich.background$imped_fc)

#Figure S5
FigureS5 <- ggplot(for_cor, aes(x = imped_fc, y = fc, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_point(aes(size = pos, alpha = 0.72), pch = 21, show.legend = FALSE) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  scale_size(range = c(5,10)) +
  scale_fill_manual(values = cols[c(8,4)]) +
  scale_colour_manual(values = cols[c(8,4)])  +
  ylab("Fold change (height)") +
  xlab("Fold change (body impedence)") +
  coord_flip() +
  xlim(0,4) +
  ylim(0,4) +
  annotate("text", x = 0.5, y = 3.75, label = "r=0.575", size=6) +
  geom_text_repel( data=subset(for_cor, is_annotate==1), aes(label=label), size=4, col = c("black"))

###Figure S6: NK cell eQTL at ERAP2.###
#read in ERAP2 cis eQTL data
erap2.cis <- read.table("FigureS6A.data.txt", header = T)

#Figure S6A
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.33%*%10^-76"))
cols <- brewer.pal(8,"Set2")
p1 = ggplot(erap2.cis, aes(x=factor(round(rs1363974.geno,0)), y=ERAP2.exp.PCcorr))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.05, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("AA", "AG", "GG")) + aes(fill = factor(round(rs1363974.geno,0)), col=factor(round(rs1363974.geno,0))) + scale_fill_manual(values = cols[c(2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2)])
p5 = p4 + ylab(expression(italic(ERAP2)~ expression)) + xlab("rs1363974")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
FigureS6A = p6 + geom_text(x=2, y=2, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#read in forest plot for erap2 eQTL across immune cells
forest.erap2 <- read.table("FigureS6B.data.txt", header = T)
#construct data.frame for forest plot
label <- c("NK", "Neut", "Mono", "CD4", "CD8")
mean  <- forest.erap2$BETA
lower <- forest.erap2$BETA-(1.96*forest.erap2$SE)
upper <- forest.erap2$BETA+(1.96*forest.erap2$SE)
df <- data.frame(label, mean, lower, upper)

#Figure S6B
df$label <- factor(df$label, levels=rev(df$label))
FigureS6B <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=1, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = cols[c(2,2,2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2,2,2)]) +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() + ylim(-0.5,2.0)  +
  ylab("beta") + annotate("text", x = 5.5, y = 0.75, label = "rs1363974:ERAP2", size=5) + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15))

#read in data to plot colocalisation at erap2
co.plot <- read.table("FigureS6CD.data.txt", header = T)
#divide R2 into bins
co.plot$bin_r2 <- 1
co.plot$bin_r2[which(co.plot$r2>0.2 & co.plot$r2 <= 0.5)] <- 2
co.plot$bin_r2[which(co.plot$r2>0.5 & co.plot$r2 <= 0.8)] <- 3
co.plot$bin_r2[which(co.plot$r2>0.8)] <- 4

#Figure S6C: neutrophil count GWAS: ERAP2 eQTL colocalisation plot
cols2 <- brewer.pal(11,"Spectral")
FigureS6C <- ggplot(co.plot, aes(x=-log10(neut.gwas.pval), y=-log10(eqtl.pval))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,3,1)])) +
  ylab("ERAP2 eQTL (-logP)") + 
  xlab("Neutrophil % GWAS (-logP)") + 
  annotate("text", x = 0, y = 82.5, label = "PP4=0.91", size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) # use a white background

#Figure S6D: lymphocyte count GWAS: ERAP2 eQTL colocalisation plot
FigureS6D <- ggplot(co.plot, aes(x=-log10(lymph.gwas.pval), y=-log10(eqtl.pval))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,3,1)])) +
  ylab("ERAP2 eQTL (-logP)") + 
  xlab("Lymphocyte % GWAS (-logP)") + 
  annotate("text", x = 0, y = 82.5, label = "PP4=0.94", size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) # use a white background

#Construct FigureS6
FigureS6 <- ((FigureS6A|FigureS6B)/(FigureS6C|FigureS6D))

###Figure S8: Study power.###
#cis analysis - varying nTests according to the effective number of indepdendent SNPs per cis window.
#e.g. 10 SNPs per window assumes 180,000 tests given 18,000 genes tested
mafs <- seq(0,0.5,0.01)
#1 SNP per window
for (i in c(1:51)){
  MAF=mafs[i]
  power.cis[i] <- powerEQTL.SLR(
    MAF=MAF,
    slope = 0.13,
    n = 245,
    power = NULL,
    sigma.y = 0.13,
    FWER = 0.05,
    nTests = 18000)[[1]]
}
for_plot1 <- data.frame(cbind(mafs, power.cis, "cis - 1"))
#10SNPs per window
power.cis <- c()
for (i in c(1:51)){
  MAF=mafs[i]
  power.cis[i] <- powerEQTL.SLR(
    MAF=MAF,
    slope = 0.13,
    n = 245,
    power = NULL,
    sigma.y = 0.13,
    FWER = 0.05,
    nTests = 180000)[[1]]
}
for_plot10 <- data.frame(cbind(mafs, power.cis, "cis - 10"))
#100 SNPs per window
for (i in c(1:51)){
  MAF=mafs[i]
  power.cis[i] <- powerEQTL.SLR(
    MAF=MAF,
    slope = 0.13,
    n = 245,
    power = NULL,
    sigma.y = 0.13,
    FWER = 0.05,
    nTests = 1800000)[[1]]
}
for_plot100 <- data.frame(cbind(mafs, power.cis, "cis - 100"))
#1000 SNPs per window
for (i in c(1:51)){
  MAF=mafs[i]
  power.cis[i] <- powerEQTL.SLR(
    MAF=MAF,
    slope = 0.13,
    n = 245,
    power = NULL,
    sigma.y = 0.13,
    FWER = 0.05,
    nTests = 18000000)[[1]]
}
for_plot1000 <- data.frame(cbind(mafs, power.cis, "cis - 1000"))
#5000 SNPs per window
for (i in c(1:51)){
  MAF=mafs[i]
  power.cis[i] <- powerEQTL.SLR(
    MAF=MAF,
    slope = 0.13,
    n = 245,
    power = NULL,
    sigma.y = 0.13,
    FWER = 0.05,
    nTests = 90000000)[[1]]
}
for_plot5000 <- data.frame(cbind(mafs, power.cis, "cis - 5000"))
#in trans we assume a MAF of 10% and 1,000,000 independent SNPs genome-wide
powerEQTL.SLR(
  MAF=0.1,
  slope = 2,
  n = 245,
  power = NULL,
  sigma.y = 2,
  FWER = 0.05,
  nTests = 18000000000)
for_plot.trans <- data.frame(cbind(mafs, power.cis, "trans"))

#concatenate outputs
for_plot <- rbind(for_plot1, for_plot10, for_plot100, for_plot1000, for_plot5000, for_plot.trans)
for_plot$mafs <- as.numeric(as.character(for_plot$mafs))
for_plot$power.cis <- as.numeric(as.character(for_plot$power.cis))
for_plot$V3 <-factor(for_plot$V3)

#Figure S8
cols <- brewer.pal(11, "Spectral")
FigureS8 = ggplot(for_plot, aes(mafs, power.cis)) + 
  geom_line(aes(colour = factor(V3)), size=2) +
  theme_bw() +
  xlim(0, 0.2) +
  scale_fill_manual(values = cols[c(5:1,10)]) + 
  scale_colour_manual(values = cols[c(5:1,10)]) +
  labs(colour='Independent SNPs \nper testing window') +
  ylab("Power") +
  xlab("MAF") +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0)) +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=15), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12))
FigureS8

###Figure S9: Study sample population structure.###
#read in data
pca <- read.table("FigureS9.data.txt", header = T)
levels(pca$super_pop)[5] <- "NK eQTL samples"

#Figure S9
cols <- brewer.pal(12, "Set3")
FigureS9 = ggplot(na.omit(pca), aes(PC2, PC1)) + 
  geom_point(aes(colour = factor(super_pop)), size = 3) +
  scale_fill_manual(values = cols[c(1:6)]) + 
  scale_colour_manual(values = cols[c(1:6)]) +
  labs(colour='Population') +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))



