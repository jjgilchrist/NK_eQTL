library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

cols <- brewer.pal(8,"Set2")

#construct dataframe of cis eQTL numbers (primary and conditional) - adding  the proportion that colocalise with at least one GWAS trait
df <- data.frame(eqtl <- rep(c("primary", "secondary","tertiary","quaternary"), each = 2),
rtc <- rep(c("no", "yes"), 4),
count <- c((3951-1242), 1242, (528-158), 158, (63-16), 16, 3,0)
)

df$eqtl <- factor(df$eqtl, levels = c("primary", "secondary", "tertiary", "quaternary"))

annotation <- data.frame(
  x = c(1,2,3,4),
  y = c(3951+200,528+200,63+200,3+200),
  label = c("1,242/3,951", "158/528", "16/63", "0/3")
)

#Panel 1A: number of cis eQTL by conditional rank and the proportion with a colocalised (RTC GWAS trait
Figure1A <- ggplot(df, aes(x = eqtl, y = count))+
  geom_col(aes(fill = rtc), width = 0.7)+
  scale_fill_manual(values = cols[c(8,4)]) +
  scale_colour_manual(values = cols[c(8,4)]) +
  theme_bw() + xlab("Conditional eQTL Rank") +
  theme(legend.position = "none", axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_text(data=annotation, aes( x=x, y=y, label=label),
            color=cols[c(4)], 
            size=5 , fontface="bold" )

#Read in eQTL sharing estimates from moloc
moloc.out <- read.table("Figure1B.data.txt", header = F)
moloc.out$sig <- 0
moloc.out$sig[which(moloc.out$V3<0.2)] <- 1
moloc.out$sig[which(moloc.out$V3>0.8)] <- 2

cols <- brewer.pal(8,"Set2")

#Panel 1B: Distribution of cis eQTL unique to NK cells or shared with at least one of neutrophils, monocytes, CD4+ and CD8+ T cells.
Figure1B<-ggplot(moloc.out, aes(x=(1-V3), fill=factor(sig), colour=factor(sig))) + 
  geom_histogram(binwidth=0.012, alpha=0.72) +
  scale_fill_manual(values = cols[c(8,1,2)]) +
  scale_colour_manual(values = cols[c(8,1,2)]) +
  theme_bw() + xlab("Probability eQTL is specific to NK cells") + theme(legend.position = "none") +
  annotate("text", x = -0.02, y = 775, label = "Shared", size=10, hjust=0) + ylim(NA,800) +
  annotate("text", x = 1.02, y = 775, label = "Unique", size=10, hjust=1) +
  geom_hline(yintercept=0, color=cols[8]) +
theme(axis.text=element_text(size=12),
      axis.title=element_text(size=14))

#Read in total cis eQTL:GWAS enrichment results (as estimated by coloc)
#TABLE HEADERS:
#trait/label = GWAS tested
#coloc = #cis eQTL colocalising with associations for that trait, 
#pos = #cis eQTL within a 250kb window of a trait association
#h_coloc = 	#cis eQTL colocalising with associations for the background trait (height)
#h_pos = 	#cis eQTL cwithin a 250kb window of an association for the background trait (height)
#fc = fold change
#p = p-value (Fisher's exact)
total_enrich <- read.table("Figure1Ca.data.txt", header = T)
total_enrich$p.adj <- p.adjust(total_enrich$p, method = "fdr")

#Read in cis eQTL (unique to NK cells):GWAS enrichment results (as estimated by coloc)
#TABLE HEADERS: as for total_enrich
unique_enrich <- read.table("Figure1Cb.data.txt", header = T)
unique_enrich$p.adj <- p.adjust(unique_enrich$p, method = "fdr")

#Prepare data for facet plot
unique_enrich <- unique_enrich[order(unique_enrich$p), ]
unique_enrich$threshold <- c(rep(1,4),rep(0,54))
unique_enrich$label <- as.character(unique_enrich$label)
unique_enrich$is_annotate <- c(rep(1,4),rep(0,54))

total_enrich <- total_enrich[order(total_enrich$p), ]
total_enrich$threshold <- c(rep(1,8),rep(0,65))
total_enrich$label <- as.character(total_enrich$label)
total_enrich$is_annotate <- c(rep(1,8),rep(0,65))

total_enrich$expt <- "total"
unique_enrich$expt <- "unique"
facet_enrich <- rbind(total_enrich, unique_enrich)

Figure1C <- ggplot(facet_enrich, aes(x = -log(p,10), y = fc, colour=as.factor(threshold), fill=as.factor(threshold))) +
  geom_point(aes(size = pos, alpha = 0.72), pch = 21, show.legend = FALSE) +
  facet_grid(. ~ expt) +
  theme_bw() +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
    axis.title=element_text(size=14))+
  scale_size(range = c(5,10)) +
  scale_fill_manual(values = cols[c(8,4)]) +
  scale_colour_manual(values = cols[c(8,4)])  +
  ylab("Fold change") +
  xlab("-log10(P)") +
  coord_flip() +
  geom_text_repel( data=subset(facet_enrich, is_annotate==1), aes(label=label), size=4, col = c("black"))


#Read in CD226 expression data (corrected for 32 PCs of expression data) and rs1788098 genotype doses
cd226 <- read.table("Figure1D.data.txt", header = T)
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.92%*%10^-19"))

#Panel 1D: Effect of rs1788098 genotype on CD226 expression in NK cells.
p1 = ggplot(cd226, aes(x=factor(round(rs1788098,0)), y=cd226.exp.pc32_corr))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.015, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(round(rs1788098,0)), col=factor(round(rs1788098,0))) + scale_fill_manual(values = cols[c(1,1,1)]) + scale_colour_manual(values = cols[c(1,1,1)])
p5 = p4 + ylab(expression(italic(CD226)~ expression)) + xlab("rs1788098")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + ylim(NA, 0.5)
Figure1D = p6 + geom_text(x=2, y=0.45, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#Read in estimates of effect of rs1788098 genotype on CD226 expression across NK cell, neutrophils, monocytes, CD4+ & CD8+ T cells
forest <- read.table("Figure1E.data.txt", header = T, row.names = 1)

#Calculate 95% confidence intervals
label <- c("NK", "Neut", "Mono", "CD4", "CD8")
mean  <- forest$BETA
lower <- forest$BETA-(1.96*forest$SE)
upper <- forest$BETA+(1.96*forest$SE)

df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=rev(df$label))

#Panel 1E: The effect of rs1788098 genotype on CD226 expression is specific to NK cells. Significant eQTL effects are highlighted (green).
Figure1E <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=5, size = 2) + aes(fill = label, col=label) + scale_fill_manual(values = cols[c(8,8,8,8,1)]) + scale_colour_manual(values = cols[c(8,8,8,8,1)]) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("beta") + annotate("text", x = 5.5, y = -0.08, label = "rs1788098:CD226", size=5) + scale_x_discrete(expand = expand_scale(add = 1)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) 

#Read in data to compare local eQTL association at CD226 with DM1 GWAS
co.plot <- read.table("Figure1F.data.txt",header = T)

#Divide r2 into bins
co.plot$bin_r2 <- 1
co.plot$bin_r2[which(co.plot$r2>0.2 & co.plot$r2 <= 0.5)] <- 2
co.plot$bin_r2[which(co.plot$r2>0.5 & co.plot$r2 <= 0.8)] <- 3
co.plot$bin_r2[which(co.plot$r2>0.8)] <- 4

#Panel 1F: The CD226 eQTL in NK cells colocalises with a risk locus for type-1 diabetes. SNPs are coloured according to strength of LD (CEU population) to the peak eSNP (rs1788098).
cols2 <- brewer.pal(11,"Spectral")
Figure1F <- ggplot(co.plot, aes(x=-log10(gwas_p), y=-log10(eqtl_p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,3,1)])) +
  ylab("CD226 eQTL (-logP)") + 
  xlab("Type 1 Diabetes GWAS (-logP)") + 
  annotate("text", x = 0, y = 18, label = "PP4=0.97", size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12))

#Read in data describing rs1788097 (proxy for rs1788098) association with GWAS of autoimmune diseases and haematological indices
forest.traits <- read.table("Figure1G.data.txt", header = T)

#calculate 95% confidence intervals
label <- forest.traits$trait
mean  <- forest.traits$beta
lower <- forest.traits$beta-(1.96*forest.traits$se)
upper <- forest.traits$beta+(1.96*forest.traits$se)
df <- data.frame(label, mean, lower, upper)
df$label <- factor(df$label, levels=rev(df$label))

#Panel 1G: Association of rs1788097 (exact proxy for rs1788098 in European populations, r2 = 1) with autoimmune diseases and haematological indices. GWAS-significant (p < 5e-8) associations are highlighted (pink).
Figure1G <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=1, size = 1.5) + aes(fill = label, col=label) + scale_fill_manual(values = cols[rev(c(4,4,4,4,8,8,8,8,4,4,4,4,4,4,4,8,8,8,4,8))]) + scale_colour_manual(values = cols[rev(c(4,4,4,4,8,8,8,8,4,4,4,4,4,4,4,8,8,8,4,8))]) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("beta/log(OR)") + annotate("text", x = 20, y = 0.05, label = "rs1788097", size=6, hjust=0.5) + scale_x_discrete(expand = expand_scale(add = 1.5)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=10),
                     axis.title=element_text(size=15)) # use a white background

#Read in data to compare local unconditionned eQTL association at IRF5 with SLE GWAS
co.plot <- read.table("Figure1HI.data.txt", header = T)

#Divide r2 into bins
co.plot$bin_r2 <- 1
co.plot$bin_r2[which(co.plot$r2>0.2 & co.plot$r2 <= 0.5)] <- 2
co.plot$bin_r2[which(co.plot$r2>0.5 & co.plot$r2 <= 0.8)] <- 3
co.plot$bin_r2[which(co.plot$r2>0.8)] <- 4

#highlight rs7796963
co.plot <- co.plot[order(co.plot$eqtl_p),]
co.plot$annotate <- c(1, (rep(0,length(co.plot$eqtl_p)-1)))

#Panel 1H: Regional association plot of the primary IRF5 eQTL in NK cells.
Figure1H <- ggplot(co.plot, aes(x=bp/1000, y=-log10(eqtl_p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,3,1)])) +
  ylab("IRF5 primary eQTL (-logP)") + 
  xlab("Chr 7 BP (kb)") + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  geom_text_repel( data=subset(co.plot, annotate==1), aes(label=snp), size=6, col = c("black"), min.segment.length = unit(0, 'lines'),
                   nudge_y = 8) +
  ylim(NA, 90)

#Panel 1I: Plot demonstrating that the primary IRF5 eQTL does not colocalise with a GWAS risk locus in the IRF5 region for systemic lupus erythematosus.
Figure1I <- ggplot(co.plot, aes(x=-log10(gwas_p), y=-log10(eqtl_p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,3,1)])) +
  ylab("IRF5 primary eQTL (-logP)") + 
  xlab("SLE GWAS (-logP)") + 
  annotate("text", x = 0, y = 82.5, label = "PP4=0.00", size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) 

#Read in data to compare local eQTL association at IRF5 (conditionned on rs7796963) with SLE GWAS
co.plot <- read.table("Figure1JK.data.txt", header = T)

#Divide r2 into bins
co.plot$bin_r2 <- 1
co.plot$bin_r2[which(co.plot$r2>0.2 & co.plot$r2 <= 0.5)] <- 2
co.plot$bin_r2[which(co.plot$r2>0.5 & co.plot$r2 <= 0.8)] <- 3
co.plot$bin_r2[which(co.plot$r2>0.8)] <- 4

#Panel 1J: Conditionning on the peak IRF5 eSNP reveals an independent, secondary eQTL for IRF5 in NK cells.
Figure1J <- ggplot(co.plot, aes(x=bp/1000, y=-log10(eqtl_p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,1)])) +
  ylab("IRF5 secondary eQTL (-logP)") + 
  xlab("Chr 7 BP (kb)") + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  annotate("text", x = 128400, y = 14, label = "Conditioned: rs7796963", size=5, hjust=0) +
  ylim(NA, 14)

#Panel 1K: The secondary IRF5 eQTL does colocalise with the IRF5 association for systemic lupus erythematosus.
Figure1K <- ggplot(co.plot, aes(x=-log10(gwas_p), y=-log10(eqtl_p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[8],cols2[c(5,1,1)])) + scale_colour_manual(values = c(cols[8],cols2[c(5,1,1)])) +
  ylab("IRF5 secondary eQTL (-logP)") + 
  xlab("SLE GWAS (-logP)") + 
  annotate("text", x = 0, y = 14.5, label = "PP4=0.95", size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12))

Figure1 <- (((Figure1A|Figure1B|Figure1C))/((Figure1D|Figure1E|Figure1F|Figure1G))/((Figure1H|Figure1I|Figure1J|Figure1K)))+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))

