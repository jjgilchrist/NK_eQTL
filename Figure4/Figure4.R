library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
cols <- brewer.pal(8,"Set2")

#read in protein phenotype data
protein.data <- read.table("Figure4.protein_data.txt", header = T)
protein.data$facet <- "protein"
#read in RNA phenotype data: expression is corrected for 32 PCs
rna.data<- read.table("Figure4.rna_data.txt", header = T)

#logit transform FACS proportions
protein.data$CD57_logit <- log((protein.data$CD57_pos/100)/(1-(protein.data$CD57_pos/100)))
protein.data$CD226_logit <- log((protein.data$CD226_pos/100)/(1-(protein.data$CD226_pos/100)))
protein.data$KIR_logit <- log((protein.data$panKIR/100)/(1-(protein.data$panKIR_pos/100)))

#protein:genotype correlations
summary(lm(protein.data$CD57_logit~protein.data$rs77478906))
summary(lm(protein.data$CD226_logit~protein.data$rs1788098))
summary(lm(protein.data$KIR_logit~protein.data$kir2ds4del))

#Figure 4 B3GAT1/CD57 RNA eQTL
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.90%*%10^-11"))
rna.data$facet <- "RNA"
facet <- "RNA/protein"
p1 = ggplot(rna.data, aes(x=factor(rs77478906), y=ENSG00000109956))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TC", "CC")) + aes(fill = factor(round(rs77478906,0)), col=factor(round(rs77478906,0))) + scale_fill_manual(values = cols[c(2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2)])
p5 = p4 +ylab(expression(italic(B3GAT1)~ expression)) + xlab("rs77478906")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + ylim(NA, 1.1)
p7 <- p6+ facet_wrap( ~ facet, ncol=1, scales="free") +theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
  size = 15, color = "white", face = "bold"))
Figure4.CD57.RNA = p7 + geom_text(x=2, y=1.05, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#Figure 4 B3GAT1/CD57 protein QTL
p_labels = data.frame(expt = c("facet"), label = c("italic(P)==3.85%*%10^-5"))
Figure4.CD57.protein = ggplot(protein.data, aes(x=factor(round(rs77478906,0)), y=CD57_pos)) +
  geom_dotplot(binaxis="y", binwidth=1.2, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TC", "CC")) + aes(fill = factor(round(rs77478906,0)), col=factor(round(rs77478906,0))) + scale_fill_manual(values = cols[c(2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2)]) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) +
  ylab("%CD57+ NK cells") + xlab("rs77478906") +
  facet_wrap( ~ facet, ncol=2, scales="free") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold")) + ylim(NA, 80) +
  geom_text(x=2, y=78, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)
Figure4.CD57.protein


#Figure 4 CD226 RNA eQTL
p_labels = data.frame(expt = c("base"), label = c("italic(P)==5.32%*%10^-21"))
p1 = ggplot(rna.data, aes(x=factor(rs1788098), y=ENSG00000150637))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.015, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(round(rs1788098,0)), col=factor(round(rs1788098,0))) + scale_fill_manual(values = cols[c(3,3,3)]) + scale_colour_manual(values = cols[c(3,3,3)])
p5 = p4 +ylab(expression(italic(CD226)~ expression)) + xlab("rs1788098")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + ylim(NA, 0.5)
p7 <- p6+ facet_wrap( ~ facet, ncol=1, scales="free") +theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
  size = 15, color = "white", face = "bold"))
Figure4.CD226.RNA = p7 + geom_text(x=2, y=0.475, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)
Figure4.CD226.RNA

#Figure 4 CD226 protein QTL
p_labels = data.frame(expt = c("facet"), label = c("italic(P)==0.0301"))
Figure4.CD226.protein = ggplot(protein.data[-which(is.na(protein.data$CD226_pos)),], aes(x=factor(round(rs1788098,0)), y=CD226_pos)) +
  geom_dotplot(binaxis="y", binwidth=1.2, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(round(rs1788098,0)), col=factor(round(rs1788098,0))) + scale_fill_manual(values = cols[c(3,3,3)]) + scale_colour_manual(values = cols[c(3,3,3)]) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) +
  ylab("%CD226+ NK cells") + xlab("rs1788098") +
  facet_wrap( ~ facet, ncol=2, scales="free") +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold")) + ylim(NA, 100) +
  geom_text(x=2, y=98, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)
Figure4.CD226.protein

#Figure 4 KIR2D protein QTL
p_labels = data.frame(expt = c("facet"), label = c("italic(P)==0.0459"))
Figure4.kir.protein = ggplot(protein.data, aes(x=factor(round(kir2ds4del,0)), y=panKIR_pos)) +
  geom_dotplot(binaxis="y", binwidth=1.5, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("0", "1", "2")) + aes(fill = factor(round(kir2ds4del,0)), col=factor(round(kir2ds4del,0))) + scale_fill_manual(values = cols[c(4,4,4)]) + scale_colour_manual(values = cols[c(4,4,4)]) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) +
  ylab("%KIR2D+ NK cells") + xlab("KIR2DS4del") +
  ylim(NA, 85) +
  geom_text(x=2, y=83.5, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)
Figure4.kir.protein

#Construct Figure 4
Figure4 <- (Figure4.CD57.RNA|Figure4.CD57.protein|Figure4.CD226.RNA|Figure4.CD226.protein|Figure4.kir.protein) + plot_layout(widths = c(1,1,1,1,2))

