library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(patchwork)
library(ggbio)
library(data.table)
library(GenomicRanges)


#Read in data for trans manhattan plot
trans.nom.less <- read.table("Figure2Aa.data.txt", header = T)

#Order chromosomes
trans.nom.less$chr <- factor(trans.nom.less$chr, levels = c("chr1", 
                                                  "chr2", 
                                                  "chr3", 
                                                  "chr4", 
                                                  "chr5", 
                                                  "chr6", 
                                                  "chr7", 
                                                  "chr8", 
                                                  "chr9", 
                                                  "chr10", 
                                                  "chr11", 
                                                  "chr12", 
                                                  "chr13", 
                                                  "chr14", 
                                                  "chr15", 
                                                  "chr16", 
                                                  "chr17", 
                                                  "chr18", 
                                                  "chr19", 
                                                  "chr20", 
                                                  "chr21", 
                                                  "chr22", 
                                                  "chrX", 
                                                  "chrY"))

levels(trans.nom.less$chr) <- c(1:22, "X", "Y")

#Read in significant loci to highlight:
trans.fdr <- read.table("Figure2Ab.data.txt", header = T)

#identify significant (FDR<0.05) loci
trans.nom.less$is.sig <- 0
trans.fdr$is.sig <- 1
trans.nom.less$is.sig <- trans.fdr$is.sig[match(trans.nom.less$pair, trans.fdr$pair)]

#identify cis mediators
trans.nom.less$cis.mediator <- NA
trans.nom.less$cis.mediator <- trans.fdr$cis_mediator[match(trans.nom.less$pair, trans.fdr$pair)]

#identify trans gene for labels
trans.nom.less$trans_gene <- NA
trans.nom.less$trans_gene <- trans.fdr$trans_gene[match(trans.nom.less$pair, trans.fdr$pair)]

#identify peak eSNP
trans.nom.less$min_for_peak <- NA
trans.nom.less$min_for_peak <- trans.fdr$min_for_peak[match(trans.nom.less$pair, trans.fdr$pair)]

trans.nom.less$cis.mediator <- as.character(trans.nom.less$cis.mediator)
trans.nom.less$trans_gene <- as.character(trans.nom.less$trans_gene)

don <- trans.nom.less %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(trans.nom.less, ., by=c("chr"="chr")) %>%
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Panel 2A: Manhattan plot depicting signi cant trans eQTL in NK cells. Physical position (x axis) represents location of the target gene.
Figure2A <- ggplot(don, aes(x=BPcum, y=-log10(p))) +
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 24 )) +
  # Add highlighted points
  geom_point(data=subset(don, cis.mediator=="other"), color=cols[1], size=1.8) +
  geom_point(data=subset(don, cis.mediator=="MC1R"), color=cols[2], size=1.8) +
  geom_point(data=subset(don, cis.mediator=="GNLY"), color=cols[3], size=1.8) +
  geom_point(data=subset(don, cis.mediator=="UVSSA"), color=cols[6], size=1.8) +
  geom_point(data=subset(don, cis.mediator=="KIR"), color=cols[4], size=1.8) +
  # Add labels
  geom_text_repel( data=subset(subset(don, cis.mediator=="MC1R" | cis.mediator=="GNLY" | cis.mediator=="UVSSA" | cis.mediator=="KIR"), min_for_peak==1), aes(label=trans_gene), size=5, min.segment.length = unit(0, 'lines'), nudge_y = 5) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("chromosome") +
  # Custom the theme:
  theme_bw() + ylim(NA, 110) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=12),
    axis.title=element_text(size=15)
  )

#Read in GNLY trans eQTL data
gnly <- read.table("Figure2BD.data.txt", header = T)

p_labels = data.frame(expt = c("base"), label = c("italic(P)==8.86%*%10^-23"))

#Panel 2B: Effect of rs1866140 genotype on GNLY expression in NK cells.
p1 = ggplot(gnly, aes(x=factor(round(geno,0)), y=resid))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.015, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("GG", "GT", "TT")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(3,3,3)]) + scale_colour_manual(values = cols[c(3,3,3)])
p5 = p4 + ylab(expression(italic(GNLY)~ expression)) + xlab("rs1866140") + ylim(NA,0.75)
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
Figure2B = p6 + geom_text(x=2, y=0.55, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#construct dataframe for facet plot
for_graph <- as.data.frame(rbind(cbind(gnly$geno, gnly$resid.N4BP2, "N4BP2"), 
                   cbind(gnly$geno, gnly$resid.KIAA1586, "KIAA1586"), 
                   cbind(gnly$geno, gnly$resid.ARHGAP30, "ARHGAP30"), 
                   cbind(gnly$geno, gnly$resid.EDNRB, "EDNRB"),
                   cbind(gnly$geno, gnly$resid.NOP56, "NOP56")))
colnames(for_graph) <- c("geno", "expn", "gene")
for_graph$geno <- as.numeric(as.character(for_graph$geno))
for_graph$expn <- as.numeric(as.character(for_graph$expn))

p_labels = data.frame(gene = c("N4BP2", "KIAA1586", "ARHGAP30", "EDNRB", "NOP56"), 
                      label1 = c("italic(P)==9.24%*%10^-22", "italic(P)==4.86%*%10^-12", "italic(P)==1.86%*%10^-12", "italic(P)==6.01%*%10^-15", "italic(P)==3.74%*%10^-14"),
                      label2 = c("PP4=1.00", "PP4=1.00", "PP4=1.00", "PP4=1.00", "PP4=0.97"))

#Panel 2D: Effect of rs1866140 genotype on GNLY regulatory network genes in trans.
p1 = ggplot(for_graph, aes(x=factor(round(geno,0)), y=expn))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.01, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("GG", "GT", "TT")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(3,3,3,3,3)]) + scale_colour_manual(values = cols[c(3,3,3,3,3)])
p5 = p4 + ylab("Trans expression") + xlab("rs1866140")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
p7 = p6 + facet_grid(.~gene) + geom_text(x=2, y=0.5, aes(label=label1), data=p_labels, parse=TRUE, inherit.aes=F, size = 4)
p8 = p7 + theme(
  strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
  axis.title=element_text(size=14))
Figure2D = p8 + geom_text(x=2, y=-0.75, aes(label=label2), data=p_labels, parse=F, inherit.aes=F, size = 4)



#Read in MC1R trans eQTL data
mc1r <- read.table("Figure2EG.data.txt", header = T)

p_labels = data.frame(expt = c("base"), label = c("italic(P)==5.80%*%10^-35"))

#Panel 2E: Effect of rs1866140 genotype on GNLY regulatory network genes in trans.
p1 = ggplot(mc1r, aes(x=factor(round(geno,0)), y=resid))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.015, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TC", "CC")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2)])
p5 = p4 + ylab(expression(italic(MC1R)~ expression)) + xlab("rs11642267")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + ylim(NA, 1.7)
Figure2E = p6 + geom_text(x=2, y=1.6, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#construct dataframe for facet plot
for_graph <- as.data.frame(rbind(cbind(mc1r$geno, mc1r$resid.DRD3, "DRD3"), 
                                 cbind(mc1r$geno, mc1r$resid.SNORD85, "SNORD85"), 
                                 cbind(mc1r$geno, mc1r$resid.FAM169B, "FAM169B"), 
                                 cbind(mc1r$geno, mc1r$resid.JADE1, "JADE1"),
                                 cbind(mc1r$geno, mc1r$resid.PDHA2, "PDHA2")))
colnames(for_graph) <- c("geno", "expn", "gene")
for_graph$geno <- as.numeric(as.character(for_graph$geno))
for_graph$expn <- as.numeric(as.character(for_graph$expn))

p_labels = data.frame(gene = c("DRD3", "SNORD85", "FAM169B", "JADE1", "PDHA2"), 
                      label1 = c("italic(P)==1.54%*%10^-34", "italic(P)==1.49%*%10^-25", "italic(P)==1.47%*%10^-17", "italic(P)==2.18%*%10^-20", "italic(P)==1.35%*%10^-11"),
                      label2 = c("PP4=0.91", "PP4=0.89", "PP4=0.87", "PP4=0.60", "PP4=0.87"))

#Panel 2G: Effect of rs117406136 genotype on MC1R regulatory network genes in trans.
p1 = ggplot(for_graph, aes(x=factor(round(geno,0)), y=expn))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.0075, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TC", "CC")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(2,2,2,2,2)]) + scale_colour_manual(values = cols[c(2,2,2,2,2)])
p5 = p4 + ylab("Trans expression") + xlab("rs11642267")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
p7 = p6 + facet_grid(.~gene) + geom_text(x=2, y=0.575, aes(label=label1), data=p_labels, parse=TRUE, inherit.aes=F, size = 4)
p8 = p7 + theme(
  strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
  axis.title=element_text(size=14)) + ylim(-0.4, 0.6)
Figure2G = p8 + geom_text(x=2, y=-0.35, aes(label=label2), data=p_labels, parse=F, inherit.aes=F, size = 4)

#Read in UVSSA trans eQTL data
uvssa <- read.table("Figure2HJ.data.txt", header = T)

p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.65%*%10^-40"))

#Panel 2H: Effect of rs111632154 genotype on UVSSA expression in NK cells.
p1 = ggplot(uvssa, aes(x=factor(round(geno,0)), y=resid))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.015, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TG", "GG")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(6,6,6)]) + scale_colour_manual(values = cols[c(6,6,6)])
p5 = p4 + ylab(expression(italic(UVSSA)~ expression)) + xlab("rs111632154")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + ylim(NA, 1.7)
Figure2H = p6 + geom_text(x=2, y=1.6, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)

#construct dataframe for facet plot
for_graph <- as.data.frame(rbind(cbind(uvssa$geno, uvssa$resid.SPTBN4, "SPTBN4"), 
                                 cbind(uvssa$geno, uvssa$resid.AVP, "AVP"), 
                                 cbind(uvssa$geno, uvssa$resid.SIM2, "SIM2"), 
                                 cbind(uvssa$geno, uvssa$resid.PINLYP, "PINLYP")))

colnames(for_graph) <- c("geno", "expn", "gene")
for_graph$geno <- as.numeric(as.character(for_graph$geno))
for_graph$expn <- as.numeric(as.character(for_graph$expn))

p_labels = data.frame(gene = c("SPTBN4", "AVP", "SIM2", "PINLYP"), 
                      label1 = c("italic(P)==4.48%*%10^-26", "italic(P)==2.00%*%10^-23", "italic(P)==3.69%*%10^-12", "italic(P)==2.87%*%10^-11"),
                      label2 = c("PP4=1.00", "PP4=1.00", "PP4=1.00", "PP4=1.00"))

#Panel 2J: Effect of rs111632154 genotype on UVSSA regulatory network genes in trans.
p1 = ggplot(for_graph, aes(x=factor(round(geno,0)), y=expn))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.0075, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TG", "GG")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[c(6,6,6,6)]) + scale_colour_manual(values = cols[c(6,6,6,6)])
p5 = p4 + ylab("Trans expression") + xlab("rs111632154")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
p7 = p6 + facet_grid(.~gene) + geom_text(x=2, y=0.575, aes(label=label1), data=p_labels, parse=TRUE, inherit.aes=F, size = 4)
p8 = p7 + theme(
  strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
  axis.title=element_text(size=14)) + ylim(-0.4, 0.6)
Figure2J = p8 + geom_text(x=2, y=-0.35, aes(label=label2), data=p_labels, parse=F, inherit.aes=F, size = 4)

#Read in KIR trans eQTL data
kir <- read.table("Figure2L.data.txt", header = T)

#construct dataframe for facet plot
for_graph <- as.data.frame(rbind(cbind(kir$geno, kir$resid.TAB1, "TAB1"), 
                                 cbind(kir$geno, kir$resid.GPRASP1, "GPRASP1"), 
                                 cbind(kir$geno, kir$resid.LRMP, "LRMP"), 
                                 cbind(kir$geno, kir$resid.HOXD4, "HOXD4"), 
                                 cbind(kir$geno, kir$resid.PSD, "PSD"), 
                                 cbind(kir$geno, kir$resid.HMOX2, "HMOX2"), 
                                 cbind(kir$geno, kir$resid.AC128709.4, "AC128709.4"), 
                                 cbind(kir$geno, kir$resid.THADA, "THADA"), 
                                 cbind(kir$geno, kir$resid.MT1L, "MT1L"),
                                 cbind(kir$geno, kir$resid.LINC01235, "LINC01235")))

colnames(for_graph) <- c("geno", "expn", "gene")
for_graph$geno <- as.numeric(as.character(for_graph$geno))
for_graph$expn <- as.numeric(as.character(for_graph$expn))

p_labels = data.frame(gene = c("TAB1", 
                               "GPRASP1", 
                               "LRMP", 
                               "HOXD4",
                               "PSD",
                               "HMOX2",
                               "AC128709.4",
                               "THADA",
                               "MT1L",
                               "LINC01235"), 
                      label1 = c("italic(P)==2.06%*%10^-13", 
                                 "italic(P)==3.58%*%10^-24", 
                                 "italic(P)==4.17%*%10^-16", 
                                 "italic(P)==5.07%*%10^-19", 
                                 "italic(P)==4.64%*%10^-26", 
                                 "italic(P)==5.80%*%10^-24", 
                                 "italic(P)==8.93%*%10^-21", 
                                 "italic(P)==4.48%*%10^-28", 
                                 "italic(P)==3.99%*%10^-29",
                                 "italic(P)==4.94%*%10^-13"))

#Panel 2L: Effect of KIR2DS4del genotype on KIR regulatory network genes in trans.
p1 = ggplot(for_graph, aes(x=factor(round(geno,0)), y=expn))
p3 = p1 + geom_dotplot(binaxis="y", binwidth=0.007, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2)
p4 = p3 + scale_x_discrete(breaks=c(0,1,2), labels = c("0", "1", "2")) + aes(fill = factor(round(geno,0)), col=factor(round(geno,0))) + scale_fill_manual(values = cols[rep(4,9)]) + scale_colour_manual(values = cols[rep(4,9)])
p5 = p4 + ylab("Trans expression") + xlab("KIR2DS4del")
p6 = p5 + theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15))
p7 = p6 + facet_wrap(vars(gene), ncol = 10) + geom_text(x=2, y=0.575, aes(label=label1), data=p_labels, parse=TRUE, inherit.aes=F, size = 4)
p8 = p7 + theme(
  strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold.italic"), legend.position = "none", axis.text=element_text(size=12),
  axis.title=element_text(size=14)) + ylim(-0.4, 0.6)
Figure2L = p8

#Create data for circos plot for KIR trans network
bed2 <- as.data.frame(cbind(c("chr22",
                              "chrX",
                              "chr12",
                              "chr2",
                              "chr10",
                              "chr16",
                              "chr3",
                              "chr2",
                              "chr16",
                              "chr9"), 
                            c(39824061,
                              101912972,
                              25157456,
                              177017741,
                              104162476,
                              4560155,
                              197184502,
                              43657396,
                              56652438,
                              13406379), 
                            c(39824062,
                              101912973,
                              25157457,
                              177017742,
                              104162477,
                              4560156,
                              197184503,
                              43657397,
                              56652439,
                              13406380)
))


colnames(bed2) <- c("chr", "start", "end")
bed2$chr <- as.character(bed2$chr)
bed2$start <- as.numeric(as.character(bed2$start))
bed2$end <- as.numeric(as.character(bed2$end))

bed2$chr.b <- c(22,
                "X",
                12,
                2,
                10,
                16,
                3,
                2,
                16,
                9)

bed3 <- bed2[order(bed2$chr.b, bed2$start),]
eQTL.gr <- with(bed3,GRanges(chr.b, IRanges(start, end)))

data("hg19Ideogram", package = "biovizBase")
seqlevelsStyle(hg19Ideogram) <- "NCBI"
seqlevels(hg19Ideogram, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
gr2 <- GRanges(c("19"), IRanges(55344131, width = 1))
seqlevels(gr2, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(gr2) <- seqlengths(hg19Ideogram)
values(eQTL.gr)$to.gr <- gr2
seqlevels(eQTL.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(eQTL.gr) <- seqlengths(hg19Ideogram)


#Panel 2K: Circos plot depicting a trans regulatory network of 10 genes, mediated mediated by KIR2DS4del.
Figure2K <- ggplot() + layout_circle(hg19Ideogram, geom = "ideo", fill = "gray70",
                                       radius = 30, trackWidth = 4)
Figure2K <- Figure2K + layout_circle(hg19Ideogram, geom = "scale", size = 2, radius = 35, trackWidth= 2)
Figure2K <- Figure2K + layout_circle(hg19Ideogram, geom = "text", aes(label = seqnames), vjust = 0, radius = 38, trackWidth = 7)
Figure2K <- Figure2K + layout_circle(eQTL.gr, geom = "link", linked.to = "to.gr", color = cols[4], radius=27, trackwidth=6, size=2)
Figure2K <- Figure2K + annotate("text", x=0, y=0, label = "KIR2DS4del", size = 10)

#Create data for circos plot for GNLY trans network
bed2 <- as.data.frame(cbind(c("chr1",
                              "chr4",
                              "chr13",
                              "chr6",
                              "chr20"), 
                            c(161016901,
                              40156212,
                              78474714,
                              56919561,
                              2637727), 
                            c(161016902,
                              40156213,
                              78474715,
                              56919562,
                              2637728)
))


colnames(bed2) <- c("chr", "start", "end")
bed2$chr <- as.character(bed2$chr)
bed2$start <- as.numeric(as.character(bed2$start))
bed2$end <- as.numeric(as.character(bed2$end))
bed2$chr.b <- c(1,
                4,
                13,
                6,
                20)
bed3 <- bed2[order(bed2$chr.b, bed2$start),]
eQTL.gr <- with(bed3,GRanges(chr.b, IRanges(start, end)))
data("hg19Ideogram", package = "biovizBase")
seqlevelsStyle(hg19Ideogram) <- "NCBI"
seqlevels(hg19Ideogram, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
gr2 <- GRanges(c("2"), IRanges(85921353, width = 1))
seqlevels(gr2, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(gr2) <- seqlengths(hg19Ideogram)
values(eQTL.gr)$to.gr <- gr2
seqlevels(eQTL.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(eQTL.gr) <- seqlengths(hg19Ideogram)

#Panel 2C: Circos plot depicting a trans regulatory network of 5 genes, mediated in cis by GNLY expression.
Figure2C <- ggplot() + layout_circle(hg19Ideogram, geom = "ideo", fill = "gray70",
                                       radius = 30, trackWidth = 4)
Figure2C <- Figure2C+ layout_circle(hg19Ideogram, geom = "scale", size = 2, radius = 35, trackWidth= 2)
Figure2C <- Figure2C + layout_circle(hg19Ideogram, geom = "text", aes(label = seqnames), vjust = 0, radius = 38, trackWidth = 7)
Figure2C <- Figure2C + layout_circle(eQTL.gr, geom = "link", linked.to = "to.gr", color = cols[3], radius=27, trackwidth=6, size=2)
Figure2C <- Figure2C + annotate("text", x=0, y=0, label = "GNLY", size = 10)

#Create data for circos plot for MC1R trans network
bed2 <- as.data.frame(cbind(c("chr15",
                              "chr4",
                              "chr3",
                              "chr1",
                              "chr4"), 
                            c(98980540,
                              129783182,
                              113847724,
                              31441012,
                              96761239), 
                            c(98980541,
                              129783183,
                              113847725,
                              31441013,
                              96761240)
))
colnames(bed2) <- c("chr", "start", "end")
bed2$chr <- as.character(bed2$chr)
bed2$start <- as.numeric(as.character(bed2$start))
bed2$end <- as.numeric(as.character(bed2$end))
bed2$chr.b <- c(15,
                4,
                3,
                1,4)
bed3 <- bed2[order(bed2$chr.b, bed2$start),]
eQTL.gr <- with(bed3,GRanges(chr.b, IRanges(start, end)))
data("hg19Ideogram", package = "biovizBase")
seqlevelsStyle(hg19Ideogram) <- "NCBI"
seqlevels(hg19Ideogram, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
gr2 <- GRanges(c("16"), IRanges(89974605, width = 1))
seqlevels(gr2, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(gr2) <- seqlengths(hg19Ideogram)
values(eQTL.gr)$to.gr <- gr2
seqlevels(eQTL.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(eQTL.gr) <- seqlengths(hg19Ideogram)

#Panel 2F: Circos plot depecting a trans regulatory network of 4 genes, mediated in cis by MC1R expression.
Figure2F <- ggplot() + layout_circle(hg19Ideogram, geom = "ideo", fill = "gray70",
                                        radius = 30, trackWidth = 4)
Figure2F <- Figure2F+ layout_circle(hg19Ideogram, geom = "scale", size = 2, radius = 35, trackWidth= 2)
Figure2F <- Figure2F + layout_circle(hg19Ideogram, geom = "text", aes(label = seqnames), vjust = 0, radius = 38, trackWidth = 7)
Figure2F <- Figure2F + layout_circle(eQTL.gr, geom = "link", linked.to = "to.gr", color = cols[2], radius=27, trackwidth=6, size=2)
Figure2F <- Figure2F + annotate("text", x=0, y=0, label = "MC1R", size = 10)

#Create data for circos plot for UVSSA trans network
bed2 <- as.data.frame(cbind(c("chr19",
                              "chr20",
                              "chr21",
                              "chr19"), 
                            c(40972148,
                              3063202,
                              38071433,
                              44080952), 
                            c(40972149,
                              3063203,
                              38071434,
                              44080953)
))
colnames(bed2) <- c("chr", "start", "end")
bed2$chr <- as.character(bed2$chr)
bed2$start <- as.numeric(as.character(bed2$start))
bed2$end <- as.numeric(as.character(bed2$end))
bed2$chr.b <- c(19,
                20,
                21,
                19)
bed3 <- bed2[order(bed2$chr.b, bed2$start),]
eQTL.gr <- with(bed3,GRanges(chr.b, IRanges(start, end)))
data("hg19Ideogram", package = "biovizBase")
seqlevelsStyle(hg19Ideogram) <- "NCBI"
seqlevels(hg19Ideogram, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
gr2 <- GRanges(c("4"), IRanges(1012300, width = 1))
seqlevels(gr2, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(gr2) <- seqlengths(hg19Ideogram)
values(eQTL.gr)$to.gr <- gr2
seqlevels(eQTL.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(eQTL.gr) <- seqlengths(hg19Ideogram)

#Panel 2I: Circos plot depicting a trans regulatory network of 4 genes, mediated in cis by UVSSA expression.
Figure2I <- ggplot() + layout_circle(hg19Ideogram, geom = "ideo", fill = "gray70",
                                        radius = 30, trackWidth = 4)
Figure2I <- Figure2I+ layout_circle(hg19Ideogram, geom = "scale", size = 2, radius = 35, trackWidth= 2)
Figure2I <- Figure2I + layout_circle(hg19Ideogram, geom = "text", aes(label = seqnames), vjust = 0, radius = 38, trackWidth = 7)
Figure2I <- Figure2I + layout_circle(eQTL.gr, geom = "link", linked.to = "to.gr", color = cols[6], radius=27, trackwidth=6, size=2)
Figure2I <- Figure2I + annotate("text", x=0, y=0, label = "UVSSA", size = 10)

#construct whole figure
Figure2 <- (((Figure2A))/((Figure2B|Figure2C|Figure2D))/((Figure2E|Figure2F|Figure2G))/((Figure2H|Figure2I|Figure2J))/((Figure2K|Figure2L))) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
