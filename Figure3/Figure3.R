library(ggplot2)
library(patchwork)
library(ggrepel)
library(dplyr)
library(RColorBrewer)

#read in TWAS summary statisitics
uc <- read.table("UC_full.txt", header = T)
ra <- read.table("RA_full.txt", header = T)
pbc <- read.table("PBC_full.txt", header = T)
crohns <- read.table("CD_full.txt", header = T)
sle <- read.table("SLE_full.txt", header = T)

#data for Crohn's disease manahattan plot
twas <- crohns
twas$gene <- as.character(twas$gene)
twas$fdr <- p.adjust(twas$TWAS.P, method = "fdr")

#highlight hits of interest
twas$is.sig <- 0
twas$is.sig[which(twas$fdr<0.05 & twas$cond.sig==1)] <- 1
twas$novel.cond.sig <- 0
twas$novel.cond.sig[which(twas$cond.sig==1 & twas$novel==1)] <- 1
twas$nk.uniq <- 0
twas$nk.uniq[which(twas$twas.fdr.sig==1 & twas$unique==1 & twas$cond.sig==1 & twas$ppa<0.2 & twas$COLOC.PP4>0.8)] <- 1

don <- twas %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(P1)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(twas, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, P1) %>%
  mutate( BPcum=P1+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Plot Crohn's disease TWAS Manhattan
manh <- ggplot(don, aes(x=BPcum, y=-log10(TWAS.P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=2.5) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 11 )) +
  
  # Add highlighted points
  geom_point(data=subset(don, novel.cond.sig==1), color=cols2[6], size=4) +

  # Add labls
  geom_text_repel( data=subset(don, is.sig==1), aes(label=gene), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +

  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +
  xlab("chromosome") +
  ylab("-log10(p)") +
  ylim(NA,max(-log10(na.omit(don$TWAS.P)))+1) +
  ggtitle("Crohns disease") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    plot.title=element_text(size=25, face = "bold")
  )

cd.manh <- manh

#data for Primary Biliary Cirrhosis manahattan plot
twas <- pbc
twas$gene <- as.character(twas$gene)
twas$fdr <- p.adjust(twas$TWAS.P, method = "fdr")
#highlight hits of interest
twas$is.sig <- 0
twas$is.sig[which(twas$fdr<0.05 & twas$cond.sig==1)] <- 1
twas$novel.cond.sig <- 0
twas$novel.cond.sig[which(twas$cond.sig==1 & twas$novel==1)] <- 1
twas$nk.uniq <- 0
twas$nk.uniq[which(twas$twas.fdr.sig==1 & twas$unique==1 & twas$cond.sig==1 & twas$ppa<0.2 & twas$COLOC.PP4>0.8)] <- 1

don <- twas %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(P1)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(twas, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, P1) %>%
  mutate( BPcum=P1+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Plot PBC TWAS Manhattan
manh <- ggplot(don, aes(x=BPcum, y=-log10(TWAS.P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=2.5) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 11 )) +
  
  # Add highlighted points
  geom_point(data=subset(don, nk.uniq==1), color=cols2[4], size=7) +
  geom_point(data=subset(don, novel.cond.sig==1), color=cols2[6], size=4) +
  
  # Add labels
  geom_text_repel( data=subset(don, is.sig==1), aes(label=gene), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +
  
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("chromosome") +
  ylab("-log10(p)") +
  ylim(NA,max(-log10(na.omit(don$TWAS.P)))+1) +
  ggtitle("Primary biliary cirrhosis") +

  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    plot.title=element_text(size=25, face = "bold")
  )

pbc.manh <- manh

#data for SLE manahattan plot
twas <- sle
twas$gene <- as.character(twas$gene)
twas$fdr <- p.adjust(twas$TWAS.P, method = "fdr")

#highlights hits of interest
twas$is.sig <- 0
twas$is.sig[which(twas$fdr<0.05 & twas$cond.sig==1)] <- 1
twas$novel.cond.sig <- 0
twas$novel.cond.sig[which(twas$cond.sig==1 & twas$novel==1)] <- 1
twas$nk.uniq <- 0
twas$nk.uniq[which(twas$twas.fdr.sig==1 & twas$unique==1 & twas$cond.sig==1 & twas$ppa<0.2 & twas$COLOC.PP4>0.8)] <- 1

don <- twas %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(P1)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(twas, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, P1) %>%
  mutate( BPcum=P1+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Plot SLE TWAS Manhattan
manh <- ggplot(don, aes(x=BPcum, y=-log10(TWAS.P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=2.5) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 11 )) +
  
  # Add highlighted points
  geom_point(data=subset(don, nk.uniq==1), color=cols2[4], size=7) +
  geom_point(data=subset(don, novel.cond.sig==1), color=cols2[6], size=4) +
  
  # Add labels
  geom_text_repel( data=subset(don, is.sig==1), aes(label=gene), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +

  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("chromosome") +
  ylab("-log10(p)") +
  ylim(NA,max(-log10(na.omit(don$TWAS.P)))+1) +
  ggtitle("Systemic lupus erythematosus") +

  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    plot.title=element_text(size=25, face = "bold")
  )

sle.manh <- manh

#data for rheumatoid arthritis manahattan plot
twas <- ra
twas$gene <- as.character(twas$gene)
twas$fdr <- p.adjust(twas$TWAS.P, method = "fdr")

#highlight hits of interest
twas$is.sig <- 0
twas$is.sig[which(twas$fdr<0.05 & twas$cond.sig==1)] <- 1
twas$novel.cond.sig <- 0
twas$novel.cond.sig[which(twas$cond.sig==1 & twas$novel==1)] <- 1
twas$nk.uniq <- 0
twas$nk.uniq[which(twas$twas.fdr.sig==1 & twas$unique==1 & twas$cond.sig==1 & twas$ppa<0.2 & twas$COLOC.PP4>0.8)] <- 1

don <- twas %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(P1)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(twas, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, P1) %>%
  mutate( BPcum=P1+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Plot rheumatoid arthritis TWAS Manhattan
manh <- ggplot(don, aes(x=BPcum, y=-log10(TWAS.P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=2.5) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 11 )) +
  
  # Add highlighted points
  geom_point(data=subset(don, novel.cond.sig==1), color=cols2[6], size=4) +
  #geom_point(data=subset(don, nk.uniq==1), color=cols2[4], size=5) +
  
  # Add labels
  geom_text_repel( data=subset(don, is.sig==1), aes(label=gene), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +
  
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("chromosome") +
  ylab("-log10(p)") +
  ylim(NA,max(-log10(na.omit(don$TWAS.P)))+1) +
  ggtitle("Rheumatoid arthritis") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    plot.title=element_text(size=25, face = "bold")
  )
ra.manh <- manh

#data for ulcerative colitis manahattan plot
twas <- uc
twas$gene <- as.character(twas$gene)
twas$fdr <- p.adjust(twas$TWAS.P, method = "fdr")

#highlight hits of interest
twas$is.sig <- 0
twas$is.sig[which(twas$fdr<0.05 & twas$cond.sig==1)] <- 1
twas$novel.cond.sig <- 0
twas$novel.cond.sig[which(twas$cond.sig==1 & twas$novel==1)] <- 1
twas$nk.uniq <- 0
twas$nk.uniq[which(twas$twas.fdr.sig==1 & twas$unique==1 & twas$cond.sig==1 & twas$ppa<0.2 & twas$COLOC.PP4>0.8)] <- 1

don <- twas %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(P1)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(twas, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, P1) %>%
  mutate( BPcum=P1+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")

#Plot ulcerative colitis TWAS Manhattan
manh <- ggplot(don, aes(x=BPcum, y=-log10(TWAS.P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=2.5) +
  scale_color_manual(values = rep(c(cols[8], cols2[8]), 11 )) +
  
  # Add highlighted points
    geom_point(data=subset(don, novel.cond.sig==1), color=cols2[6], size=4) +

  
  # Add labels
  geom_text_repel( data=subset(don, is.sig==1), aes(label=gene), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 1) +
  
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("chromosome") +
  ylab("-log10(p)") +
  ylim(NA,max(-log10(na.omit(don$TWAS.P)))+1) +
  ggtitle("Ulcerative colitis") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=16),
    plot.title=element_text(size=25, face = "bold")
  )
uc.manh <- manh

Figure3 <- (((uc.manh))/((ra.manh))/((sle.manh))/((pbc.manh))/((cd.manh)))

