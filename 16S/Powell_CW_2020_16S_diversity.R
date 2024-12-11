#Maya Powell
#March 28th - April 3, 2023
#August 17th 2023, making pub quality figures
#August 13th 2024 FINISHING GRAPHS AND STATS
#16S DIVERSITY analysis
#Based on scripts from Anastasia Dulskiy and Nicola Kreifall and Hannah Aichelman
setwd("~/Documents/Castillo Lab/CW_2020/CW_2020_16S")
# Setup

## Packages
library(ggplot2)
library(cowplot)
#BiocManager::install("phyloseq")
library(phyloseq)
#BiocManager::install("car")
library(car)
library(ggpubr)
library(vegan) #vegan 2.6-4
#BiocManager::install("dada2")
library(dada2)
library(RColorBrewer)
library(wesanderson)
library(dplyr)
#install.packages("lme4")
library("lme4")
#install.packages("glmmTMB")
library("glmmTMB")

##COMING BACK TO THIS ANALYSIS???
##SEE LINE 80 TO READ IN DIVERSITY DATAFRAME!!!!!! WOOOHOOO!!

## Read in data
#for this analysis - should be done on raw, untrimmed dataset
#therefore, using ps.less! (raw, untrimmed data with just 9 samples with less than 1000 read counts removed)

#load in Phyloseq objects
#ps.cleanest <- readRDS("CW_2020_16S_ps.cleanest.RDS")
#ps.cleanest #27653 taxa and 144 samples, RAW DATA untrimmed, unrarefied
#remove duplicates from lane1 - they are very close (not stat sig dif) and lane 2 has higher coverage so keep those
#ps.cleanest.nd = subset_samples(ps.cleanest, id!="N10" & id!="N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.cleanest.nd <- data.frame(sample_data(ps.cleanest.nd))

ps.less <- readRDS("CW_2020_16S_ps.less.RDS")
ps.less #27653 taxa and 135 samples, samples with counts <1,000 removed (9 samples removed)
#in addition, removing the astreoides samples
ps.less.nd = subset_samples(ps.less, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
ps.less.nd.no.astreoides = subset_samples(ps.less.nd, id!="J9" & id!="J10" & id!="J11" & id!="J12" & id!="K1" & id!="K2")
samdf.less.nd <- data.frame(sample_data(ps.less.nd.no.astreoides))

#ps.rare <- readRDS("CW_2020_16S_ps.rare.RDS")
#ps.rare #23759 taxa and 98 samples, rarefied to 10,000 (46 samples removed)
#ps.rare.nd = subset_samples(ps.rare, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.rare.nd <- data.frame(sample_data(ps.rare.nd))

ps.trim <- readRDS("CW_2020_16S_ps.trim.RDS")
ps.trim #835 taxa and 134 samples, trimmed witN MCMC OTU for low abundance taxa (10 samples removed)
ps.trim.nd = subset_samples(ps.trim, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
ps.trim.nd.no.astreoides = subset_samples(ps.trim.nd, id!="J9" & id!="J10" & id!="J11" & id!="J12" & id!="K1" & id!="K2")
samdf.trim.nd <- data.frame(sample_data(ps.trim.nd.no.astreoides))

#ps.trim.rare <- readRDS("CW_2020_16S_ps.trim.rare.RDS")
#ps.trim.rare #835 taxa and 120 samples, trimmed dataset rarefied to 3498 (24 samples removed)
#ps.trim.rare.nd = subset_samples(ps.trim.rare, id!="H10" & id!="H11" & id!="H2" & id!="H3" & id!="H4" & id!="H5" & id!="H6" & id!="H7" & id!="H8" & id!="H9")
#samdf.trim.rare.nd <- data.frame(sample_data(ps.trim.rare.nd))

# Diversity
#[Notes from phyloseq author](https://rdrr.io/bioc/phyloseq/man/estimate_richness.html)
#Visualize alpha-diversity - Should be done on raw, untrimmed dataset

#Generate diversity metrics
#run through this code to estimate these based on each different dataset above
df <- data.frame(estimate_richness(ps.less.nd.no.astreoides, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df$id <- rownames(df)
df.div <- merge(df,samdf.less.nd,by="id") #add sample data
#df.div <- df.div %>% select(-collection_latitude,-collection_longitude,-collection_year,-collection_month,-collection_depth,-collection_date)

#add evenness = shannon diversity divided by species richness
df.div$even <- df.div$Shannon/(log(df.div$Observed))

df.div$r_b_m_y <- paste(df.div$reef_bay,df.div$m_y)

#SAVE DATAFRAME
write.csv(df.div, file = "CW_2020_diversity_df.ps.less.csv", row.names = FALSE)

#READ BACK IN HERE!! START HERE IF COMING BACK TO THIS SCRIPT!!
df.div <- read.csv("CW_2020_diversity_df.ps.less.csv")

#melt to get all variables together
df.div.melt <- reshape2::melt(df.div)
library(utils)
df.div.melt$site_zone <- as.factor(df.div.melt$site_zone)
lapply(df.div.melt, levels)

#Separate by species
df.div.ss <- subset(df.div,host_species=="siderea")
df.div.sr <- subset(df.div,host_species=="radians")
df.div.pp <- subset(df.div,host_species=="porites")
#df.div.pa <- subset(df.div,host_species=="astreoides")

#Separate by site
#df.div.SWB <- subset(df.div,site_zone=="SW_bay")
#df.div.SWR <- subset(df.div,site_zone=="SW_reef")
#df.div.SMB <- subset(df.div,site_zone=="SM_bay")
#df.div.SMR <- subset(df.div,site_zone=="SM_reef")

###ALPHA DIVERSITY METRICS###
#GRAPHS ONLY - stats down below
#for graphs, using these colors
#Santa Martha Bay
"#E78AC3"
"#F9AD2A"
"coral3"
#Santa Martha Reef
"#8DA0CB"
"paleturquoise"
#Spaanse Water Bay
"#FC8D62"
"lightsalmon1"
#Spaanse Water Reef
"#66C2A5"
"deepskyblue4"
#scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+

### Shannon
#shannon by site
###Shannon by site & timepoint
#all species together
gg.site.sha <- ggplot(df.div,aes(x=site_zone,y=Shannon,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5)+
  theme_classic(base_size = 10)+
  scale_fill_brewer(palette = "Set2", direction = - 1)+
  ylab("Shannon index")+
  xlab("Site")+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())+
  facet_wrap(m_y~host_species,scales = "free_x")
gg.site.sha
#ggsave(gg.site.sha,file="16S.all.sha.pdf",h=20,w=10)

#Shannon by species across everything
gg.spp.sha <- ggplot(df.div,aes(x=host_species,y=Shannon,fill=host_species))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size=22)+
  geom_jitter(alpha=0.5)+
  scale_fill_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))+
  ylab("Shannon Index")+
  xlab("Coral Species")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.spp.sha
ggsave(gg.spp.sha,file="16S.species.sha.pdf",h=15,w=10)

##Separated by species

#siderastrea siderea SHANNON PLOT
gg.site.sha.ss <- ggplot(df.div.ss,aes(x=site_zone,y=Shannon,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 25)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun = max, vjust = -0.5, size=8)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"), labels = c("Santa Martha Bay","Santa Martha Reef","Spaanse Water Bay", "Spaanse Water Reef"))+
  ylab("Shannon index")+
  xlab("Site")+
  ylim(0, 7)+
  #ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  #facet_wrap(m_y~.,scales = "free_x")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
gg.site.sha.ss
#ggsave(gg.site.sha.ss,file="16S.ss.sha.pdf",h=6,w=8)

#siderastrea radians shannon plot
gg.site.sha.sr <- ggplot(df.div.sr,aes(x=site_zone,y=Shannon,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3","lightsalmon1"))+
  ylab("Shannon index")+
  xlab("Site")+
  ylim(0, 7)+
  theme_classic(base_size = 25)+
  #ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.sha.sr
#ggsave(gg.site.sha.sr,file="16S.sr.sha.pdf",h=8,w=6)

#porites porites shannon graph
gg.site.sha.pp <- ggplot(df.div.pp,aes(x=site_zone,y=Shannon,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3", "paleturquoise"))+
  ylab("Shannon index")+
  xlab("Site")+
  ylim(0, 7)+
  theme_classic(base_size = 25)+
  #ggtitle(substitute(paste(italic("Porites porites"))))+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.sha.pp
#ggsave(gg.site.sha.pp,file="16S.pp.sha.pdf",h=8,w=6)

#put all 3 figures together
#Shannon - all species
#srppsha <- ggarrange(gg.site.sha.sr,gg.site.sha.pp,nrow = 1)
#all_shannon_plot <- ggarrange(gg.site.sha.ss, srppsha,nrow = 1)
#ggsave(all_shannon_plot, file = "all_shannon_plot.ps.less.pdf", h=8, w=16)

#gg.site.sha.pa <- ggplot(df.div.pa,aes(x=site_zone,y=Shannon,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_color_manual(values=c("#FC8D62","#66C2A5"))+
#   stat_summary(geom = 'text', label = c("a","b"), fun.y = max, vjust = -0.5, size =6)+
#   ylab("Shannon index")+
#   xlab("Site")+
#   theme_classic(base_size = 22)+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   #facet_wrap(m_y~.,scales = "free_x")
# gg.site.sha.pa
# #ggsave(gg.site.sha.pa,file="16S.pa.sha.pdf",h=8,w=6)

### INVERSE SIMPSON INDEX
#all species together by site
# gg.site.simp <- ggplot(df.div,aes(x=site_zone,y=InvSimpson,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   theme_classic(base_size=10)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("InvSimpson index")+
#   xlab("Site")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(m_y~host_species,scales = "free_x")
# gg.site.simp
# #ggsave(gg.site.simp,file="16S.all.simp.pdf",h=20,w=10)

#all sites together by species
# gg.spp.simp <- ggplot(df.div,aes(x=host_species,y=InvSimpson,color=host_species))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("InvSimpson index")+
#   xlab("Species")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(site_zone~.,scales = "free_x")
# gg.spp.simp
#ggsave(gg.spp.simp,file="16S.species.simp.pdf",h=15,w=10)

##Plots separated by species - INVERSE SIMPSON

#species comparison - simposon
gg.spp.simp <- ggplot(df.div,aes(x=host_species,y=InvSimpson,fill=host_species))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size=22)+
  geom_jitter(alpha=0.5)+
  scale_fill_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))+
  ylab("Inverse Simpson's Index")+
  xlab("Coral Species")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
#facet_wrap(site_zone~.,scales = "free_x")
gg.spp.simp

#siderastrea siderea inverse simpson plot
gg.site.simp.ss <- ggplot(df.div.ss,aes(x=site_zone,y=InvSimpson,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 25)+
  #stat_summary(geom = 'text', label = c("ab","a","a","b"), fun.y = max, vjust = -0.5, size=8)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"), labels = c("Santa Martha Bay","Santa Martha Reef","Spaanse Water Bay", "Spaanse Water Reef"))+
  ylab("Inverse Simpson Index")+
  xlab("Site")+
  ylim(0, 290)+
  ggtitle(substitute(paste(italic(" "))))+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.simp.ss
#ggsave(gg.site.simp.ss,file="16S.ss.simp.pdf",h=6,w=8)

#siderastrea radians inverse simpson plot
gg.site.simp.sr <- ggplot(df.div.sr,aes(x=site_zone,y=InvSimpson,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3","lightsalmon1"))+
  ylab("Inverse Simpson Index")+
  xlab("Site")+
  ylim(0, 290)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.simp.sr
#ggsave(gg.site.simp.sr,file="16S.sr.simp.pdf",h=8,w=6)

#porites porites inverse simpson plot
gg.site.simp.pp <- ggplot(df.div.pp,aes(x=site_zone,y=InvSimpson,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3", "paleturquoise"))+
  ylab("Inverse Simpson Index")+
  xlab("Site")+
  ylim(0, 290)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.simp.pp
#ggsave(gg.site.simp.pp,file="16S.pp.simp.pdf",h=8,w=6)

#put plots together
#srppsimp <- ggarrange(gg.site.simp.sr,gg.site.simp.pp,nrow=1)
#all_simpson_plot <- ggarrange(gg.site.simp.ss, srppsimp, nrow=1)
#ggsave(all_simpson_plot, file = "all_simpson_plot.ps.less.pdf",h=8,w=16)

# gg.site.simp.pa <- ggplot(df.div.pa,aes(x=site_zone,y=InvSimpson,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_color_manual(values=c("#FC8D62","#66C2A5"))+
#   #stat_summary(geom = 'text', label = c("a","b"), fun.y = max, vjust = -0.5, size = 6)+
#   ylab("InvSimpson index")+
#   xlab("Site")+
#   theme_classic(base_size = 22)+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# #facet_wrap(m_y~.,scales = "free_x")
# gg.site.simp.pa
#ggsave(gg.site.simp.pa,file="16S.pa.simp.pdf",h=8,w=6)

### Richness ###
# #all species together - plot by site
# gg.site.rich <- ggplot(df.div,aes(x=site_zone,y=Observed,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   theme_classic(base_size=10)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("Richness")+
#   xlab("Site")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(m_y~host_species,scales = "free_x")
# gg.site.rich
# #ggsave(gg.site.rich,file="16S.all.rich.pdf",h=20,w=10)
# 
# #all sites toether - plot by species
# gg.spp.rich <- ggplot(df.div,aes(x=host_species,y=Observed,color=host_species))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("Richness")+
#   xlab("Species")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(site_zone~.,scales = "free_x")
# gg.spp.rich
#ggsave(gg.spp.rich,file="16S.species.rich.pdf",h=15,w=10)

##RICHNESS GRAPHS/PLOTS BY SPECIES

#species comparison
gg.spp.rich <- ggplot(df.div,aes(x=host_species,y=Observed,fill=host_species))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size=22)+
  geom_jitter(alpha=0.5)+
  scale_fill_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))+
  ylab("Richness")+
  xlab("Coral Species")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
#facet_wrap(site_zone~.,scales = "free_x")
gg.spp.rich

#siderastrea siderea richness
gg.site.rich.ss <- ggplot(df.div.ss,aes(x=site_zone,y=Observed,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 25)+
  #stat_summary(geom = 'text', label = c("ab","a","a","b"), fun.y = max, vjust = -1, size=8)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"), labels = c("Santa Martha Bay","Santa Martha Reef","Spaanse Water Bay", "Spaanse Water Reef"))+
  ylab("Richness")+
  xlab("Site")+
  ylim(0, 2300)+
  ggtitle(substitute(paste(italic(" "))))+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.rich.ss
ggsave(gg.site.rich.ss,file="16S.ss.rich.pdf",h=6,w=8)

#siderastrea radians richness
gg.site.rich.sr <- ggplot(df.div.sr,aes(x=site_zone,y=Observed,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  #stat_summary(geom = 'text', label = c("a","a"), fun.y = max, size=6)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3","lightsalmon1"))+
  ylab("Richness")+
  xlab("Site")+
  ylim(0, 2250)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.rich.sr
#ggsave(gg.site.rich.sr,file="16S.sr.rich.pdf",h=8,w=6)

#porites porites richness
gg.site.rich.pp <- ggplot(df.div.pp,aes(x=site_zone,y=Observed,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  #stat_summary(geom = 'text', label = c("a","a"), fun.y = max, vjust = -1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3", "paleturquoise"))+
  ylab("Richness")+
  xlab("Site")+
  ylim(0, 2250)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.rich.pp
#ggsave(gg.site.rich.pp,file="16S.pp.rich.pdf",h=8,w=6)

#put plots together
#srpprich <- ggarrange(gg.site.rich.sr,gg.site.rich.pp,nrow=1)
#all_richness_plot <- ggarrange(gg.site.rich.ss, srpprich, nrow=1)
#ggsave(all_richness_plot, file = "all_richness_plot.ps.less.pdf",h=8,w=16)

#porites astreoides richness
# gg.site.rich.pa <- ggplot(df.div.pa,aes(x=site_zone,y=Observed,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_color_manual(values=c("#FC8D62","#66C2A5"))+
#   stat_summary(geom = 'text', label = c("a","b"), fun.y = max, vjust = -0.5, size =6)+
#   ylab("Richness")+
#   xlab("Site")+
#   theme_classic(base_size = 22)+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# #facet_wrap(m_y~.,scales = "free_x")
# gg.site.rich.pa
#ggsave(gg.site.rich.pa,file="16S.pa.rich.pdf",h=8,w=6)

### Evenness
#all species together - site comparisons evenness
# gg.site.even <- ggplot(df.div,aes(x=site_zone,y=even,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("Evenness")+
#   xlab("Site")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(m_y~host_species,scales = "free_x")
# gg.site.even
#ggsave(gg.site.even,file="16S.all.even.pdf",h=20,w=10)

#all sites together - species comparisons evenness
# gg.spp.even <- ggplot(df.div,aes(x=host_species,y=even,color=host_species))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_colour_brewer(palette = "Set2", direction = - 1)+
#   ylab("Evenness")+
#   xlab("Species")+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_wrap(site_zone~.,scales = "free_x")
# gg.spp.even
#ggsave(gg.spp.even,file="16S.species.even.pdf",h=15,w=10)

#species comparison
gg.spp.even <- ggplot(df.div,aes(x=host_species,y=even,fill=host_species))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size=22)+
  geom_jitter(alpha=0.5)+
  scale_fill_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))+
  ylab("Evenness")+
  xlab("Coral Species")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
#facet_wrap(site_zone~.,scales = "free_x")
gg.spp.even


#siderastrea siderea evenness plot
gg.site.even.ss <- ggplot(df.div.ss,aes(x=site_zone,y=even,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 25)+
  #stat_summary(geom = 'text', label = c("ab","a","a","b"), fun.y = max, vjust = -0.5,size=8)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  ylab("Evenness")+
  xlab("Site")+
  ylim(0, 1)+
  ggtitle(substitute(paste(italic(" "))))+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.even.ss
#ggsave(gg.site.even.ss,file="16S.ss.even.pdf",h=8,w=8)

#siderastrea radians evenness plot
gg.site.even.sr <- ggplot(df.div.sr,aes(x=site_zone,y=even,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  #stat_summary(geom = 'text', label = c("a","a"), fun.y = max, size=6)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3","lightsalmon1"))+
  ylab("Evenness")+
  xlab("Site")+
  ylim(0, 1)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.even.sr
#ggsave(gg.site.even.sr,file="16S.sr.even.pdf",h=8,w=6)

#porites porites evenness plot
gg.site.even.pp <- ggplot(df.div.pp,aes(x=site_zone,y=even,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual(values = c("coral3", "paleturquoise"))+
  ylab("Evenness")+
  xlab("Site")+
  ylim(0, 1)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="none",axis.text.x = element_blank(),axis.ticks.x=element_blank())
  #facet_wrap(m_y~.,scales = "free_x")
gg.site.even.pp
#ggsave(gg.site.even.pp,file="16S.pp.even.pdf",h=8,w=6)

#put plots together
#srppeven <- ggarrange(gg.site.even.sr,gg.site.even.pp,nrow=1)
#all_evenness_plot <- ggarrange(gg.site.even.ss, srppeven, nrow=1)
#ggsave(all_evenness_plot, file = "all_evenness_plot.ps.less.pdf",h=8,w=16)

shasimprichsid <- ggarrange(gg.site.sha.ss,gg.site.simp.ss,gg.site.rich.ss, legend = FALSE, nrow=1, labels = c("A","B","C"), font.label = list(size = 30, color = "black"))
evenfaithsid <- ggarrange(gg.site.even.ss,gg.site.faith.ss, nrow=1, widths = c(0.6, 1), labels = c("D","E"), font.label = list(size = 30, color = "black"))
all_sid_alpha_plot <- ggarrange(shasimprichsid,evenfaithsid, nrow=2)
#all_sid_alpha_plot <- annotate_figure(all_sid_alpha_plot, top = text_grob("Siderastrea siderea", color = "black", hjust = 1.15, face = "italic", size = 60))
ggsave(all_sid_alpha_plot, file = "all_sid_alpha_plot.ps.less.pdf",h=13,w=17)

shasimprichrad <- ggarrange(gg.site.sha.sr,gg.site.simp.sr,gg.site.rich.sr, legend = FALSE, nrow=1, labels = c("A","B","C"), font.label = list(size = 30, color = "black"))
evenfaithrad <- ggarrange(gg.site.even.sr,gg.site.faith.sr, nrow=1, widths = c(0.5, 0.7), labels = c("D","E"), font.label = list(size = 30, color = "black"))
all_rad_alpha_plot <- ggarrange(shasimprichrad,evenfaithrad, nrow=2)
#all_rad_alpha_plot <- annotate_figure(all_rad_alpha_plot, top = text_grob("Siderastrea radians", color = "black", hjust = 1.15, face = "italic", size = 60))
ggsave(all_rad_alpha_plot, file = "all_rad_alpha_plot.ps.less.pdf",h=13,w=17)

#all_rad_alpha_plot <- ggarrange(gg.site.sha.sr, gg.site.simp.sr, gg.site.rich.sr, gg.site.even.sr, nrow=2, ncol=2)
#all_rad_alpha_plot <- annotate_figure(all_rad_alpha_plot, fig.lab = "Figure 1", fig.lab.face = "italic",fig.lab.size = 40)
                                      #top = text_grob(substitute(paste(italic("Siderastrea siderea"))), color = "black", face = "bold", size = 40))
#ggsave(all_rad_alpha_plot, file = "all_rad_alpha_plot.ps.less.pdf",h=12,w=12)

#all_por_alpha_plot <- ggarrange(gg.site.sha.pp, gg.site.simp.pp, gg.site.rich.pp, gg.site.even.pp, nrow=2, ncol=2)
#ggsave(all_por_alpha_plot, file = "all_por_alpha_plot.ps.less.pdf",h=12,w=12)

shasimprichpor <- ggarrange(gg.site.sha.pp,gg.site.simp.pp,gg.site.rich.pp, legend = FALSE, nrow=1, labels = c("A","B","C"), font.label = list(size = 30, color = "black"))
evenfaithpor <- ggarrange(gg.site.even.pp,gg.site.faith.pp, nrow=1, widths = c(0.5, 0.7), labels = c("D","E"), font.label = list(size = 30, color = "black"))
all_por_alpha_plot <- ggarrange(shasimprichpor,evenfaithpor, nrow=2)
#all_por_alpha_plot <- annotate_figure(all_por_alpha_plot, top = text_grob("Porites porites", color = "black", hjust = 1.6, face = "italic", size = 60))
ggsave(all_por_alpha_plot, file = "all_por_alpha_plot.ps.less.pdf",h=13,w=17)

#all plots by coral species together
shasimprichspp <- ggarrange(gg.spp.sha,gg.spp.simp,gg.spp.rich, nrow=1, labels = c("A","B","C"), font.label = list(size = 30, color = "black"), legend = "none")
evenfaithspp <- ggarrange(gg.spp.even,gg.spp.pd, nrow=1, widths = c(1, 1), labels = c("D","E"), font.label = list(size = 30, color = "black"), legend = "right", common.legend = TRUE)
all_spp_alpha_plot <- ggarrange(shasimprichspp,evenfaithspp, nrow=2)
#all_spp_alpha_plot <- annotate_figure(all_spp_alpha_plot, top = text_grob("Coral Species", color = "black", hjust = 1.15, face = "bold", size = 60))
ggsave(all_spp_alpha_plot, file = "all_spp_alpha_plot.pdf",h=13,w=17)

# gg.site.even.pa <- ggplot(df.div.pa,aes(x=site_zone,y=even,color=site_zone))+
#   geom_boxplot(outlier.shape=NA,size=1)+
#   geom_jitter(alpha=0.5)+
#   scale_color_manual(values=c("#FC8D62","#66C2A5"))+
#   stat_summary(geom = 'text', label = c("a","b"), fun.y = max, vjust = -0.5, size =6)+
#   ylab("Evenness")+
#   xlab("Site")+
#   theme_classic(base_size = 22)+
#   theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# #facet_wrap(m_y~.,scales = "free_x")
# gg.site.even.pa
#ggsave(gg.site.even.pa,file="16S.pa.even.pdf",h=8,w=6)

## Phylogenetic diversity (Faith's D)

#Tutorial from dada2 author: https://f1000research.com/articles/5-1492/v2

#packages for Faith's D
#install.packages('devtools')
library(devtools)
#BiocManager::install("DESeq2")
library("DESeq2")
#BiocManager::install("genefilter")
library("genefilter")
#BiocManager::install("biomformat")
library("biomformat")
BiocManager::install('twbattaglia/btools')
library(btools)

#Generate fasta file (did this in previous script, so don't need it here)
#(I'm not running the following chunk every time, because only need to generate the file once)

#rare.otu <- as.matrix(ps.rare@otu_table)
# rare.taxa <- data.frame(ps.rare@tax_table)
# rownames(rare.taxa)==colnames(rare.otu)
# colnames(rare.otu) <- rare.taxa$V8
# ids <- rownames(rare.taxa)
# path="~/Documents/curacao_2020/16S_analysis/curacao_2020.cleanest.fasta"
# uniquesToFasta(rare.otu, path, ids = ids, mode = "w", width = 20000)

#Actual analysis part: 
#(Im not running the following chunk every time, because only need to generate the files once)
#seqs <- getSequences("CW_2020_16S_less.fasta")
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#saveRDS(seqs,file="CW_2020_16S_less_phylo_seqs.rds")
#and now a trimmed one too to compare
#seqs <- getSequences("CW_2020_16S_trim.fasta")
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#saveRDS(seqs,file="CW_2020_16S_trim_phylo_seqs.rds")

#Do this next part in the cluster because it takes forever

##script faith_d.R:
##Maya Powell
#Faith's D script for getting phylogenetic diversity info
#Run this on the cluster
#July 29th 2023

library(dada2)
library(phangorn)
library(DECIPHER)

seqs <- readRDS("/proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/Faith_D/CW_2020_16S_less_phylo_seqs.rds")
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# ## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
saveRDS(fitGTR, file="/proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/Faith_D/CW_2020_16S_phylo_fitgtr_less.rds")
#make sure to put absolute filepath, otherwise it won't work
#tried to run with only 10g of memory but it ran out of memory so had to increase to 30

#make phylo.sh file
#nano phylo.sh

#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=30g
#SBATCH -n 28
#SBATCH -t 2-00:00:00

#module load r/4.3.1
#module load rstudio
#Rscript faith_d.R

##exit (cntrl X) and save!

##on cluster:
#sbatch phylo.sh
#to check status: squeue -u mayapow
#you'll also be able to see slurm output file updating as you go

##saved output as: CW_2020_16S_phylo_fitgtr_less.rds

#Now back to R

library(btools)

fitGTR <- readRDS("CW_2020_16S_phylo_fitgtr_trim.rds")

#using ps.trim here bc don't have the data for ps.less right now
#currently running on the cluster 9/16/2023

#new phyloseq object:
taxa.trim <- data.frame(ps.trim.nd.no.astreoides@tax_table)
seqtab.trim <- data.frame(ps.trim.nd.no.astreoides@otu_table)
taxa.trim$sqs <- row.names(taxa.trim) 
taxa.trim$sqs == colnames(seqtab.trim)
row.names(taxa.trim) <- taxa.trim$V8
colnames(seqtab.trim) <- taxa.trim$V8
row.names(taxa.trim) == colnames(seqtab.trim)
taxa.trim <- as.matrix(taxa.trim)

ps.trim.tree <- phyloseq(otu_table(seqtab.trim, taxa_are_rows = FALSE),
                         sample_data(samdf.trim.nd),
                         tax_table(taxa.trim),
                         phy_tree(fitGTR$tree))

pd.div <- estimate_pd(ps.trim.tree)
row.names(df.div) <- df.div$id
df.div.pd <- merge(df.div,pd.div,by=0)

## saving diversity data frame ##
#save & read back in as needed
write.csv(df.div.pd,file="df.pd.div.csv") #saving
df.div.pd <- read.csv("df.pd.div.csv",row.names=1,header=TRUE) #reading back in

##post trimming:
# fitGTR.trim <- readRDS("phylo.fitgtr.rev.cleanest.trimmed.rds")
# 
# #new phyloseq object:
# taxa.rare <- data.frame(ps.rare.trim@tax_table)
# seqtab.rare <- data.frame(ps.rare.trim@otu_table)
# taxa.rare$sqs <- row.names(taxa.rare) 
# taxa.rare$sqs == colnames(seqtab.rare)
# row.names(taxa.rare) <- taxa.rare$V8
# colnames(seqtab.rare) <- taxa.rare$V8
# row.names(taxa.rare) == colnames(seqtab.rare)
# taxa.rare <- as.matrix(taxa.rare)
# 
# ps.rare.tree <- phyloseq(otu_table(seqtab.rare, taxa_are_rows = FALSE),
#                          sample_data(samdf),
#                          tax_table(taxa.rare),
#                          phy_tree(fitGTR.trim$tree))
# 
# pd.div <- estimate_pd(ps.rare.tree)
# row.names(df.div) <- df.div$id
# df.div.pd <- merge(df.div,pd.div,by=0)


#ATTEMPTING TO MAKE FAITH'S D using picante
library('picante')
library(microbiomeMarker)

#using phyloseq object:
ps.less.nd.no.astreoides

tree.ps.less <- microbiomeMarker::get_treedata_phyloseq(ps.less.nd.no.astreoides, sep = "|")
#tree.ps.less <- phy_tree(ps.less.nd.no.astreoides)
#no tree object already created within the ps object, need to make a new one

#it did not work bc I found out that the only way to create a tree really is through the phang.align and then the NJ step
#which also don't work for other people
#and the only solution I can find online is "MORE RAM" which is stupid

### Faith's D Plots {.tabset}

#subset by species for new dataframe with PD added
df.div.pd.ss <- subset(df.div.pd,host_species=="siderea")
df.div.pd.sr <- subset(df.div.pd,host_species=="radians")
df.div.pd.pp <- subset(df.div.pd,host_species=="porites")

##species comparison faith's d
gg.spp.pd <- ggplot(df.div.pd,aes(x=host_species,y=PD,fill=host_species))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size=22)+
  geom_jitter(alpha=0.5)+
  scale_fill_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))+
  ylab("Phylogenetic Diversity")+
  xlab("Coral Species")+
  theme(legend.position="right",axis.text.x = element_blank(),axis.ticks.x=element_blank())
#facet_wrap(site_zone~.,scales = "free_x")
gg.spp.pd

#siderastrea siderea faith D plot
gg.site.faith.ss <- ggplot(df.div.pd.ss,aes(x=site_zone,y=PD,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 25)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun.y = max, vjust = -0.5,size=8)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual("Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"), labels=c('Santa Martha Bay', 'Santa Martha Reef', 'Spaanse Water Bay','Spannse Water Reef'))+
  ylab("Phylogenetic Diversity")+
  xlab("Site")+
  ylim(5, 55)+
  ggtitle(substitute(paste(italic(" "))))+
  theme(legend.position="right",legend.key.size = unit(2, 'cm'),axis.text.x = element_blank(),axis.ticks.x = element_blank())
  #theme(legend.position="right",legend.key.size = unit(2, 'cm'),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#facet_wrap(m_y~.,scales = "free_x")
gg.site.faith.ss
ggsave(gg.site.faith.ss,file="16S.ss.faith.pdf",h=8,w=8)

#siderastrea radians faithness plot
gg.site.faith.sr <- ggplot(df.div.pd.sr,aes(x=site_zone,y=PD,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  #stat_summary(geom = 'text', label = c("a","a"), fun.y = max, size=6)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual("Site", values = c("coral3","lightsalmon1"),labels=c('Santa Martha Bay','Spaanse Water Bay'))+
  ylab("Phylogenetic Diversity")+
  xlab("Site")+
  ylim(5, 55)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="right", legend.key.size = unit(2, 'cm'), axis.text.x = element_blank(), axis.ticks.x=element_blank())
#facet_wrap(m_y~.,scales = "free_x")
gg.site.faith.sr
#ggsave(gg.site.faith.sr,file="16S.sr.faith.pdf",h=8,w=6)

#porites porites faithness plot
gg.site.faith.pp <- ggplot(df.div.pd.pp,aes(x=site_zone,y=PD,fill=site_zone))+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(alpha=0.5,width = .2)+
  scale_fill_manual("Site", values = c("coral3", "paleturquoise"),labels=c('Santa Martha Bay', 'Santa Martha Reef'))+
  ylab("Phylogenetic Diversity")+
  xlab("Site")+
  ylim(5, 55)+
  ggtitle(substitute(paste(italic(" "))))+
  theme_classic(base_size = 25)+
  theme(legend.position="right", legend.key.size = unit(2, 'cm'), axis.text.x = element_blank(), axis.ticks.x=element_blank())
#facet_wrap(m_y~.,scales = "free_x")
gg.site.faith.pp
#ggsave(gg.site.faith.pp,file="16S.pp.faith.pdf",h=8,w=6)


#MP 25 June 2023
#Separating into groups by samples
#IDs
#samples in bay with majority D vs majority C
#samples in reef with majority D vs majority C

#SMB D
#14, 74, 82, 93, 35, 37, 40, 8, 82, 18, 33, 35, 37, 44, 70, 8, 82, 99
#SMB C
#10, 86, 70, 39
#"10SSSMBMarch2020"|"86SSSMBMarch2020"|"70SSSMBNovember2020"|"39SSSMBNovember2021"
#70 flip flops in between majorities!! probs D bc of copy number thing

#SMR D
#3, 62, 48, 15, 16, 19
#"3SSSMRMarch2020"|"62SSSMRMarch2020"|"48SSSMRNovember2020"|"15SSSMRNovember2021"|"16SSSMRNovember2021"|"19SSSMRNovember2020"
#SMR C
#61, 63, 64, 65, 42, 43, 63, 64, 72, 17, 18, 42, 63, 64, 72

#SWB D
#71, 78, 83, 88, 95, 56, 71, 77, 83, 95, XX, 59, 71, 77, 83, 94, 95
#SWB C
#77 march 2020 (which also flips between times and is probs D)
#"77SSSWBMarch2020"

#SWR D
#1, 49, 96nov2020,85 (96 march 2020 is on edge but slightly more C??? very even tbh so probs D bc CN)
#"1SSDBMarch2020"|"49SSDBNovember2020"|"96SSDBNovember2020"|"85SSDBNovember2021"
#SWR C
#13, 2, 85, 96, 97, 10, 41, 85, 97, 38, 69, 96 97

#MAJORITY C
#"10SSSMBMarch2020"|"86SSSMBMarch2020"|"70SSSMBNovember2020"|"39SSSMBNovember2021"
#"77SSSWBMarch2020"

#MAJORITY D
#"3SSSMRMarch2020"|"62SSSMRMarch2020"|"48SSSMRNovember2020"|"15SSSMRNovember2021"|"16SSSMRNovember2021"|"19SSSMRNovember2020"
##"1SSDBMarch2020"|"49SSDBNovember2020"|"96SSDBNovember2020"|"85SSDBNovember2021"

#SWB 1 C sample, rest D
df.div.ss.swb.cd <- df.div.ss.SWB %>%
  mutate(major_clade = case_when(
    full_sample_id == "77SSSWBMarch2020" ~ "Majority C",
    full_sample_id != "77SSSWBMarch2020" ~ "Majority D"))

#SWR 4 D samples, rest C
df.div.ss.swr.cd <- df.div.ss.SWR %>%
  mutate(major_clade = case_when(
    full_sample_id == "1SSDBMarch2020"|full_sample_id =="49SSDBNovember2020"|full_sample_id =="96SSDBNovember2020"|full_sample_id =="85SSDBNovember2021" ~ "Majority D",
    full_sample_id != "1SSDBMarch2020"|full_sample_id != "49SSDBNovember2020"|full_sample_id != "96SSDBNovember2020"|full_sample_id != "85SSDBNovember2021" ~ "Majority C"))

#SMB 4 C samples, rest D
df.div.ss.smb.cd <- df.div.ss.SMB %>%
  mutate(major_clade = case_when(
    full_sample_id == "100SSSMBMarch2020"|full_sample_id =="86SSSMBMarch2020"|full_sample_id =="70SSSMBNovember2020"|full_sample_id =="39SSSMBNovember2021" ~ "Majority C",
    full_sample_id != "100SSSMBMarch2020"|full_sample_id !="86SSSMBMarch2020"|full_sample_id !="70SSSMBNovember2020"|full_sample_id !="39SSSMBNovember2021" ~ "Majority D"))
#10 was labeled as 100 for some reason? check labeling, i think it should be 100

#SMR - 3 D samples, rest C
df.div.ss.smr.cd <- df.div.ss.SMR %>%
  mutate(major_clade = case_when(
    full_sample_id == "3SSSMRMarch2020"|full_sample_id =="62SSSMRMarch2020"|full_sample_id =="15SSSMRNovember2021" ~ "Majority D",
    full_sample_id != "3SSSMRMarch2020"|full_sample_id !="62SSSMRMarch2020"|full_sample_id !="15SSSMRNovember2021"  ~ "Majority C"))
#sample IDs that are not in the dataset (bc trimmed earlier, remember weird bummer time of sequencing low coverage #tragic #disaster):
#"48SSSMRNovember2020"|"16SSSMRNovember2021"|"19SSSMRNovember2020"

#TOTAL 5 Bay C samples, 7 Reef D samples
#now for some little analyses of alpha diversity (still need to do faith's D in general)

#put the dataframes all back together:
df.div.ss.cd <- rbind(df.div.ss.smb.cd,df.div.ss.smr.cd,df.div.ss.swb.cd,df.div.ss.swr.cd)

#ss.c <- subset(df.div.ss.cd, major_clade=="Majority C") #34 samples
#ss.d <- subset(df.div.ss.cd, major_clade=="Majority D") #40 samples

#Shannon
shapiro.test(df.div.ss.cd$Shannon) #ok
leveneTest(Shannon~major_clade,data=df.div.ss.cd) #ok
sha.div.ss.cd <- aov(Shannon~major_clade,data=df.div.ss.cd)
summary(sha.div.ss.cd) #ns 
TukeyHSD(sha.div.ss.cd)

gg.cd.sha.ss <- ggplot(df.div.ss.cd,aes(x=major_clade,y=Shannon,fill=major_clade))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 22)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun = max, vjust = -0.5, size=6)+
  geom_jitter(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  ylab("Shannon index")+
  xlab("Symbiodiniaceae Genus")+
  #facet_wrap(m_y~.,scales = "free_x")
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg.cd.sha.ss
ggsave(gg.cd.sha.ss,file="cd.ss.sha.pdf",h=8,w=8)

#Simpsons
shapiro.test(df.div.ss.cd$InvSimpson) #sig - bad!
leveneTest(InvSimpson~major_clade,data=df.div.ss.cd) #ok
simp.div.ss.cd <- aov(InvSimpson~major_clade,data=df.div.ss.cd)
summary(simp.div.ss.cd) #ns 
TukeyHSD(simp.div.ss.cd)

gg.cd.simp.ss <- ggplot(df.div.ss.cd,aes(x=major_clade,y=InvSimpson,fill=major_clade))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 22)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun = max, vjust = -0.5, size=6)+
  geom_jitter(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  ylab("Simpson index")+
  xlab("Symbiodiniaceae Genus")+
  #facet_wrap(m_y~.,scales = "free_x")
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg.cd.simp.ss
ggsave(gg.cd.simp.ss,file="cd.ss.simp.pdf",h=8,w=8)

#Evenness
shapiro.test(df.div.ss.cd$even) #sig - bad!
leveneTest(even~major_clade,data=df.div.ss.cd) #ok
even.div.ss.cd <- aov(even~major_clade,data=df.div.ss.cd)
summary(even.div.ss.cd) #ns 
TukeyHSD(even.div.ss.cd)

gg.cd.even.ss <- ggplot(df.div.ss.cd,aes(x=major_clade,y=even,fill=major_clade))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 22)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun = max, vjust = -0.5, size=6)+
  geom_jitter(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  ylab("Evenness")+
  xlab("Symbiodiniaceae Genus")+
  #facet_wrap(m_y~.,scales = "free_x")
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg.cd.even.ss
ggsave(gg.cd.even.ss,file="cd.ss.even.pdf",h=8,w=8)

#Richness
shapiro.test(df.div.ss.cd$Observed) #sig - bad!
leveneTest(Observed~major_clade,data=df.div.ss.cd) #ok
Observed.div.ss.cd <- aov(Observed~major_clade,data=df.div.ss.cd)
summary(Observed.div.ss.cd) #ns 
TukeyHSD(Observed.div.ss.cd)

gg.cd.Observed.ss <- ggplot(df.div.ss.cd,aes(x=major_clade,y=Observed,fill=major_clade))+
  geom_boxplot(outlier.shape=NA,size=1)+
  theme_classic(base_size = 22)+
  #stat_summary(geom = 'text', label = c("a","a","a","b"), fun = max, vjust = -0.5, size=6)+
  geom_jitter(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  ylab("Richness")+
  xlab("Symbiodiniaceae Genus")+
  #facet_wrap(m_y~.,scales = "free_x")
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg.cd.Observed.ss
ggsave(gg.cd.Observed.ss,file="cd.ss.Observed.pdf",h=8,w=8)

###ALL STATS WOOOOOOOOH###

#coral species stats
#shannon
shapiro.test(df.div$Shannon) #ok
leveneTest(Shannon~host_species,data=df.div) #ok
sha.div.spp <- aov(Shannon~host_species,data=df.div)
summary(sha.div.spp) 
#no significant difference
#              Df Sum Sq Mean Sq F value Pr(>F)
#host_species   2    0.5  0.2518   0.187   0.83
#Residuals    120  161.8  1.3486   

#simpson
shapiro.test(df.div$InvSimpson) #very sig
shapiro.test(log(df.div$InvSimpson)) #not sig! use log
leveneTest(log(InvSimpson)~host_species,data=df.div) #ok
simp.div.spp <- aov(log(InvSimpson)~host_species,data=df.div)
summary(simp.div.spp) 
#no significant difference
#              Df Sum Sq Mean Sq F value Pr(>F)
#host_species   2   1.42   0.711   0.657   0.52
#Residuals    120 129.82   1.082   

#richness
shapiro.test(df.div$Observed) #very sig
shapiro.test(log(df.div$Observed)) #still sig, use non-para
leveneTest(log(Observed)~host_species,data=df.div) #ok
kruskal_test(Observed~host_species,data=df.div)
#NS
#.y.          n statistic    df     p method        
#  1 Observed   123    0.0316     2 0.984 Kruskal-Wallis

#dunn_test(df.div.sr, even~site_zone, p.adjust.method = "bonferroni")
#rich.div.spp <- aov(log(Observed)~host_species,data=df.div)
#summary(rich.div.spp) 

#evenness
shapiro.test(df.div$even) #sig
shapiro.test(log(df.div$even)) #very sig - use non para
leveneTest(log(even)~host_species,data=df.div) #ok
kruskal_test(even~host_species,data=df.div)
#ns
#.y.       n statistic    df     p method
#1 even    123      1.65     2 0.437 Kruskal-Wallis
#simp.div.spp <- aov(log(even)~host_species,data=df.div)
#summary(simp.div.spp) 

#faith's d
shapiro.test(df.div.pd$PD) #sig
shapiro.test(log(df.div.pd$PD)) #very sig - use non para
leveneTest(log(PD)~host_species,data=df.div.pd) #ok
#pd.div.spp <- aov(log(InvSimpson)~host_species,data=df.div)
#summary(pd.div.spp) 
kruskal_test(PD~host_species,data=df.div.pd)
#ns
#.y.       n statistic    df     p method 
#1 PD      122      3.70     2 0.157 Kruskal-Wallis

##Shannon STATS!!!###
#siderastrea siderea - shannon stats
shapiro.test(df.div.ss$Shannon) #ok
leveneTest(Shannon~reef_bay/area*m_y,data=df.div.ss) #ok
sha.div.ss <- aov(Shannon~reef_bay/area+m_y,data=df.div.ss)
summary(sha.div.ss) 
#                   Df Sum Sq Mean Sq F value   Pr(>F)    
#reef_bay           1   3.51   3.510   3.059 0.085243 .  
#m_y                2  16.52   8.262   7.200 0.001542 ** 
#reef_bay:area      2  20.79  10.397   9.061 0.000353 ***
#reef_bay:m_y       2   1.70   0.850   0.741 0.480969    
#reef_bay:area:m_y  4   0.78   0.195   0.170 0.953036    
#Residuals         62  71.14   1.147 
TukeyHSD(sha.div.ss)
#ps.less:

#zone/area
#reef:SW-reef:SM -1.5142160 -2.4568843 -0.5715476 0.0004301
#reef:SW-bay:SW  -1.2541304 -2.2105611 -0.2976998 0.0052703
#reef:SW-bay:SM  -1.1432212 -2.0515991 -0.2348432 0.0079663

#timepoint:
#November_2020-March_2020    -1.1571524 -1.9163591 -0.39794574 0.0014984
#November_2021-March_2020    -0.8107399 -1.5262489 -0.09523089 0.0226564

##siderastrea radians - Shannon STATS!!!###
shapiro.test(df.div.sr$Shannon) #ok
leveneTest(Shannon~reef_bay/area*m_y,data=df.div.sr) #ok
sha.div.sr <- aov(Shannon~site_zone*m_y,data=df.div.sr)
summary(sha.div.sr)
#nothing significant

#porites porites - shannon stats
shapiro.test(df.div.pp$Shannon) #ok
leveneTest(Shannon~site_zone*m_y,data=df.div.pp) #ok
sha.div.pp <- aov(Shannon~site_zone*m_y,data=df.div.pp)
summary(sha.div.pp)
#nothing significant!

#porites astreoides - shannon stats
shapiro.test(df.div.pa$Shannon) #ok
leveneTest(Shannon~site_zone,data=df.div.pa) #ok
sha.div.pa <- aov(Shannon~site_zone,data=df.div.pa)
summary(sha.div.pa) 
TukeyHSD(sha.div.pa)
#SWR vs SWB sig * (but  only 5 data points)

#INVERSE SIMPSON STATS
#siderastrea siderea - Inverse Simpson STATS!!!###
shapiro.test(df.div.ss$InvSimpson) #sig! try log
shapiro.test(log(df.div.ss$InvSimpson)) #ns yay
leveneTest(log(InvSimpson)~(reef_bay/area)*m_y,data=df.div.ss) #ok!
simp.div.ss <- aov(log(InvSimpson)~reef_bay/area+m_y,data=df.div.ss)
summary(simp.div.ss) 
#                  Df Sum Sq Mean Sq F value Pr(>F)   
#reef_bay           1   1.77   1.765   1.811 0.1833   
#m_y                2   7.60   3.799   3.897 0.0254 * 
#reef_bay:area      2  11.04   5.522   5.665 0.0055 **
#reef_bay:m_y       2   2.64   1.321   1.355 0.2654   
#reef_bay:area:m_y  4   1.16   0.289   0.297 0.8790   
#Residuals         62  60.43   0.975  
TukeyHSD(simp.div.ss)
#site significance
#reef:SW-reef:SM -1.0895853 -1.9583929 -0.22077759 0.0082422
#reef:SW-bay:SW  -0.9587991 -1.8402908 -0.07730742 0.0278477

#timepoint significance
#reef:November_2020-bay:March_2020     -1.299746968 -2.5116079 -0.08788608 0.0285908

#interactive significance #not interesting too bad!
#reef:SW:November_2020-reef:SM:March_2020    -2.18636300 -4.1219099 -0.25081614 0.0142950

#siderastrea radians - Inverse Simpson STATS!!!###
shapiro.test(df.div.sr$InvSimpson) #sig! try log
shapiro.test(log(df.div.sr$InvSimpson)) #ns yay
leveneTest(log(InvSimpson)~site_zone*m_y,data=df.div.sr) #ok!
simp.div.sr <- aov(log(InvSimpson)~site_zone*m_y,data=df.div.sr)
summary(simp.div.sr) 
#ns

#porites porites - Inverse simpson STATS
#only keeps 3 samples for SM_reef for final
#df.div.pp$lane <- c("1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2")
shapiro.test(df.div.pp$InvSimpson) #sig! try log
shapiro.test(log(df.div.pp$InvSimpson)) #ns yay
leveneTest(log(InvSimpson)~site_zone*m_y,data=df.div.pp) #ok!
simp.div.pp <- aov(log(InvSimpson)~site_zone*m_y,data=df.div.pp)
summary(simp.div.pp) 
#ns

#porites astreoides - Inverse Simpson stats
shapiro.test(df.div.pa$InvSimpson) #sig! try log
shapiro.test(log(df.div.pa$InvSimpson)) #ns yay
leveneTest(log(InvSimpson)~site_zone*m_y,data=df.div.pa) #ok!
simp.div.pa <- aov(log(InvSimpson)~site_zone,data=df.div.pa)
summary(simp.div.pa)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#site_zone    1  6.927   6.927   24.17 0.00795 **
#Residuals    4  1.146   0.287  
#only two sites so this is the p-value for comparison too

#RICHNESS STATS
#siderastrea siderea - richness STATS!!!###
shapiro.test(df.div.ss$Observed) #sig! try log
shapiro.test(log(df.div.ss$Observed)) #SIG - use non para
leveneTest(log(Observed)~reef_bay/area*m_y,data=df.div.ss) #ok!
rich.div.ss <- aov(log(Observed)~reef_bay/area*m_y,data=df.div.ss)
summary(rich.div.ss) 
#                  Df Sum Sq Mean Sq F value   Pr(>F)    
#reef_bay           1   1.80   1.802   2.007  0.16161    
#m_y                2  23.23  11.614  12.931 2.02e-05 ***
#reef_bay:area      2  10.90   5.449   6.067  0.00392 ** 
#reef_bay:m_y       2   0.30   0.149   0.166  0.84732    
#reef_bay:area:m_y  4   0.07   0.018   0.020  0.99920    
#Residuals         62  55.68   0.898                      
TukeyHSD(rich.div.ss)

#site significant difs
#reef:SW-bay:SM  -0.88380842 -1.6874733 -0.08014353 0.0256191
#reef:SW-reef:SM -1.09789353 -1.9318960 -0.26389106 0.0050586

#timepoint sig difs
#November_2020-March_2020    -1.4127466 -2.08443594 -0.7410573 0.0000123
#November_2021-March_2020    -0.8109043 -1.44393317 -0.1778754 0.0086433

#S siderea actual stats for paper:
shapiro.test(df.div.ss$Observed) #very sig
shapiro.test(log(df.div.ss$Observed)) #still sig, use non-para
leveneTest(log(Observed)~reef_bay/area*m_y,data=df.div.ss) #ok!
kruskal_test(Observed~r_b_m_y,data=df.div.ss)
#reef_bay
#.y.          n statistic    df     p method  
#1 Observed    74      2.04     1 0.153 Kruskal-Wallis
#m_y
#Observed    74      19.6     2 0.0000559 Kruskal-Wallis
#site_zone (aka reef_bay:area)
#Observed    74      11.9     3 0.00786 Kruskal-Wallis
#r_b_m_y
#Observed    74      22.1     5 0.000507 Kruskal-Wallis
#m_y_s_z (aka reef_bay:area:m_y)
#Observed    74      31.6    11 0.000899 Kruskal-Wallis

dunn_test(df.div.ss, Observed~m_y_s_z, p.adjust.method = "bonferroni")
#m_y
#  .y.      group1        group2           n1    n2 statistic         p     p.adj p.adj.signif
#  1 Observed March_2020    November_2020    24    22     -4.39 0.0000115 0.0000345 ****        
#  2 Observed March_2020    November_2021    24    28     -2.73 0.00633   0.0190    *           
#  3 Observed November_2020 November_2021    22    28      1.88 0.0603    0.181     ns 

#site_zone (aka reef_bay:area)
#  .y.      group1  group2     n1    n2 statistic       p  p.adj p.adj.signif
#1 Observed SM_bay  SM_reef    21    18     0.663 0.507   1      ns          
#2 Observed SM_bay  SW_bay     21    17     0.158 0.874   1      ns          
#3 Observed SM_bay  SW_reef    21    18    -2.59  0.00970 0.0582 ns          
#4 Observed SM_reef SW_bay     18    17    -0.477 0.634   1      ns          
#5 Observed SM_reef SW_reef    18    18    -3.13  0.00174 0.0105 *           
#6 Observed SW_bay  SW_reef    17    18    -2.61  0.00907 0.0544 ns  

#r_b_m_y
#   .y.      group1             group2           n1    n2 statistic       p   p.adj p.adj.signif
# * <chr>    <chr>              <chr>         <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
#  1 Observed bay March_2020     bay November…    12    11    -3.46  5.50e-4 8.25e-3 **          
#  2 Observed bay March_2020     bay November…    12    15    -2.40  1.65e-2 2.48e-1 ns          
#  3 Observed bay March_2020     reef March_2…    12    12    -1.34  1.81e-1 1   e+0 ns          
#  4 Observed bay March_2020     reef Novembe…    12    11    -4.06  4.96e-5 7.44e-4 ***         
#  5 Observed bay March_2020     reef Novembe…    12    13    -2.88  3.98e-3 5.96e-2 ns          
#  6 Observed bay November_2020  bay November…    11    15     1.29  1.96e-1 1   e+0 ns          
#  7 Observed bay November_2020  reef March_2…    11    12     2.15  3.19e-2 4.78e-1 ns          
#  8 Observed bay November_2020  reef Novembe…    11    11    -0.590 5.55e-1 1   e+0 ns          
#  9 Observed bay November_2020  reef Novembe…    11    13     0.706 4.80e-1 1   e+0 ns          
# 10 Observed bay November_2021  reef March_2…    15    12     0.987 3.24e-1 1   e+0 ns          
# 11 Observed bay November_2021  reef Novembe…    15    11    -1.93  5.39e-2 8.08e-1 ns          
# 12 Observed bay November_2021  reef Novembe…    15    13    -0.592 5.54e-1 1   e+0 ns          
# 13 Observed reef March_2020    reef Novembe…    12    11    -2.75  5.98e-3 8.98e-2 ns          
# 14 Observed reef March_2020    reef Novembe…    12    13    -1.52  1.30e-1 1   e+0 ns          
# 15 Observed reef November_2020 reef Novembe…    11    13     1.32  1.87e-1 1   e+0 ns   

#siderastrea radians - richness STATS!!!###
shapiro.test(df.div.sr$Observed) #sig! try log
shapiro.test(log(df.div.sr$Observed)) #ns yay
leveneTest(log(Observed)~site_zone*m_y,data=df.div.sr) #ok!
rich.div.sr <- aov(log(Observed)~site_zone*m_y,data=df.div.sr)
summary(rich.div.sr) 
#ns

#porites porites - richness stats
shapiro.test(df.div.pp$Observed) #sig! try log
shapiro.test(log(df.div.pp$Observed)) #ns yay
leveneTest(log(Observed)~site_zone*m_y,data=df.div.pp) #ok!
rich.div.pp <- aov(log(Observed)~site_zone*m_y,data=df.div.pp)
summary(rich.div.pp) 
#               Df Sum Sq Mean Sq F value Pr(>F)  
#site_zone      1  0.007   0.007   0.010 0.9228  
#m_y            1  5.927   5.927   8.073 0.0149 *
#site_zone:m_y  1  0.002   0.002   0.003 0.9561  
#Residuals     12  8.811   0.734 
TukeyHSD(rich.div.pp)
#site - all ns
#timepoint
#November_2021-November_2020 -1.520532 -2.716383 -0.3246813 0.0169501
#no interactive significance

#porites astreoides - richness stats
shapiro.test(df.div.pa$Observed) #ok
leveneTest(Observed~site_zone,data=df.div.pa) #ok
aov.pa.rich <- aov(Observed~site_zone,data=df.div.pa)
summary(aov.pa.rich) #sig
TukeyHSD(aov.pa.rich)

#EVENNESS STATS:
#siderastrea siderea evenness STATS!!!

shapiro.test(df.div.ss$even) #sig! try log
shapiro.test(log(df.div.ss$even)) #VERY SIG - do non para
leveneTest(even~site_zone*m_y,data=df.div.ss) #ok!
even.div.ss <- aov(even~reef_bay/area*m_y,data=df.div.ss)
summary(even.div.ss) 
#normal aov finds just site-zone sig dif
#                  Df Sum Sq Mean Sq F value Pr(>F)  
#reef_bay           1 0.0445 0.04449   2.945 0.0911 .
#m_y                2 0.0166 0.00831   0.550 0.5797  
#reef_bay:area      2 0.1416 0.07081   4.687 0.0127 *
#reef_bay:m_y       2 0.0773 0.03866   2.559 0.0855 .
#reef_bay:area:m_y  4 0.0317 0.00793   0.525 0.7178  
#Residuals         62 0.9367 0.01511    
TukeyHSD(even.div.ss)
#site
#reef:SW-reef:SM -0.1230077387 -0.23117882 -0.014836653 0.0196659
#reef:SW-bay:SW  -0.1235623574 -0.23331267 -0.013812049 0.0213173

#no timepoint or interactive are significant

#siderastrea radians evenness stats
shapiro.test(df.div.sr$even) #ns
leveneTest(even~site_zone,data=df.div.sr) #ok
#kruskal_test(even~site_zone,data=df.div.sr)
#dunn_test(df.div.sr, even~site_zone, p.adjust.method = "bonferroni")
even.div.sr <- aov(even~site_zone*m_y,data=df.div.sr)
summary(even.div.sr) #ns

#porites porites evenness stats
shapiro.test(df.div.pp$even) #ok
leveneTest(even~site_zone,data=df.div.pp) #ok
even.div.pp <- aov(even~site_zone*m_y,data=df.div.pp)
summary(even.div.pp)#ns

#porites astreoides evenness stats
shapiro.test(df.div.pa$even) #ok
leveneTest(even~site_zone,data=df.div.pa) #ok
aov.pa.even <- aov(even~site_zone,data=df.div.pa)
summary(aov.pa.even) #sig
TukeyHSD(aov.pa.even)

###Faith's D Phylogenetic Diversity STATS:
#siderastrea siderea PD STATS!!!

shapiro.test(df.div.pd.ss$PD) #almost not sig, just use
shapiro.test(log(df.div.pd.ss$PD)) #VERY SIG - do non para or just use reg
leveneTest(PD~site_zone*m_y,data=df.div.pd.ss) #ok!
PD.div.ss <- aov(PD~reef_bay/area*m_y,data=df.div.pd.ss)
summary(PD.div.ss) 

#normal aov finds site-zone and time sig dif
#Df Sum Sq Mean Sq F value  Pr(>F)   
#reef_bay           1    292   292.0   2.631 0.10986   
#m_y                2   1344   672.0   6.055 0.00396 **
#reef_bay:area      2   1090   545.0   4.911 0.01047 * 
#reef_bay:m_y       2    102    51.2   0.462 0.63246   
#reef_bay:area:m_y  4     34     8.4   0.076 0.98936   
#Residuals         62   6880   111.0    
TukeyHSD(PD.div.ss)

#site
#reef:SW-reef:SM -0.1230077387 -0.23117882 -0.014836653 0.0196659
#reef:SW-bay:SW  -0.1235623574 -0.23331267 -0.013812049 0.0213173

#non-para - use these for final paper!
kruskal_test(PD~site_zone,data=df.div.pd.ss)
#  .y.       n statistic    df      p method  
#1 PD       74      11.0     3 0.0117 Kruskal-Wallis
dunn_test(df.div.pd.ss, PD~m_y, p.adjust.method = "bonferroni")
#> dunn_test(df.div.pd.ss, PD~site_zone, p.adjust.method = "bonferroni")
# A tibble: 6 × 9
#.y.   group1  group2         n1    n2 statistic   p  p.adj p.adj.signif
#* <chr> <chr>   <chr>   <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
#  1 PD    SM_bay  SM_reef    21    18     0.223 0.824   1      ns          
#  2 PD    SM_bay  SW_bay     21    17    -0.450 0.652   1      ns          
#  3 PD    SM_bay  SW_reef    21    18    -2.81  0.00497 0.0298 *           
#  4 PD    SM_reef SW_bay     18    17    -0.646 0.518   1      ns          
#  5 PD    SM_reef SW_reef    18    18    -2.92  0.00348 0.0209 *           
#  6 PD    SW_bay  SW_reef    17    18    -2.23  0.0255  0.153  ns 
dunn_test(df.div.pd.ss, PD~m_y, p.adjust.method = "bonferroni")
#1 PD    March_2020    November_2020    24    22     -3.00 0.00267 0.00801 **          
#2 PD    March_2020    November_2021    24    28     -1.16 0.246   0.739   ns          
#3 PD    November_2020 November_2021    22    28      1.98 0.0477  0.143   ns   
#march 2020 significant from nov 2020

#no timepoint or interactive are significant

#siderastrea radians PDness stats
shapiro.test(df.div.pd.sr$PD) #ns
leveneTest(PD~site_zone,data=df.div.pd.sr) #ok
PD.div.sr <- aov(PD~area*m_y,data=df.div.pd.sr)
summary(PD.div.sr) #ns

#porites porites PDness stats
shapiro.test(df.div.pd.pp$PD) #ok
leveneTest(PD~site_zone,data=df.div.pd.pp) #ok
PD.div.pp <- aov(PD~reef_bay*m_y,data=df.div.pd.pp)
summary(PD.div.pp) 
TukeyHSD(PD.div.pp)

#teehee trying glms bc I think all of my stats are wrong xoxo gopissgirl

#updates 7/15/24 with help from Rachael on how to actually do the stats!!
#here's the plan:
#1. Run glmms using glmmTMB - use all interactive comparisons to check distributions
#2. Check residuals using simulateResiduals
#3. For distributions - use gaussian, and if less than 4 outliers and meets assumptions - keep
#otherwise then try log of response, gamma distribution, or try t_family
#note that t_family does great w/outliers but has a really hard time with degrees of freedom - ask Rachael
#4. For parameter selection - no matter what, include genotype as random effect
#then select using AICc
#useful help page: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html!!
#example: 
type3 <- list(m_y = contr.sum, reef_bay = contr.sum, area = contr.sum)
glm.3e.1 <- glmmTMB(Shannon ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian(), contrasts = )
#CHECKING ASSUMPTIONS AGAIN 7/15/24 with Rachael
#using DHARMa package
library("DHARMa") #0.4.6
sim.resids <- simulateResiduals(glm.3e.1, plot = TRUE)
outliers(sim.resids) #only 2 - check which genotypes these are and compare with alpha div data
Anova(glm.3e.1, type = "III")
#N = 74
#using type 3 anova becuase I have an unbalanced design and interested in interactions
#check with book tmrw but I think this is right?
#read somewhere that saying type II or III in lm output is garbage bc you didn't actually specify it in the model
marginal = emmeans(glm.3e.1, ~ m_y*reef_bay*area)
pairs(marginal)
emmeans(glm.3e.1)
#likelihood ratio test - do comparison of the less interactive models

library(stats)
#install.packages("MuMIn")
library(MuMIn)
#library(lmer)

#normality of all data for species comparisons
hist(df.div$Shannon) #normal
shapiro.test(df.div$Shannon) #normal - use gaussian dist
hist(df.div$InvSimpson) #right skewed
shapiro.test(df.div$InvSimpson) #sig/not normal
hist(log(df.div$InvSimpson)) #normal
shapiro.test(log(df.div$InvSimpson)) #normal - use gaussian for log(invsimp)
hist(df.div$Observed) #right skewed
shapiro.test(df.div$Observed) #very sig/not normal
hist(log(df.div$Observed)) #kinda normal tbh or maybe bimodal?
shapiro.test(log(df.div$Observed)) #still sig - need different distribution
hist(df.div$even) #left skewed
shapiro.test(df.div$even) #sig/not normal - need different distribution
hist(log(df.div$even)) #even more left skewed
shapiro.test(log(df.div$even)) #even more sig
hist(df.div.pd$PD) #left skewed
shapiro.test(df.div.pd$PD) #so close almost normal 0.02264
hist(log(df.div.pd$PD)) #even more left skewed
shapiro.test(log(df.div.pd$PD)) #more sig!

#choosing distributions for all data species comparisons:
#gaussian for Shannon
#gaussian for log of Inverse Simpsons
#something new for Observed - and probably log of observed
#something new for evenness
#something new for PD - almost normal, maybe just dif link?

##GLMs for all species shannon

glm.0 <- glmmTMB(Shannon ~ 1, data=df.div, family=gaussian())
glm.1a <- glmmTMB(Shannon ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian())

sim.resids <- simulateResiduals(glm.1a, plot = TRUE) #looks great

AICc(glm.0,glm.1a)
#best = glm.0, predicted by null model
#Anova(glm.1a, type = "III")

#GLMs of log(InvSimpson) for all species
glm.0 <- glmmTMB(log(InvSimpson) ~ 1, data=df.div, family=gaussian())
glm.1a <- glmmTMB(log(InvSimpson) ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian())

sim.resids <- simulateResiduals(glm.1a, plot = TRUE) #looks great
AICc(glm.0,glm.1a) #null model best

#for richness
glm.0 <- glmmTMB(log(Observed) ~ 1, data=df.div, family=gaussian())
glm.1a <- glmmTMB(log(Observed) ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian())

sim.resids <- simulateResiduals(glm.1a, plot = TRUE) #looks great
AICc(glm.0,glm.1a)
#best = glm.0, predicted by null model

#something new for evenness
glm.a <- glmmTMB(even ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian())
glm.b <- glmmTMB(even ~ host_species+(1|site:no_sp_lo), data=df.div, family=Gamma())
glm.d <- glmmTMB(even ~ host_species+(1|site:no_sp_lo), data=df.div)
glm.f <- glmmTMB(even ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian(link = "inverse"))

AICc(glm.a, glm.b, glm.d, glm.e, glm.f)
#glm.a is best
glm.0 <- glmmTMB(even ~ 1, data=df.div, family=gaussian())
glm.1a <- glmmTMB(even ~ host_species+(1|site:no_sp_lo), data=df.div, family=gaussian())

sim.resids <- simulateResiduals(glm.1a, plot = TRUE) #looks great
AICc(glm.0,glm.1a)
#best = glm.0, predicted by null model

#something new for PD - almost normal, maybe just dif link?
glm.a <- glmmTMB(PD ~ host_species+(1|site:no_sp_lo), data=df.div.pd, family=gaussian())
glm.b <- glmmTMB(PD ~ host_species+(1|site:no_sp_lo), data=df.div.pd, family=Gamma())
glm.d <- glmmTMB(PD ~ host_species+(1|site:no_sp_lo), data=df.div.pd)
glm.f <- glmmTMB(PD ~ host_species+(1|site:no_sp_lo), data=df.div.pd, family=gaussian(link = "inverse"))

AICc(glm.a, glm.b, glm.d, glm.f)
#glm.a is best

glm.0 <- glmmTMB(PD ~ 1, data=df.div.pd, family=gaussian())
glm.1a <- glmmTMB(PD ~ host_species+(1|site:no_sp_lo), data=df.div.pd, family=gaussian())

sim.resids <- simulateResiduals(glm.1a, plot = TRUE) #looks great
AICc(glm.0,glm.1a)
#best = glm.0, predicted by null model

#qqPlot(df.div$Observed)
#ggdensity(df.div$Observed) #right skewed

##normality of all data for ssid
hist(df.div.ss$Shannon) #normal enough
shapiro.test(df.div.ss$Shannon) #normal - use gaussian
hist(df.div.ss$InvSimpson) #right skewed
shapiro.test(df.div.ss$InvSimpson) #sig
hist(log(df.div.ss$InvSimpson)) #normalish
shapiro.test(log(df.div.ss$InvSimpson)) #normal - use gaussian
hist(df.div.ss$Observed) #right skewed
shapiro.test(df.div.ss$Observed) #sig
hist(log(df.div.ss$Observed)) #bimodal?
shapiro.test(log(df.div.ss$Observed)) #sig - use different distribution
hist(df.div.ss$even) #left skewed
shapiro.test(df.div.ss$even) #sig - use different distribution
hist(df.div.pd.ss$PD) #weird kinda flat
shapiro.test(df.div.pd.ss$PD) #so close almost normal 0.025
hist(log(df.div.pd.ss$PD)) #even more left skewed
shapiro.test(log(df.div.pd.ss$PD)) #more sig!

#choosing distributions for ssid:
#gaussian for Shannon
#gaussian for log of Inverse Simpsons
#something new for Observed - and probably log of observed
#something new for evenness
#something new for PD - almost normal

#GLMs for ssid shannon
glm.0 <- glmmTMB(Shannon ~ 1+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1a <- glmmTMB(Shannon ~ reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1b <- glmmTMB(Shannon ~ m_y+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1c <- glmmTMB(Shannon ~ area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.2a <- glmmTMB(Shannon ~ m_y+reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2b <- glmmTMB(Shannon ~ m_y+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2c <- glmmTMB(Shannon ~ reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.3a <- glmmTMB(Shannon ~ m_y+reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
#edited this model - make sure to check when running through
#ok so during convo with Sarah and Verena they said I should use the model with just reef_bay
#and not include area, but for me, area explains a lot of difference, so I am including it here
#bc it is important otherwise the best model is just m_y only and no affect of any location differences
glm.3b <- glmmTMB(Shannon ~ m_y+reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3c <- glmmTMB(Shannon ~ m_y*reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3d <- glmmTMB(Shannon ~ reef_bay+m_y*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3e <- glmmTMB(Shannon ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.1c, glm.2a,glm.2b,glm.2c, glm.3a, glm.3b, glm.3c,glm.3d,glm.3e)
#best = glm.3b AICc = 225.9361
sim.resids <- simulateResiduals(glm.3b, plot = TRUE) #looks good, only 1 outlier!
#no significant problems detected
outliers(sim.resids)
#30 is the outlier

Anova(glm.3b, type = "III")

marginal = emmeans(glm.3b, ~ reef_bay*area)
pairs(marginal)

#habitat with site nested m_y+area*reef_bay+(1|site:no_sp_lo)
#sarah's model:
#model5.8<-glmer(POR$POR_SHIFT_d13C.he ~habitat*time+(1|site:geno.unique),family= Gamma(link = "log"), data=POR)

#                 Chisq Df Pr(>Chisq)    
#(Intercept)   284.5212  1  < 2.2e-16 ***
#m_y            22.9844  2  1.021e-05 ***
#reef_bay        0.5581  1    0.45505    
#area            0.0044  1    0.94724    
#reef_bay:area   6.1848  1    0.01289 * 

#GLMs of log(InvSimpson) for SSID
glm.0 <- glmmTMB(log(InvSimpson) ~ 1+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1a <- glmmTMB(log(InvSimpson) ~ reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1b <- glmmTMB(log(InvSimpson) ~ m_y+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1c <- glmmTMB(log(InvSimpson) ~ area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.2a <- glmmTMB(log(InvSimpson) ~ m_y+reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2b <- glmmTMB(log(InvSimpson) ~ m_y+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2c <- glmmTMB(log(InvSimpson) ~ site_zone+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.3a <- glmmTMB(log(InvSimpson) ~ m_y+reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3b <- glmmTMB(log(InvSimpson) ~ m_y+reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3c <- glmmTMB(log(InvSimpson) ~ m_y*reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3d <- glmmTMB(log(InvSimpson) ~ reef_bay+m_y*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3e <- glmmTMB(log(InvSimpson) ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.1c, glm.2a,glm.2b,glm.2c, glm.3a, glm.3b, glm.3c,glm.3d,glm.3e)
#glm3b is best again wooh AICc 215.2185
sim.resids <- simulateResiduals(glm.3b, plot = TRUE) #no significant problems - 1 outlier
outliers(sim.resids) #30 same one
Anova(glm.3b, type = "III")
#(Intercept)   179.2614  1    < 2e-16 ***
#m_y             9.3662  2    0.00925 ** 
#reef_bay        0.6770  1    0.41063    
#area            0.2185  1    0.64017    
#reef_bay:area   5.4691  1    0.01936 *  

marginal = emmeans(glm.3b, ~ reef_bay*area)
pairs(marginal)

marginal = emmeans(glm.3b, ~ m_y)
pairs(marginal)

#GLMs for ssid Observed
#expecting this model to be the best so using it to compare distributions/links
glm.a <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.b <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "inverse"))
glm.c <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "identity"))
glm.d <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "log"))
glm.e <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=inverse.gaussian())
glm.f <- glmmTMB(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian(link = "log"))
#glm.g <- glm.nb(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss)
#glm.h <- glm.nb(log(Observed) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss)

AICc(glm.a, glm.b,glm.c,glm.d,glm.f)
#glm.a = best - AICc = 205.1460

#GLMs of log(Observed) for SSID
glm.0 <- glmmTMB(log(Observed) ~ 1+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1a <- glmmTMB(log(Observed) ~ reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1b <- glmmTMB(log(Observed) ~ m_y+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1c <- glmmTMB(log(Observed) ~ area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.2a <- glmmTMB(log(Observed) ~ m_y+reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2b <- glmmTMB(log(Observed) ~ m_y+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2c <- glmmTMB(log(Observed) ~ site_zone+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.3a <- glmmTMB(log(Observed) ~ m_y+reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3b <- glmmTMB(log(Observed) ~ m_y+reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3c <- glmmTMB(log(Observed) ~ m_y*reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3d <- glmmTMB(log(Observed) ~ reef_bay+m_y*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3e <- glmmTMB(log(Observed) ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.1c, glm.2a,glm.2b,glm.2c, glm.3a, glm.3b, glm.3c,glm.3d,glm.3e)
#glm3b is best again wooh AICc 205.1460
sim.resids <- simulateResiduals(glm.3b, plot = TRUE) #no significant issues, one outlier!
outliers(sim.resids) #57

Anova(glm.3b, type = "III")
#(Intercept)   719.8392  1  < 2.2e-16 ***
#m_y            29.8453  2  3.305e-07 ***
#reef_bay        0.5572  1   0.455397    
#area            0.0340  1   0.853686    
#reef_bay:area   6.6554  1   0.009886 ** 

marginal = emmeans(glm.3b, ~ m_y)
pairs(marginal)

marginal = emmeans(glm.3b, ~ reef_bay*area)
pairs(marginal)

#GLMs evenness
#expecting this model to be the best so using it to compare distributions/links
glm.a <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian(link = "identity"))
glm.b <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "inverse"))
glm.c <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "identity"))
glm.d <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=Gamma(link = "log"))
glm.e <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian(link = "inverse"))
glm.f <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian(link = "log"))
glm.g <- glm.nb(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss) #def not nb for this
glm.h <- glmmTMB(even ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.ss, family=inverse.gaussian())

AICc(glm.a, glm.b,glm.c,glm.d,glm.e,glm.f)
#glm.a = best - AICc = -89.34095

#GLMs of even for SSID
glm.0 <- glmmTMB(even ~ 1+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1a <- glmmTMB(even ~ reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1b <- glmmTMB(even ~ m_y+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.1c <- glmmTMB(even ~ area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.2a <- glmmTMB(even ~ m_y+reef_bay+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2b <- glmmTMB(even ~ m_y+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.2c <- glmmTMB(even ~ reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

glm.3a <- glmmTMB(even ~ m_y+reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3b <- glmmTMB(even ~ m_y+reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3c <- glmmTMB(even ~ m_y*reef_bay+area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3d <- glmmTMB(even ~ reef_bay+m_y*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())
glm.3e <- glmmTMB(even ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.ss, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.1c, glm.2a,glm.2b,glm.2c, glm.3a, glm.3b, glm.3c,glm.3d,glm.3e)
#best = glm.0: AICc -92.31296
#null model - no signficant interactions for S. siderea evenness

#Faith's D Phylogenetic diversity for S. siderea
#expecting this model to be the best so using it to compare distributions/links
glm.a <- glmmTMB(PD ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian(link = "identity"))
glm.b <- glmmTMB(log(PD) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=Gamma(link = "inverse"))
glm.c <- glmmTMB(log(PD) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=Gamma(link = "identity"))
glm.d <- glmmTMB(log(PD) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=Gamma(link = "log"))
glm.e <- glmmTMB(log(PD) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian(link = "identity"))
glm.f <- glmmTMB(log(PD) ~ m_y+reef_bay/area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian(link = "log"))

AICc(glm.a, glm.b,glm.c,glm.d,glm.e,glm.f)
#glm.e = best - AICc = 98.8967

#GLMs of PD for SSID
glm.0 <- glmmTMB(log(PD) ~ 1+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.1a <- glmmTMB(log(PD) ~ reef_bay+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.1b <- glmmTMB(log(PD) ~ m_y+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.1c <- glmmTMB(log(PD) ~ area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())

glm.2a <- glmmTMB(log(PD) ~ m_y+reef_bay+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.2b <- glmmTMB(log(PD) ~ m_y+area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.2c <- glmmTMB(log(PD) ~ reef_bay+area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())

glm.3a <- glmmTMB(log(PD) ~ m_y+reef_bay+area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.3b <- glmmTMB(log(PD) ~ m_y+reef_bay*area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.3c <- glmmTMB(log(PD) ~ m_y*reef_bay+area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.3d <- glmmTMB(log(PD) ~ reef_bay+m_y*area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())
glm.3e <- glmmTMB(log(PD) ~ m_y*reef_bay*area+(1|site:no_sp_lo), data=df.div.pd.ss, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.1c, glm.2a,glm.2b,glm.2c, glm.3a, glm.3b, glm.3c,glm.3d,glm.3e)
#best = glm.3b AICc = 98.89670

sim.resids <- simulateResiduals(glm.3b, plot = TRUE) #no significant issues, one outlier!
outliers(sim.resids) #59, 66

Anova(glm.3b, type = "III")
#(Intercept)   914.3698  1  < 2.2e-16 ***
#m_y            16.3963  2  0.0002752 ***
#reef_bay        0.0949  1  0.7580737    
#area            0.5847  1  0.4444819    
#reef_bay:area   3.5478  1  0.0596256 . 

marginal = emmeans(glm.3b, ~ m_y)
pairs(marginal)

marginal = emmeans(glm.3b, ~ reef_bay*area)
pairs(marginal)

###normality of all data for srad
hist(df.div.sr$Shannon) #normal enough
shapiro.test(df.div.sr$Shannon) #0.063 normal - use gaussian
hist(df.div.sr$InvSimpson) #right skewed
shapiro.test(df.div.sr$InvSimpson) #sig
hist(log(df.div.sr$InvSimpson)) #normalish
shapiro.test(log(df.div.sr$InvSimpson)) #normal - use gaussian of log
hist(df.div.sr$Observed) #right skewed
shapiro.test(df.div.sr$Observed) #sig
hist(log(df.div.sr$Observed)) #bimodal?
shapiro.test(log(df.div.sr$Observed)) #normal - use gaussian of log
hist(df.div.sr$even) #left skewed
shapiro.test(df.div.sr$even) #normal - use gaussian
hist(df.div.pd.sr$PD) #weird kinda flat
shapiro.test(df.div.pd.sr$PD) #normal - use gaussian
hist(log(df.div.pd.sr$PD)) #even more left skewed
shapiro.test(log(df.div.pd.sr$PD)) #also normal but just use regular yay

#choosing distributions for srad:
#gaussian for Shannon
#gaussian for log of Inverse Simpsons
#gaussian for log of Observed
#gaussian for even
#gaussian for faith's d

##GLMs for srad all variables
glm.0 <- glmmTMB(log(PD) ~ 1+(1|site:no_sp_lo), data=df.div.pd.sr, family=gaussian())
glm.1a <- glmmTMB(log(PD) ~ area+(1|site:no_sp_lo), data=df.div.pd.sr, family=gaussian())
glm.1b <- glmmTMB(log(PD) ~ m_y+(1|site:no_sp_lo), data=df.div.pd.sr, family=gaussian())

glm.2a <- glmmTMB(log(PD) ~ area+m_y+(1|site:no_sp_lo), data=df.div.pd.sr, family=gaussian())

glm.3a <- glmmTMB(log(PD) ~ area*m_y+(1|site:no_sp_lo), data=df.div.pd.sr, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.2a,glm.3a)
#best glm.0 for Shannon, log(InvSimpson), log(Observed), even
#Faith's D! log - glm.2a is best!!

sim.resids <- simulateResiduals(glm.2a, plot = TRUE) #no significant issues, one outlier!
outliers(sim.resids) #no outliers COME ON!!!!

Anova(glm.2a, type = "III")
# Response: log(PD)
# Chisq Df Pr(>Chisq)    
# (Intercept) 983.1468  1    < 2e-16 ***
#   area          3.8748  1    0.04902 *  
#   m_y           8.5816  2    0.01369 *

marginal = emmeans(glm.2a, ~ m_y)
pairs(marginal)

marginal = emmeans(glm.2a, ~ area)
pairs(marginal)

###normality of all data for porites
hist(df.div.pp$Shannon) #normal enough
shapiro.test(df.div.pp$Shannon) #ormal - use gaussian
hist(df.div.pp$InvSimpson) #right skewed
shapiro.test(df.div.pp$InvSimpson) #sig
hist(log(df.div.pp$InvSimpson)) #normalish
shapiro.test(log(df.div.pp$InvSimpson)) #normal - use gaussian of log
hist(df.div.pp$Observed) #right skewed
shapiro.test(df.div.pp$Observed) #sig
hist(log(df.div.pp$Observed)) #bimodal?
shapiro.test(log(df.div.pp$Observed)) #normal - use gaussian of log
hist(df.div.pp$even) #left skewed
shapiro.test(df.div.pp$even) #normal - use gaussian
hist(df.div.pd.pp$PD) #weird kinda flat
shapiro.test(df.div.pd.pp$PD) #normal - use gaussian
hist(log(df.div.pd.pp$PD)) #even more left skewed
shapiro.test(log(df.div.pd.pp$PD)) #also normal but just use regular yay

##choosing distributions for ppor:
#gaussian for Shannon
#gaussian for log of Inverse Simpsons
#gaussian for log of Observed
#gaussian for even
#gaussian for faith's d

#GLMs for ppor shannon
glm.0 <- glmmTMB(even ~ 1+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())
glm.1a <- glmmTMB(even ~ reef_bay+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())
glm.1b <- glmmTMB(even ~ m_y+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())

glm.1b <- glmmTMB(even ~ m_y+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())

glm.2a <- glmmTMB(even ~ reef_bay+m_y+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())

glm.3a <- glmmTMB(even ~ reef_bay*m_y+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())

AICc(glm.0,glm.1a,glm.1b,glm.2a,glm.3a)
#best glm.0 for Shannon, log(InvSimpson), even #bummer but thats ok
#best glm.2a for log(Observed)
glm.2a <- glmmTMB(log(Observed) ~ reef_bay+m_y+(1|site:no_sp_lo), data=df.div.pp, family=gaussian())
sim.resids <- simulateResiduals(glm.2a, plot = TRUE) #no issues/outliers

Anova(glm.2a, type = "III")
# Response: log(Observed)
# Chisq Df Pr(>Chisq)    
# (Intercept) 205.2935  1  < 2.2e-16 ***
#   reef_bay      0.3712  1  0.5423373    
#   m_y          11.8117  1  0.0005886 ***

marginal = emmeans(glm.2a, ~ m_y)
pairs(marginal)

#best glm.1b for log(PD) - timepoint = more important here
glm.1b <- glmmTMB(log(PD) ~ m_y+(1|site:no_sp_lo), data=df.div.pd.pp, family=gaussian())
sim.resids <- simulateResiduals(glm.1b, plot = TRUE) #few points so idk seems ok visually

Anova(glm.1b, type = "III")
# Response: log(PD)
#               Chisq Df Pr(>Chisq)    
# (Intercept) 407.298  1  < 2.2e-16 ***
# m_y          12.505  1  0.0004059 ***

marginal = emmeans(glm.1b, ~ m_y)
pairs(marginal)

