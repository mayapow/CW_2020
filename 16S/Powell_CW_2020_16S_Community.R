#Maya Powell
#April, 2023 & 2024 for final pub stuff
#16S COMMUNITY COMPOSITION
#Based on scripts from many people online (see packages for reference if you want lol)
#and also some real live people: Ana Dulskiy (UNC Castillo lab), Nicola Kreifall (BU Davies Lab), Hannah Aichelman (BU Davies Lab) and Steph Smith (UNC Septer Lab)

# Setup
#set working directory
setwd("~/Documents/Castillo Lab/CW_2020/CW_2020_16S")
#load in libraries
library(rlang)
library(stringr)
library(dplyr)
library(stats)
library(ggplot2)
library(ggpubr)
library(vegan)
library(cowplot)
library(tidyverse)
#remotes::install_github("Jtrachsel/funfuns")
library("funfuns")
#BiocManager::install("microbiome")
#remotes::install_github("r-lib/rlang")
library(phyloseq)
library(microbiome)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(microbiomeutilities)
library(viridis)
#devtools::install_github("david-barnett/microViz")
library("microViz")

#ARE YOU COMING BACK TO THIS SCRIPT????
#CHECK LINE 112 AND START THERE!!!!!!!!!! PLEASE!!! SAVE YOURSELF!!!

## Read in data
load("CW_2020_16S_taxa2.Rdata")
#load in Phyloseq objects
ps.cleanest <- readRDS("CW_2020_16S_ps.cleanest.RDS")
ps.cleanest #27653 taxa and 144 samples, RAW DATA untrimmed, unrarefied
#remove duplicates from lane2 - they are very close (not stat sig dif) and lane 2 gets trimmed out earlier
ps.cleanest.nd = subset_samples(ps.cleanest, id!="N10" & id!="N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.cleanest.nd <- data.frame(sample_data(ps.cleanest.nd))

ps.less <- readRDS("CW_2020_16S_ps.less.RDS")
ps.less #27653 taxa and 135 samples, samples witN counts <1,000 removed (9 samples removed)
ps.less.nd = subset_samples(ps.less, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.less.nd <- data.frame(sample_data(ps.less.nd))

ps.rare <- readRDS("CW_2020_16S_ps.rare.RDS")
ps.rare #23759 taxa and 98 samples, rarefied to 10,000 (46 samples removed)
ps.rare.nd = subset_samples(ps.rare, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.rare.nd <- data.frame(sample_data(ps.rare.nd))

ps.trim <- readRDS("CW_2020_16S_ps.trim.RDS")
ps.trim #835 taxa and 134 samples, trimmed witN MCMC OTU for low abundance taxa (10 samples removed)
ps.trim.nd = subset_samples(ps.trim, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.trim.nd <- data.frame(sample_data(ps.trim.nd))

ps.trim.rare <- readRDS("CW_2020_16S_ps.trim.rare.RDS")
ps.trim.rare #835 taxa and 120 samples, trimmed dataset rarefied to 3498 (24 samples removed)
ps.trim.rare.nd = subset_samples(ps.trim.rare, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.trim.rare.nd <- data.frame(sample_data(ps.trim.rare.nd))

#Rename ASVs to be more informative
#go through and do this with all different ones
tax <- as.data.frame(ps.trim.rare.nd@tax_table@.Data)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
tax_table(ps.trim.rare.nd) <- as.matrix(tax.clean)
#do this for each one, then save and subset below, then only need to do step above once

#resaving and reading in nd dataframes with cleaned up tax tables
#can always redo from beginning
saveRDS(ps.cleanest.nd, "CW_2020_16S_ps.cleanest.nd.RDS")
saveRDS(ps.less.nd, "CW_2020_16S_ps.less.nd.RDS")
saveRDS(ps.rare.nd, "CW_2020_16S_ps.rare.nd.RDS")
saveRDS(ps.trim.nd, "CW_2020_16S_ps.trim.nd.RDS")
saveRDS(ps.trim.rare.nd, "CW_2020_16S_ps.trim.rare.nd.RDS")

#coming back to this script - START HERE and just read in phyloseq objects
#these have cleaned up tax tables and no duplicate samples
#and now these have all the environmental data added! MP 4/20/2023
ps.cleanest.nd <- readRDS("CW_2020_16S_ps.cleanest.nd.RDS")
ps.less.nd <- readRDS("CW_2020_16S_ps.less.nd.RDS")
ps.less.nd <- subset_samples(ps.less.nd,host_species!="astreoides")
ps.rare.nd <- readRDS("CW_2020_16S_ps.rare.nd.RDS")
ps.trim.nd <- readRDS("CW_2020_16S_ps.trim.nd.RDS")
ps.trim.rare.nd <- readRDS("CW_2020_16S_ps.trim.rare.nd.RDS")

#adding in site_nice column
ps.less.nd <- ps.less.nd %>% ps_mutate(site_nice =
                                         case_when(site == "SMB" ~ "Santa Martha Bay", 
                                                   site == "SMR" ~ "Santa Martha Reef",
                                                   site == "SWB" ~ "Spaanse Water Bay", 
                                                   site == "DB" ~ "Spaanse Water Reef"))
#add column for Habitat:Timepoint
ps.less.nd@sam_data$m_rb <- paste(ps.less.nd@sam_data$m_y,ps.less.nd@sam_data$reef_bay)

#View(data.frame(ps.less.nd@sam_data))
#add environmental data here too! and then resave - only doing this once and then don't need it anymore
head(sample_data(ps.trim.nd))

#split by species for each:
#do this for each type
ps.cl.ss <- subset_samples(ps.less.nd,host_species=="siderea")
ps.cl.sr <- subset_samples(ps.less.nd,host_species=="radians")
ps.cl.pp <- subset_samples(ps.less.nd,host_species=="porites")
#ps.cl.pa <- subset_samples(ps.trim.nd,host_species=="astreoides")

#make relative phyloseq objects
ps.all.rel <- transform_sample_counts(ps.less.nd, function(x) x / sum(x))
ps.ss.rel <- transform_sample_counts(ps.cl.ss, function(x) x / sum(x))
ps.sr.rel <- transform_sample_counts(ps.cl.sr, function(x) x / sum(x))
ps.pp.rel <- transform_sample_counts(ps.cl.pp, function(x) x / sum(x))
#ps.pa.rel <- transform_sample_counts(ps.cl.pa, function(x) x / sum(x))
#ps.all.rel <- subset_samples(ps.all.rel,host_species!="astreoides")

#JUST READ THESE IN BELOW IF ALREADY CREATED!
seq.all.ps.less <- data.frame(otu_table(ps.less.nd))
seq.ss.rel <- data.frame(otu_table(ps.ss.rel))
seq.sr.rel <- data.frame(otu_table(ps.sr.rel))
seq.pp.rel <- data.frame(otu_table(ps.pp.rel))
#seq.pa.rel <- data.frame(otu_table(ps.pa.rel))

#write.csv(seq.all.ps.rare, "seq.all.ps.rare.csv")
#seq.all.ps.rare.t <- as.data.frame(t(seq.all.ps.rare))
#seq.all.ps.rare.t <- seq.all.ps.rare.t %>%
#  mutate(ids=row.names(seq.all.ps.rare.t),
#         .before=A1)
#seq.all.ps.rare.t
#write.table(seq.all.ps.rare.t, file="seq.all.ps.rare.t.txt", sep="\t", row.names = FALSE)
#rownames(seq.all.ps.rare.t) <- NULL

saveRDS(seq.all.ps.less, "seq.all.rel.ps.less.RDS")
saveRDS(seq.ss.rel, "seq.ss.rel.ps.less.RDS")
saveRDS(seq.sr.rel, "seq.sr.rel.ps.less.RDS")
saveRDS(seq.pp.rel, "seq.pp.rel.ps.less.RDS")
#saveRDS(seq.pa.rel, "seq.pa.rel.ps.less.RDS")

#read in here
seq.all.rel <- readRDS("seq.all.rel.ps.less.RDS")
seq.ss.rel <- readRDS("seq.ss.rel.ps.less.RDS")
seq.sr.rel <- readRDS("seq.sr.rel.ps.less.RDS")
seq.pp.rel <- readRDS("seq.pp.rel.ps.less.RDS")
#seq.pa.rel <- readRDS("seq.pa.rel.ps.;ess.RDS")

#MP redoing things Aug 23rd 2023 to make them look good and maybe even the stats will be better?
##All
#ordination plots
samdf.all <- data.frame(sample_data(ps.all.rel))
samdf.all$time_zone <- paste(samdf.all$reef_bay,samdf.all$m_y)
row.names(samdf.trim) <- samdf.trim$id
seq.all.rel <- seq.all.rel[,!colnames(seq.all.rel) %in% core.all.ids ]
#remake phyloseq object - already relative
ps.acc <- phyloseq(otu_table(seq.all.rel, taxa_are_rows=FALSE), 
                   sample_data(samdf.trim), 
                   tax_table(taxa2))
ps.acc #27646 taxa accessory ps.cleanest, ps.trim 826 taxa
all.ord <- plot_ordination(ps.all.rel,ordinate(ps.all.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse()+
  theme_cowplot()
all.ord

#coral species
spp.ord <- plot_ordination(ps.all.rel,ordinate(ps.all.rel,"PCoA", "bray"),color="host_species", shape = "m_y")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  #ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites spp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))
spp.ord
ggsave(spp.ord, file="spp.ord.ps.less.pdf",w=10,h=8)

#ssid
ss.ord <- plot_ordination(ps.ss.rel,ordinate(ps.ss.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"),labels = c("Santa Martha Bay", "Santa Martha Reef", "Spaanse Water Bay", "Spaanse Water Reef"))
ss.ord
ggsave(ss.ord, file="ss.ord.ps.less.pdf",w=10,h=8)

#ssid timepoint
ss.ord <- plot_ordination(ps.ss.rel,ordinate(ps.ss.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual(name = "Timepoint", values = c("goldenrod","cornflowerblue", "skyblue"),labels = c("March 2020", "November 2020", "November 2021"))
ss.ord
ggsave(ss.ord, file="ss.ord.time.ps.less.pdf",w=10,h=8)

#srad
sr.ord <- plot_ordination(ps.sr.rel,ordinate(ps.sr.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  scale_color_manual(name = "Site", values = c("coral3","lightsalmon1"),labels = c("Santa Martha Bay", "Spaanse Water Bay"))
sr.ord
ggsave(sr.ord, file="sr.ord.ps.less.pdf",w=10,h=8)

#srad timepoint
sr.ord <- plot_ordination(ps.sr.rel,ordinate(ps.sr.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  scale_color_manual(name = "Timepoint", values = c("goldenrod","cornflowerblue", "skyblue"),labels = c("March 2020", "November 2020", "November 2021"))
sr.ord
ggsave(sr.ord, file="sr.ord.time.ps.less.pdf",w=10,h=8)

#ppor
pp.ord <- plot_ordination(ps.pp.rel,ordinate(ps.pp.rel,"PCoA", "bray"),color="site_zone")+
  #stat_ellipse(aes(linetype=reef_bay),linewidth=1)+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Porites spp."))))+
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise"),labels = c("Santa Martha Bay", "Santa Martha Reef"))
pp.ord
ggsave(pp.ord, file="pp.ord.ps.less.pdf",w=10,h=8)

#ppor timepoint
pp.ord <- plot_ordination(ps.pp.rel,ordinate(ps.pp.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=1)+
  theme_classic(base_size = 22)+
  ggtitle(substitute(paste(italic("Porites sp."))))+
  scale_color_manual(name = "Timepoint", values = c("cornflowerblue", "skyblue"),labels = c("November 2020", "November 2021"))
pp.ord
ggsave(pp.ord, file="pp.ord.time.ps.less.pdf",w=10,h=8)

pcoa.all <- ggarrange(ss.ord,sr.ord,pp.ord,nrow=1,ncol=3,legend="right",common.legend = TRUE, labels = c("A","B","C"), font.label = list(size = 30))
ggsave(pcoa.all, file="pcoa.all.png",width=18,height=5)

#now stats for above
#ALL
dist.all <- vegdist(seq.all.rel)
row.names(samdf.all)==row.names(seq.all.rel)
bet.all.spp <- betadisper(dist.all,samdf.all$host_species)

anova(bet.all.spp)
plot(bet.all.spp)

permutest(bet.all.spp, pairwise = TRUE, permutations = 999) 

adonis2(dist.all ~ host_species, data=samdf.all, permutations=999)
pairwise.adonis2(dist.all ~ host_species, data = samdf.all, permutations=999)
#all species different from eachother

#ssid
dist.ss <- vegdist(seq.ss.rel)
samdf.ss <- subset(samdf.all, host_species == "siderea")
row.names(samdf.ss)==row.names(seq.ss.rel)
bet.ss <- betadisper(dist.ss,samdf.ss$m_y_s_z)
anova(bet.ss) #ps.less ns
plot(bet.ss)

permutest(bet.ss, pairwise = TRUE, permutations = 999)

adonis2(dist.ss ~ reef_bay*m_y*site_zone,data=samdf.ss, permutations=999) #sig!

#go through each iteration of pairwise comparisons bc intxns don't work
pairwise.adonis2(dist.ss ~ m_y_s_z, data=samdf.ss, permutations=999)

#srad
dist.sr <- vegdist(seq.sr.rel)
samdf.sr <- subset(samdf.all, host_species == "siderea")
row.names(samdf.sr)==row.names(seq.sr.rel)
bet.sr <- betadisper(dist.sr,samdf.sr$m_y_s_z)
anova(bet.sr) #ps.lesr ns
plot(bet.sr)

permutest(bet.sr, pairwise = TRUE, permutations = 999)

adonis2(dist.sr ~ area*m_y,data=samdf.sr, permutations=999) #sig!

#go through each iteration of pairwise comparisons bc intxns don't work
pairwise.adonis2(dist.sr ~ m_y_s_z, data=samdf.sr, permutations=999)

#ppor
dist.pp <- vegdist(seq.pp.rel)
samdf.pp <- data.frame(sample_data(ps.pp.rel))
row.names(samdf.pp)==row.names(seq.pp.rel)

bet.pp <- betadisper(dist.pp,samdf.pp$m_y_s_z)
anova(bet.pp) #ns
plot(bet.pp) #overlap ns

#permutest(bet.pp, pairwise = TRUE, permutations = 999) #ns

adonis2(dist.pp ~ reef_bay*m_y, data=samdf.pp, permutations=999) #ns

pairwise.adonis2(dist.pp ~ m_y, data=samdf.pp, permutations=999)

# Bar plots
library(viridis)
#bars grouped by symbiont type

bar.core.rel.pp <- plot_bar(ps.pp.rel.z, fill="Genus")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Site")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#FFD92F","#1F78B4","#FC8D62","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.core.rel.pp
ggsave(bar.core.rel.pp, file="core.bar.ps.trim.pp.pdf",w=6,h=8)

#siderastrea siderea
ps_phy_glom_ss <- tax_glom(ps.ss.rel, "Phylum")
ps_fam_glom_ss <- tax_glom(ps.ss.rel, "Family")
ps.phy.sz <- merge_samples(ps_phy_glom_ss, "site_zone")
ps.phy.rel.sz <- transform_sample_counts(ps.phy.sz, function(x) x / sum(x))
gg.bar.ss <- plot_bar(ps.phy.rel.sz,fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=40)+
  #facet_wrap(~site_zone, scales = "free")+
  scale_fill_viridis(option="H",discrete = TRUE)+
  #scale_fill_brewer(palette = "RdYlBu")+
  ylab("Relative Abundance")+
  xlab("Site Zone")+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
#gg.bar.ss = gg.bar.ss + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.ss
ggsave(gg.bar.ss,file="bac.phy.barplot.ss.site.pdf",h=10,w=35)
#full plot - to put in supplementary, with all data across sites and time
gg.bar.ss <- plot_bar(ps_phy_glom_ss,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.ss
#each number corresponds to specific sample here
ggsave(gg.bar.ss,file="bac.phy.barplot.ss.site.season.sampleid.pdf",h=15,w=30)    

#siderastrea radians
ps_phy_glom_sr <- tax_glom(ps.sr.rel, "Phylum")
ps_fam_glom_sr <- tax_glom(ps.sr.rel, "Family")

# gg.bar.sr <- plot_bar(ps_phy_glom_sr,"sample_full",fill="Phylum")+
#   geom_bar(stat="identity")+
#   theme_classic(base_size=22)+
#   facet_wrap(~site_zone, scales = "free")+
#   ylab("Relative Abundance")
# gg.bar.sr = gg.bar.sr + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
# gg.bar.sr
# ggsave(gg.bar.sr,file="bac.phy.barplot.sr.site.pdf",h=5,w=15)
ps_phy_glom_sr <- tax_glom(ps.sr.rel, "Phylum")
ps_fam_glom_sr <- tax_glom(ps.sr.rel, "Family")
ps.phy.sz <- merge_samples(ps_phy_glom_sr, "site_zone")
ps.phy.rel.sz <- transform_sample_counts(ps.phy.sz, function(x) x / sum(x))
gg.bar.sr <- plot_bar(ps.phy.rel.sz,fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=35)+
  #facet_wrap(~site_zone, scales = "free")+
  scale_fill_viridis(option="H",discrete = TRUE)+
  #scale_fill_brewer(palette = "RdYlBu")+
  ylab("Relative Abundance")+
  xlab("Site Zone")+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
#gg.bar.sr = gg.bar.sr + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.sr
ggsave(gg.bar.sr,file="bac.phy.barplot.sr.site.pdf",h=10,w=30)

#full plot - to put in supplementary, with all data across sites and time
gg.bar.sr <- plot_bar(ps_phy_glom_sr,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.sr
#each number corresponds to specific sample here, removed duplicate of D6 = sample ID 91 Nov 2020 SWB
ggsave(gg.bar.sr,file="bac.phy.barplot.sr.site.season.sampleid.pdf",h=15,w=15)    

#porites porites
ps_phy_glom_pp <- tax_glom(ps.pp.rel, "Phylum")
ps_fam_glom_pp <- tax_glom(ps.pp.rel, "Family")

# gg.bar.pp <- plot_bar(ps_phy_glom_pp,"sample_full",fill="Phylum")+
#   geom_bar(stat="identity")+
#   theme_classic(base_size=22)+
#   facet_wrap(~site_zone, scales = "free")+
#   ylab("Relative Abundance")
# gg.bar.pp = gg.bar.pp + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
# gg.bar.pp
# ggsave(gg.bar.pp,file="bac.phy.barplot.pp.site.pdf",h=5,w=15)
ps_phy_glom_pp <- tax_glom(ps.pp.rel, "Phylum")
#ps_fam_glom_pp <- tax_glom(ps.pp.rel, "Family")
ps.phy.sz <- merge_samples(ps_phy_glom_pp, "site_zone")
ps.phy.rel.sz <- transform_sample_counts(ps.phy.sz, function(x) x / sum(x))
gg.bar.pp <- plot_bar(ps.phy.rel.sz,fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=35)+
  #facet_wrap(~site_zone, scales = "free")+
  scale_fill_viridis(option="H",discrete = TRUE)+
  #scale_fill_brewer(palette = "RdYlBu")+
  ylab("Relative Abundance")+
  xlab("Site Zone")+
  ggtitle(substitute(paste(italic("Porites porites"))))
#gg.bar.pp = gg.bar.pp + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.pp
ggsave(gg.bar.pp,file="bac.phy.barplot.pp.site.pdf",h=10,w=27)

twotogether <- ggarrange(gg.bar.sr,gg.bar.pp,nrow = 1, legend = "none", labels = c("B","C"),font.label = list(size = 50))
composition_all_phylum_ps.less <- ggarrange(gg.bar.ss,twotogether,nrow = 1, common.legend = TRUE, legend = "bottom",labels = c("A"), font.label = list(size = 50))
ggsave(composition_all_phylum_ps.less, file="composition_all_phylum_ps.less.pdf",width=30,height=15)

#full plot - to put in supplementary, with all data acropp sites and time
gg.bar.pp <- plot_bar(ps_phy_glom_pp,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.pp
#each number corresponds to specific sample here #removed duplicates same as bacterial dataset (Ns)
ggsave(gg.bar.pp,file="bac.phy.barplot.pp.site.season.sampleid.pdf",h=15,w=15)    

gg.bar.pp <- plot_bar(ps.cl.pp,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.pp

#trying microshades
#remotes::install_github("KarstensLab/microshades")
library(microshades)
#remotes::install_github("mikemc/speedyseq")
library(speedyseq)

taxa.less <- as.data.frame(ps.less.nd@tax_table)

#all
#porites porites
ps_phy_glom_all <- tax_glom(ps.less.nd, "Phylum")
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
mdf_prep_all <- prep_mdf(ps.less.nd, subgroup_level = "Family")

# Create a color object for the specified data
color_obj_all <- create_color_dfs(mdf_prep_all, selected_groups = c("Proteobacteria","Bacteroidota","Desulfobacterota","Planctomycetota","Firmicutes"), group_level = "Phylum", subgroup_level = "Family", cvd = TRUE)
#can only do up to 5 based on the palette number
#using top 5 based on colors here
sums <- as.data.frame(taxa_sums(ps_phy_glom_all))
View(data.frame(ps.less.nd@tax_table))
#MAKE SURE TO SORT BY MOST ABUNDANT
#sq2 Proteobacteria
#sq16 Bacteroidota
#sq22 Desulfobacterota
#sq132 Planctomycetota
#sq35 Firmicutes
#sq47 Cyanobacteria

#most abundant taxa for all overall determined by above are:

# Extract
mdf_all <- color_obj_all$mdf
cdf_all <- color_obj_all$cdf

#plot all - species comparisons
ms_bar_plot_all <- plot_microshades(mdf_all, cdf_all, group_label = "Phylum - Family")
ms_bar_plot_final <- ms_bar_plot_all + scale_y_continuous(labels = scales::percent, expand = expansion(0))+
  theme_classic(base_size=30)+
  theme(legend.position = "bottom", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(host_species~.,scales = "free")+
  #labeller=labeller(site_zone = c("Santa Martha Bay", "Santa Martha Reef", "Spaanse Water Bay", "Spaanse Water Reef")))+
  xlab("Coral Species")+
  ylab("Relative Abundance")
  #ggtitle(substitute(paste(italic("All Species"))))
ms_bar_plot_final
ggsave(ms_bar_plot_final, file = "ms.all.spp.ps.less.pdf", w=30, h=10)

#subset different groups out
mdf_ss <- mdf_all %>% filter(host_species == "siderea")
mdf_sr <- mdf_all %>% filter(host_species == "radians")
mdf_pp <- mdf_all %>% filter(host_species == "porites")

#porites porites

#plot pp
ms_bar_plot_pp <- plot_microshades(mdf_pp, cdf_all, group_label = "Phylum - Family")
ms_bar_plot_pp_final <- ms_bar_plot_pp + scale_y_continuous(labels = scales::percent, expand = expansion(0))+
  theme_classic(base_size=22)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(.~m_y,scales = "free")+
  xlab("Site")+
  ylab("Relative Abundance")+
  ggtitle(substitute(paste(italic("Porites sp."))))
ms_bar_plot_pp_final
ggsave(ms_bar_plot_pp_final, file = "ms.pp.ps.less.pdf", w=30, h=10)
ggsave(ms_bar_plot_pp_final, file = "ms.pp.time.ps.less.pdf", w=30, h=10)

#siderastrea siderea
#plot ss
ms_bar_plot_ss <- plot_microshades(mdf_ss, cdf_all, group_label = "Phylum - Family")
ms_bar_plot_ss_final <- ms_bar_plot_ss + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme_classic(base_size=22)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(.~m_y,scales = "free")+
  xlab("Site")+
  ylab("Relative Abundance")+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ms_bar_plot_ss_final
ggsave(ms_bar_plot_ss_final, file = "ms.ss.ps.less.pdf", w=30, h=15)
ggsave(ms_bar_plot_ss_final, file = "ms.ss.time.ps.less.pdf", w=30, h=10)

#siderastrea radians
#plot sr
ms_bar_plot_sr <- plot_microshades(mdf_sr, cdf_all, group_label = "Phylum - Family")
ms_bar_plot_sr_final <- ms_bar_plot_sr + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme_classic(base_size=22)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  facet_wrap(.~m_y,scales = "free")+
  xlab("Site")+
  ylab("Relative Abundance")+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
ms_bar_plot_sr_final
ggsave(ms_bar_plot_sr_final, file = "ms.sr.ps.less.pdf", w=30, h=10)
ggsave(ms_bar_plot_sr_final, file = "ms.sr.time.ps.less.pdf", w=30, h=10)

#saving plots
ggsave(ms_bar_plot_sr_final, file = "microshades_bar_sr.ps.less.png",h=5,w=20)
ggsave(ms_bar_plot_ss_final, file = "microshades_bar_ss.ps.less.png",h=10,w=20)
ggsave(ms_bar_plot_pp_final, file = "microshades_bar_pp.ps.less.png",h=5,w=20)

#putting plots together
bac_pp_sr_plot <- ggarrange(ms_bar_plot_sr_final,ms_bar_plot_pp_final, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 50))
bac_all_plot <- ggarrange(ms_bar_plot_ss_final,bac_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 50))
ggsave(bac_all_plot, file="microshades_all_bar.ps.less.pdf",width=30,height=12)

#putting time plots together
bac_all_time_plot <- ggarrange(ms_bar_plot_ss_final,ms_bar_plot_sr_final,ms_bar_plot_pp_final, nrow=3, ncol = 1, legend = "none", common.legend = TRUE, labels = c("A","B","C"),font.label = list(size = 50))
ggsave(bac_all_time_plot, file="microshades_all_bar_time.ps.less.pdf",width=15,height=20)

#just ss and pp for emes poster april 2nd 2024
ms_ss_pp <- ggarrange(ms_bar_plot_ss_final,ms_bar_plot_pp_final,nrow = 2, ncol=1, heights= c(2,1), common.legend = TRUE, legend = "right", labels = c("A","B"),font.label = list(size = 50))
ggsave(ms_ss_pp, file="microshades_ss_pp_bar_ps.less.png",width=20,height=15)

## All, by Phylum
#make category of all less than 10%
BiocManager::install("miaViz")
library("miaViz")
tse.trim.nd <- makeTreeSummarizedExperimentFromPhyloseq(ps.trim.nd) 
tse.trim.nd
tse.trim.nd.rel <- transformCounts(tse.trim.nd, method = "relabundance")

plotAbundanceDensity(tse.trim.nd.rel, layout = "density", assay_name = "relabundance",
                     n = 5, colour_by="site_zone", point_alpha=1/10) +
  scale_x_log10()

ps.ten.all <- core(ps.trim.nd, detection = 0.1, prevalence = 0)
core.all.tax <- data.frame(ps.core.all@tax_table)

ps_phy_glom_all <- tax_glom(ps.trim.nd, "Phylum")
ps.phy.all.s <- merge_samples(ps_phy_glom_all, "m_y_s_z")
ps.rel.phy.all.s <- transform_sample_counts(ps.phy.all.s, function(x) x / sum(x))
bar.all.rel.all <- plot_bar(ps.rel.phy.all.s, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Site and Timepoint")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#9ECAE1","#E78AC3","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.all.rel.all
ggsave(bar.all.rel.all, file="all.bar.phylum.ps.trim.all.site.time.pdf",h=6,w=10)

#ss by phylum
ps_phy_glom_ss <- tax_glom(ps.cl.ss, "Phylum")
ps.phy.ss.s <- merge_samples(ps_phy_glom_ss, "m_y_s_z")
ps.rel.phy.ss.s <- transform_sample_counts(ps.phy.ss.s, function(x) x / sum(x))
bar.all.rel.ss <- plot_bar(ps.rel.phy.ss.s, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Site and Timepoint")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#9ECAE1","#E78AC3","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.all.rel.ss
ggsave(bar.all.rel.ss, file="all.bar.ps.trim.ss.phylum.site.time.pdf",w=10,h=6)

##sr by phylum
ps_phy_glom_sr <- tax_glom(ps.cl.sr, "Phylum")
ps.phy.sr.s <- merge_samples(ps_phy_glom_sr, "m_y_s_z")
ps.rel.phy.sr.s <- transform_sample_counts(ps.phy.sr.s, function(x) x / sum(x))
bar.all.rel.sr <- plot_bar(ps.rel.phy.sr.s, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Site and Timepoint")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#9ECAE1","#E78AC3","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.all.rel.sr
ggsave(bar.all.rel.sr, file="all.bar.ps.trim.sr.phylum.site.time.pdf",w=10,h=8)

##pp by phylum
ps_phy_glom_pp <- tax_glom(ps.cl.pp, "Phylum")
ps.phy.pp.s <- merge_samples(ps_phy_glom_pp, "m_y_s_z")
ps.rel.phy.pp.s <- transform_sample_counts(ps.phy.pp.s, function(x) x / sum(x))
bar.all.rel.pp <- plot_bar(ps.rel.phy.pp.s, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Site and Timepoint")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#9ECAE1","#E78AC3","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.all.rel.pp
ggsave(bar.all.rel.pp, file="all.bar.ps.trim.pp.phylum.site.time.pdf",w=10,h=8)

ps.sz <- merge_samples(ps.trim.nd, "site_zone")
ps.rel.sz <- transform_sample_counts(ps.sz, function(x) x / sum(x))

ps.all.tab <- psmelt(ps.rel.sz)%>%
  filter(!is.na(Abundance))%>%
  group_by(site_zone,m_y,m_y_s_z,Phylum,OTU)%>%
  summarize_at("Abundance",mean)

gg.bar.all.phy <- ggplot(ps.all.tab,aes(x=site_zone,y=Abundance,fill=Phylum))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  #theme(legend.position="none")+
  #scale_fill_brewer(palette="Set2")+
  xlab('Reef zone')+
  facet_wrap(m_y~site_zone,scales = "free_x")
gg.bar.all.phy
#ggsave(gg.bar.all,file="bac.bar.all.pdf",height=4)

#heat map more of them they are maybe different than the ones before and the ones after
BiocManager::install(ComplexHeatmap), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)

#FOR EMMA - share with her to see if there are any diazos relevant
#Make PCA with taxa arrows
#ssid
ps.cl.ss.pca <- ps.cl.ss %>% # replace ps.cl.ss with your phyloseq object
  tax_transform("clr", rank = "Phylum") %>% # change here if you want to plot by a different taxonomic rank - eg. Genus, Class, etc.
  ord_calc()
#do this twice if you want to split it by apo and sym - make one for apo one for sym
ss.pca <- ps.cl.ss.pca  %>% # replace ps.cl.ss.pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by - yours will be treatment
           plot_taxa = 1:64, # number of taxa to label on plot, I change based on how many there are (you can check in your tax table from your phyloseq object)
           size=3.5, # size of data points
           tax_vec_length = 1.7, # length of arrows corresponding to taxa vectors
           tax_lab_length = 1.8, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  ggtitle(label = "Siderastrea siderea PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-5, 5), ylim = c(-5, 5)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
ss.pca
ggsave(ss.pca, file="ss.phylum.PCA.ps.trim.pdf", height=8, width=8)

#srad
ps.cl.sr.pca <- ps.cl.sr %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Phylum") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
sr.pca <- ps.cl.sr.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 2, # length of arrows corresponding to taxa vectors
           tax_lab_length = 2.1, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="site_zone",values=c("#E78AC3","#FC8D62"))+
  ggtitle(label = "Siderastrea radians PCA") + # obv change to whatever you'd like for the title of your plot
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-4, 4), ylim = c(-4, 4)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
sr.pca
ggsave(sr.pca, file="sr.phylum.PCA.ps.trim.pdf", height=8, width=8)

#ppor
ps.cl.pp.pca <- ps.cl.pp %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Phylum") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
pp.pca <- ps.cl.pp.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 2, # length of arrows corresponding to taxa vectors
           tax_lab_length = 2.1, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB"))+
  ggtitle(label = "Porites porites PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-7, 7), ylim = c(-6, 6)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
pp.pca
ggsave(pp.pca, file="pp.phylum.PCA.ps.trim.pdf", height=8, width=8)

#past
ps.cl.pa.pca <- ps.cl.pa %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Phylum") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
pa.pca <- ps.cl.pa.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 2, # length of arrows corresponding to taxa vectors
           tax_lab_length = 2.1, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  ggtitle(label = "Porites astreoides PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-7, 8), ylim = c(-4, 5)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
pa.pca
#dev.off()
ggsave(pa.pca, file="pa.phylum.PCA.ps.cleanest.pdf", height=8, width=10)

#testing with adding taxa arrows with 0.05


#now do core PCA bc ppl still care about the core ig
#ssid
ps.core.ss.pca <- ps.core.ss %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Genus") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
ss.core.pca <- ps.core.ss.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 0.5, # length of arrows corresponding to taxa vectors
           tax_lab_length = 0.6, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  ggtitle(label = "Siderastrea siderea Core PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3, 3), ylim = c(-3, 3)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
ss.core.pca
ggsave(ss.core.pca, file="ss.core.PCA.ps.trim.pdf", height=8, width=8)

#srad core
ps.core.sr.pca <- ps.core.sr %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Genus") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
sr.core.pca <- ps.core.sr.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 0.5, # length of arrows corresponding to taxa vectors
           tax_lab_length = 0.6, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#FC8D62"))+
  ggtitle(label = "Siderastrea radians Core PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3, 3), ylim = c(-3, 3)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
sr.core.pca
ggsave(sr.core.pca, file="sr.core.PCA.ps.trim.pdf", height=8, width=8)

#ppor
ps.core.pp.pca <- ps.core.pp %>% # replace es.hsp with your phyloseq object
  tax_transform("clr", rank = "Genus") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
pp.core.pca <- ps.core.pp.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "site_zone", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 0.5, # length of arrows corresponding to taxa vectors
           tax_lab_length = 0.6, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB"))+
  ggtitle(label = "Porites porites Core PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=site_zone)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3, 3), ylim = c(-3, 3)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
pp.core.pca
ggsave(pp.core.pca, file="pp.core.PCA.ps.trim.pdf", height=8, width=8)

#mean for environmental drivers! does this actually mean anything? I think an evntl correlation plot of some sort might be better!
#ssid core
ps.cl.ss <- ps.cl.ss %>% tax_fix(unknowns = c("Unknown Family"))
bray.dists.ss.cl <- ps.cl.ss %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("bray")
perm2 <- bray.dists.ss.cl %>%
  dist_permanova(variables = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site"), seed = 321)
perm.mean.env.cl.ss <- bray.dists.ss.cl %>%
  ord_calc()%>%#constraints = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.mean.env.cl.ss
ggsave(perm.mean.env.cl.ss, file="CAP.mean.env.cl.ss.ps.trim.pdf")

#now look at range but like how can we do them together???
perm2 <- bray.dists.ss.cl %>%
  dist_permanova(variables = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site"), seed = 321)
perm.range.env.cl.ss <- perm2 %>%
  ord_calc(constraints = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.range.env.cl.ss
ggsave(perm.range.env.cl.ss, file="CAP.range.env.cl.ss.ps.trim.pdf")

#srad
ps.cl.sr <- ps.cl.sr %>% tax_fix(unknowns = c("Unknown Family"))
bray.dists.sr.cl <- ps.cl.sr %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("bray")
perm2 <- bray.dists.sr.cl %>%
  dist_permanova(variables = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site"), seed = 321)
perm.mean.env.cl.sr <- perm2 %>%
  ord_calc(constraints = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#FC8D62"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.mean.env.cl.sr
ggsave(perm.mean.env.cl.sr, file="CAP.mean.env.cl.sr.ps.trim.pdf")

#now look at range
perm2 <- bray.dists.sr.cl %>%
  dist_permanova(variables = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site"), seed = 321)
perm.range.env.cl.sr <- perm2 %>%
  ord_calc(constraints = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#FC8D62"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.range.env.cl.sr
ggsave(perm.range.env.cl.sr, file="CAP.range.env.cl.sr.ps.trim.pdf")

#ppor
ps.cl.pp <- ps.cl.pp %>% tax_fix(unknowns = c("Unknown Family"))
bray.dists.pp.cl <- ps.cl.pp %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("bray")
perm2 <- bray.dists.pp.cl %>%
  dist_permanova(variables = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site"), seed = 321)
perm.mean.env.cl.pp <- perm2 %>%
  ord_calc(constraints = c("mean_pH_site","mean_do_site","mean_par_site","mean_temp_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#8DA0CB"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.mean.env.cl.pp
ggsave(perm.mean.env.cl.pp, file="CAP.mean.env.cl.pp.ps.trim.pdf")

#now look at range
perm2 <- bray.dists.pp.cl %>%
  dist_permanova(variables = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site"), seed = 321)
perm.range.env.cl.pp <- perm2 %>%
  ord_calc(constraints = c("range_pH_site","range_temp_site", "range_sal_site","range_do_site","range_par_site")) %>%
  ord_plot(
    colour = "site_zone", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(
      type = "text", max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = site_zone), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(values=c("#E78AC3","#8DA0CB"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
perm.range.env.cl.pp
ggsave(perm.range.env.cl.pp, file="CAP.range.env.cl.pp.ps.trim.pdf")

# PCOA plots

## Plots - site {.tabset}

### Trimmed
all.ord.all <- plot_ordination(ps.all.rel,ordinate(ps.all.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse()+
  geom_point(size = 3)+
  theme_classic(base_size = 22)+
  scale_color_manual(values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))
all.ord.all
#scale_color_manual(name="Site",values=c("darkslategray3","darkslategray4","#000004"))+
#scale_shape_manual(name="Site",values=c(8,4,9))#+
#xlab("Axis 1 (42.6%)")+
#ylab("Axis 2 (24.1%)")+
#ggtitle("Rarefied")
ggsave(all.ord.all, file="all.ord.all.ps.trim.pdf",w=8,h=8)

#ssid
ss.ord.all <- plot_ordination(ps.cl.ss,ordinate(ps.cl.ss,"PCoA", "bray"),color="m_y")+
  stat_ellipse()+
  geom_point(size = 3)+
  theme_classic(base_size = 22)#+
  #scale_color_manual(values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))
ss.ord.all
ggsave(ss.ord.all, file="ss.ord.all.ps.trim.pdf",w=8,h=8)

#srad
sr.ord.all <- plot_ordination(ps.cl.sr,ordinate(ps.cl.sr,"PCoA", "bray"),color="m_y")+
  stat_ellipse()+
  geom_point(size = 3)+
  theme_classic(base_size = 22)#+
  #scale_color_manual(values=c("#E78AC3","#FC8D62"))
sr.ord.all
ggsave(sr.ord.all, file="sr.ord.all.ps.trim.pdf",w=8,h=8)

#ppor
pp.ord.all <- plot_ordination(ps.cl.pp,ordinate(ps.cl.pp,"PCoA", "bray"),color="m_y")+
  stat_ellipse()+
  geom_point(size = 3)+
  theme_classic(base_size = 22)#+
  #scale_color_manual(values=c("#E78AC3","#8DA0CB"))
pp.ord.all
ggsave(pp.ord.all, file="pp.ord.all.ps.trim.pdf",w=8,h=8)

#past
pa.ord.all <- plot_ordination(ps.cl.pa,ordinate(ps.cl.pa,"PCoA", "bray"),color="site_zone")+
  stat_ellipse()+
  geom_point(size = 3)+
  theme_classic(base_size = 22)+
  scale_color_manual(values=c("#FC8D62","#66C2A5"))
pa.ord.all
ggsave(pa.ord.all, file="pa.ord.all.ps.trim.pdf",w=8,h=8)
ggarrange(ss.ord.all,sr.ord.all,pp.ord.all,nrow=1,common.legend=T,legend="right")
#ggsave("pcoa.all.zone.pdf",width=11,height=3)

#Stats

#Help on adonis (here)[https://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html#more]
                      
#stats  - ps.trim, all taxa
dist.trim <- vegdist(seq.all.rel)
samdf.trim <- data.frame(sample_data(ps.all.rel))
row.names(samdf.trim)==row.names(seq.all.rel)
bet.all <- betadisper(dist.trim,samdf.trim$site_zone)
anova(bet.all) #ns
plot(bet.all) #not sig
permutest(bet.all, pairwise = TRUE, permutations = 999)
#ps.trim all ns 
adonis2(dist.trim ~ site_zone, data=samdf.trim, permutations=999) 
#ps.trim***
adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations=999)
adonis2(dist.trim ~ site_zone, strata=samdf.trim$m_y, data=samdf.trim, permutations=999)
pairwise.adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations = 999)
#ps.trim both ***

#ssid
dist.trim <- vegdist(seq.ss.rel)
samdf.trim <- data.frame(sample_data(ps.cl.ss))
row.names(samdf.trim)==row.names(seq.ss.rel)
bet.ss <- betadisper(dist.trim,samdf.trim$m_y)
anova(bet.ss) #ns (m_y very sig)
plot(bet.ss) #not sig
permutest(bet.ss, pairwise = TRUE, permutations = 999)
#ps.trim ss ns 
adonis2(dist.trim ~ m_y, data=samdf.trim, permutations=999) 
#ps.trim***
adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations=999)
adonis2(dist.trim ~ site_zone, strata=samdf.trim$m_y, data=samdf.trim, permutations=999)
pairwise.adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations = 999)
#ps.trim all either ** or ***

#srad
dist.trim <- vegdist(seq.sr.rel)
samdf.trim <- data.frame(sample_data(ps.cl.sr))
row.names(samdf.trim)==row.names(seq.sr.rel)
bet.sr <- betadisper(dist.trim,samdf.trim$m_y)
anova(bet.sr) #** (ns for time)
plot(bet.sr) #very little overlap! sig!
permutest(bet.sr, pairwise = TRUE, permutations = 999)
#ps.trim sr ** 
adonis2(dist.trim ~ site_zone, data=samdf.trim, permutations=999) 
#ps.trim sr ***
adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations=999)
adonis2(dist.trim ~ site_zone, strata=samdf.trim$m_y, data=samdf.trim, permutations=999)
pairwise.adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations = 999)
#ps.trim sr all either ***

#ppor
dist.trim <- vegdist(seq.pp.rel)
samdf.trim <- data.frame(sample_data(ps.cl.pp))
row.names(samdf.trim)==row.names(seq.pp.rel)
bet.pp <- betadisper(dist.trim,samdf.trim$site_zone)
anova(bet.pp) #ns
plot(bet.pp) #not sig
permutest(bet.pp, pairwise = TRUE, permutations = 999)
#ps.trim pp ns 
adonis2(dist.trim ~ site_zone, data=samdf.trim, permutations=999) 
#ps.trim ns
adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations=999)
adonis2(dist.trim ~ site_zone, strata=samdf.trim$m_y, data=samdf.trim, permutations=999)
pairwise.adonis2(dist.trim ~ site_zone/m_y, data=samdf.trim, permutations = 999)
#ps.trim ns, site_zone:m_y p=0.065 close

#stats: ps.cleanest
#stats 
dist.cleanest <- vegdist(seq.all.rel)
samdf.cleanest <- data.frame(sample_data(ps.all.rel))
row.names(samdf.cleanest)==row.names(seq.all.rel)
bet.all <- betadisper(dist.cleanest,samdf.cleanest$site_zone)
anova(bet.all) #ps.cleanest *
plot(bet.all) #not sig
permutest(bet.all, pairwise = TRUE, permutations = 999)
#ps.cleanest SWB vs SWR
adonis2(dist.cleanest ~ site_zone, data=samdf.cleanest, permutations=999) 
#ps.cleanest***
adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations=999)
adonis2(dist.cleanest ~ site_zone, strata=samdf.cleanest$m_y, data=samdf.cleanest, permutations=999)
pairwise.adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations = 999)
#ps.cleanest both ***

#ssid
dist.cleanest <- vegdist(seq.ss.rel)
samdf.cleanest <- data.frame(sample_data(ps.cl.ss))
row.names(samdf.cleanest)==row.names(seq.ss.rel)
bet.ss <- betadisper(dist.cleanest,samdf.cleanest$site_zone)
anova(bet.ss) #ps.cleanest *
plot(bet.ss) #doesn't look sig but is barely
permutest(bet.ss, pairwise = TRUE, permutations = 999)
#ps.cleanest ss SMR vs SWR
adonis2(dist.cleanest ~ site_zone, data=samdf.cleanest, permutations=999) 
#ps.cleanest***
adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations=999)
adonis2(dist.cleanest ~ site_zone, strata=samdf.cleanest$m_y, data=samdf.cleanest, permutations=999)
pairwise.adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations = 999)
#ps.cleanest all either ** or ***

#srad
dist.cleanest <- vegdist(seq.sr.rel)
samdf.cleanest <- data.frame(sample_data(ps.cl.sr))
row.names(samdf.cleanest)==row.names(seq.sr.rel)
bet.sr <- betadisper(dist.cleanest,samdf.cleanest$site_zone)
anova(bet.sr) #**
plot(bet.sr) #very little overlap! sig!
permutest(bet.sr, pairwise = TRUE, permutations = 999)
#ps.cleanest sr ** 
adonis2(dist.cleanest ~ site_zone, data=samdf.cleanest, permutations=999) 
#ps.cleanest sr ***
adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations=999)
adonis2(dist.cleanest ~ site_zone, strata=samdf.cleanest$m_y, data=samdf.cleanest, permutations=999)
pairwise.adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations = 999)
#ps.cleanest sr ***

#ppor
dist.cleanest <- vegdist(seq.pp.rel)
samdf.cleanest <- data.frame(sample_data(ps.cl.pp))
row.names(samdf.cleanest)==row.names(seq.pp.rel)
bet.pp <- betadisper(dist.cleanest,samdf.cleanest$site_zone)
anova(bet.pp) #ns
plot(bet.pp) #not sig
permutest(bet.pp, pairwise = TRUE, permutations = 999)
#ps.cleanest pp ns 
adonis2(dist.cleanest ~ site_zone, data=samdf.cleanest, permutations=999) 
#ps.cleanest ns
adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations=999)
adonis2(dist.cleanest ~ site_zone, strata=samdf.cleanest$m_y, data=samdf.cleanest, permutations=999)
pairwise.adonis2(dist.cleanest ~ site_zone/m_y, data=samdf.cleanest, permutations = 999)
#ps.cleanest ns, site_zone:m_y p=0.068 close

#stats all rare
#stats 
dist.rare <- vegdist(seq.all.rel)
samdf.rare <- data.frame(sample_data(ps.all.rel))
row.names(samdf.rare)==row.names(seq.all.rel)
bet.all <- betadisper(dist.rare,samdf.rare$site_zone)
anova(bet.all) #ns
plot(bet.all) #not sig
permutest(bet.all, pairwise = TRUE, permutations = 999)
#ps.rare all ns 
adonis2(dist.rare ~ site_zone, data=samdf.rare, permutations=999) 
#ps.rare***
adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations=999)
adonis2(dist.rare ~ site_zone, strata=samdf.rare$m_y, data=samdf.rare, permutations=999)
pairwise.adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations = 999)
#ps.rare *** or **

#ssid
dist.rare <- vegdist(seq.ss.rel)
samdf.rare <- data.frame(sample_data(ps.cl.ss))
row.names(samdf.rare)==row.names(seq.ss.rel)
bet.ss <- betadisper(dist.rare,samdf.rare$site_zone)
anova(bet.ss) #ns
plot(bet.ss) #not sig
permutest(bet.ss, pairwise = TRUE, permutations = 999)
#ps.rare ss ns 
adonis2(dist.rare ~ site_zone, data=samdf.rare, permutations=999) 
#ps.rare***
adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations=999)
adonis2(dist.rare ~ site_zone, strata=samdf.rare$m_y, data=samdf.rare, permutations=999)
pairwise.adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations = 999)
#ps.rare *, ** or *** depending on comparisons

#srad
dist.rare <- vegdist(seq.sr.rel)
samdf.rare <- data.frame(sample_data(ps.cl.sr))
row.names(samdf.rare)==row.names(seq.sr.rel)
bet.sr <- betadisper(dist.rare,samdf.rare$site_zone)
anova(bet.sr) #**
plot(bet.sr) #very little to no overlap! sig!
permutest(bet.sr, pairwise = TRUE, permutations = 999)
#ps.rare sr ** 
adonis2(dist.rare ~ site_zone, data=samdf.rare, permutations=999) 
#ps.rare sr ***
adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations=999)
adonis2(dist.rare ~ site_zone, strata=samdf.rare$m_y, data=samdf.rare, permutations=999)
pairwise.adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations = 999)
#ps.rare sr site *** site over time *

#ppor
dist.rare <- vegdist(seq.pp.rel)
samdf.rare <- data.frame(sample_data(ps.cl.pp))
row.names(samdf.rare)==row.names(seq.pp.rel)
bet.pp <- betadisper(dist.rare,samdf.rare$site_zone)
anova(bet.pp) #ns
plot(bet.pp) #not sig
permutest(bet.pp, pairwise = TRUE, permutations = 999)
#ps.rare pp ns 
adonis2(dist.rare ~ site_zone, data=samdf.rare, permutations=999) 
#ps.rare ns
adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations=999)
adonis2(dist.rare ~ site_zone, strata=samdf.rare$m_y, data=samdf.rare, permutations=999)
pairwise.adonis2(dist.rare ~ site_zone/m_y, data=samdf.rare, permutations = 999)
#ps.rare ns

#look at individual samples:
#SSID
#majority D symbionts:
ps.ss.sample <- subset_samples(ps.cl.ss, no_sp_lo=="96_siderea_SW_reef")
ps.ss.s.rel <- transform_sample_counts(ps.ss.sample, function(x) x / sum(x))

ps_glom.ss <- tax_glom(ps.cl.ss, "Phylum")
ps.ss.z <- merge_samples(ps_glom.ss, "m_y_s_z")
ps.ss.rel.z <- transform_sample_counts(ps.ss.z, function(x) x / sum(x))

gg.ss.bar.s <- plot_bar(ps.ss.rel.z,x="m_y_s_z",y="Abundance", fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Sample ID")+
  ylab("Relative Abundance")
gg.ss.bar.s
ggsave(gg.ss.bar.s,file="gg.ss.bar.97SSSWR.pdf",h=10,w=14) 


ps_glom.ss <- tax_glom(ps.cl.ss, "Phylum")
ps.ss.z <- merge_samples(ps_glom.ss, "m_y_s_z")
ps.ss.rel.z <- transform_sample_counts(ps.ss.z, function(x) x / sum(x))

gg.ss.bar.s <- plot_bar(ps.ss.rel.z, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Sample ID")+
  ylab("Relative Abundance")
gg.ss.bar.s
ggsave(gg.ss.bar.s,file="bar_phylum_ss_site_time_ps.trim.pdf",h=10,w=14) 


# ANCOM
#I COULDN'T GET THIS TO WORK BC IT'S OUTDATED
#this is from Ana's code and the ancombc method from microbiome marker package is better and easier and it does it for you

#Tutorial [here](https://github.com/FrederickHuangLin/ANCOM)

#maya fun tests!!
       
#heat maps               
devtools::install_github("microsud/microbiomeutilities")
#remotes::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(viridis)
devtools::install_github("vmikk/metagMisc")
library(metagMisc)

#make subsets of top 10 and top 20 taxa
#ps.tr.ss <- subset_samples(ps.trim.nd,host_species=="siderea")
#ps.cl.ss <- subset_samples(ps.cleanest.nd,host_species=="siderea")

ps.cl.ss.10 <- phyloseq_filter_top_taxa(ps.cl.ss, n = 10)
ps.cl.ss.20 <- phyloseq_filter_top_taxa(ps.cl.ss, n = 20)
ps.cl.sr.10 <- phyloseq_filter_top_taxa(ps.cl.sr, n = 10)
ps.cl.sr.20 <- phyloseq_filter_top_taxa(ps.cl.sr, n = 20)
ps.cl.pp.10 <- phyloseq_filter_top_taxa(ps.cl.pp, n = 10)
ps.cl.pp.20 <- phyloseq_filter_top_taxa(ps.cl.pp, n = 20)

# ps.tr.ss.10 <- phyloseq_filter_top_taxa(ps.tr.ss, n = 10)
# ps.tr.ss.20 <- phyloseq_filter_top_taxa(ps.tr.ss, n = 20)

#heatmaps of top 10 and top 20 taxa
#you can edit the colors if you want!
#heat map info: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_heatmap
#colors: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
#your sample.order will be treatment @emma
#siderea
plot_ly(z = data, colorscale = "Greys", type = "heatmap")  
ss.10.heat <- plot_heatmap(ps.cl.ss.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Genus",
             sample.order="site_nice", taxa.order = "Genus", low = "white",
             high = "red", na.value = "white", title = "Siderastrea Siderea Top 10 Genera")
ss.10.heat <- ss.10.heat + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
ss.10.heat
#radians
sr.10.heat <- plot_heatmap(ps.cl.sr.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Genus",
                           sample.order="site_nice", taxa.order = "Genus", low = "khaki",
                           high = "red", na.value = "white", title = "Siderastrea radians Top 10 Genera")
sr.10.heat <- sr.10.heat + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
sr.10.heat

#porites
pp.10.heat <- plot_heatmap(ps.cl.pp.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Genus",
                           sample.order="site_nice", taxa.order = "Genus", low = "khaki",
                           high = "red", na.value = "white", title = "Porites porites Top 10 Genera")
pp.10.heat <- pp.10.heat + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Porites spp."))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
pp.10.heat

pp.sr.heat <- ggarrange(sr.10.heat, pp.10.heat, nrow=1, ncol = 2, legend = "none", labels = c("B","C"),font.label = list(size = 40))
heat.all.top.10.genus.ps.less <- ggarrange(ss.10.heat,pp.sr.heat, nrow=2, ncol = 1, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
ggsave(heat.all.top.10.genus.ps.less, file="heatmap.all.top.10.genus.ps.less.png",width=20,height=10)

#also looking at this at the phylum level too
#siderea
ss.10.heat.phy <- plot_heatmap(ps.cl.ss.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Phylum",
                           sample.order="site_nice", taxa.order = "Phylum", low = "khaki",
                           high = "red", na.value = "white", title = "Siderastrea Siderea Top 10 Genera")
ss.10.heat.phy <- ss.10.heat.phy + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
ss.10.heat.phy
#radians
sr.10.heat.phy <- plot_heatmap(ps.cl.sr.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Phylum",
                           sample.order="site_nice", taxa.order = "Phylum", low = "khaki",
                           high = "red", na.value = "white", title = "Siderastrea radians Top 10 Genera")
sr.10.heat.phy <- sr.10.heat.phy + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
sr.10.heat.phy

#porites
pp.10.heat.phy <- plot_heatmap(ps.cl.pp.10, sample.label="site_nice", method = "NMDS", distance = "bray", taxa.label="Phylum",
                           sample.order="site_nice", taxa.order = "Phylum", low = "khaki",
                           high = "red", na.value = "white", title = "Porites porites Top 10 Genera")
pp.10.heat.phy <- pp.10.heat.phy + theme_classic(base_size = 22)+
  xlab("Site")+
  ggtitle(substitute(paste(italic("Porites spp."))))+
  facet_grid(.~site_nice, scales="free_x")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.title.x = element_blank())
pp.10.heat.phy

#arrange plots together
pp.sr.heat.phy <- ggarrange(sr.10.heat.phy, pp.10.heat.phy, nrow=1, ncol = 2, legend = "none", labels = c("B","C"),font.label = list(size = 40))
heat.phy.all.top.10.phy.ps.less <- ggarrange(ss.10.heat.phy,pp.sr.heat.phy, nrow=2, ncol = 1, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
ggsave(heat.phy.all.top.10.phy.ps.less, file="heatmap.all.top.10.phylum.ps.less.png",width=20,height=10)


#not keeping top 20
ss.20.heat <- plot_heatmap(ps.cl.ss.20, sample.label="site_zone", method = "NMDS", distance = "bray", taxa.label="Genus",
             sample.order="site_zone", taxa.order = "Genus", low = "khaki",
             high = "royalblue1", na.value = "white", title = "Siderastrea Siderea Top 20 Genera (trim)")
ss.20.heat
ggsave(ss.10.heat, file="heat_genus_top10_ps.tr.ss.pdf",h=3,w=8)

ggsave(ss.20.heat, file="heat_genus_top20_ps.tr.ss.pdf",h=3,w=8)

#ssid
ps.cl.ss.vibrio <- subset_taxa(ps.cl.ss, Family=="Vibrionaceae")
plot_heatmap(ps.cl.ss.vibrio)
#ok this one works!
plot_heatmap(ps.cl.ss.vibrio, "NMDS", "bray", "m_y_s_z", "Genus",low = "lightblue",high = "#000033", na.value = "white", sample.order = "m_y_s_z")

heat.sample.ss <- plot_heatmap(ps.cl.ss, method = "NMDS", distance = "bray",
  sample.label = "site_zone", taxa.label = "Genus", low = "#000033",
  high = "lightblue", na.value = "white", trans = log_trans(4),
  max.label = 20)
heat.sample.ss

library(RColorBrewer)
heat.site.ss <- plot_taxa_heatmap(ps.cl.ss,
  subset.top = 20,
  transformation = 'clr',
  VariableA = "site_zone",
  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))

heat.time.ss <- plot_taxa_heatmap(ps.cl.ss,
                                  subset.top = 20,
                                  transformation = 'clr',
                                  VariableA = "m_y",
                                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))

heat.time.site.ss <- plot_taxa_heatmap(ps.cl.ss,
                                  subset.top = 20,
                                  transformation = 'clr',
                                  VariableA = "m_y_s_z",
                                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))


#could only get this to save by using the manual function in R - "export"
#set dimensions as h=6 w=15, and named "heat.ss.site.clr.pdf"

#srad
heat.sample.sr <- plot_taxa_heatmap(ps.cl.sr,
  subset.top = 20,
  transformation = 'clr',
  VariableA = "site_zone",
  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
#saved h=6 w=10, and named "heat.sr.site.clr.pdf"

#ppor
heat.sample.pp <- plot_taxa_heatmap(ps.cl.pp,
  subset.top = 20,
  transformation = 'clr',
  VariableA = "site_zone",
  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
#saved h=6 w=15, and named "heat.pp.site.clr.pdf"

#maya tests of microbiome utilities heatmaps to make them better
heat.sample.ss <- plot_taxa_heatmap(ps.cl.pp,
                                    subset.top = 20,
                                    transformation = 'clr',
                                    VariableA = c("site_zone", "m_y"),
                                    cluster_cols = F,
                                    cluster_rows = F,
                                    show_colnames = F,
                                    heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
#could only get this to save by using the manual function in R - "export"
#set dimensions as h=6 w=15, and named "heat.ss.ps.trim.log10.pdf"

grad_ab <- colorRampPalette(c("orangered","khaki1", "royalblue3"))
heat.cols <- grad_ab(10)
heat.simple <- simple_heatmap(ps.cl.ss,
               group.facet = "m_y",
               group.order = NULL,
               abund.thres = 0.05,
               prev.thres = 0.75,
               level = "Genus",
               scale.color = "log10",
               na.fill = "white",
               color.fill = heat.cols,
               taxa.arrange=TRUE,
               remove.other=TRUE,
               panel.arrange="grid",
               ncol=NULL,
               nrow=NULL)
heat.simple+
  facet_wrap("site_zone"~"m_y")


#srad
heat.sample.sr <- plot_taxa_heatmap(ps.cl.sr,
  subset.top = 20,
  transformation = 'clr',
  VariableA = "site_zone",
  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
#saved h=6 w=10, and named "heat.sr.ps.trim.log10.pdf"

#ppor
heat.sample.pp <- plot_taxa_heatmap(ps.cl.pp,
  subset.top = 20,
  transformation = 'clr',
  VariableA = "site_zone",
  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
#saved h=6 w=15, and named "heat.pp.ps.trim.log10.pdf"







#for emma - to do
#looking at relative abundance of diazotrophs!
#make merged phyloseq that merges by treatment
#do this for apo & sym all together, and then do it for just apo & just sym separately
ps_phylum <- tax_glom(ps.rel, taxrank = "Phylum", NArm = FALSE)
ps_phy_glom_ss <- tax_glom(ps.ss.rel, "Phylum")
ps_fam_glom_ss <- tax_glom(ps.ss.rel, "Family")
phylum.df <- psmelt(ps_phylum)
#head(phylum.df)
#names(myphyla.df)  # to choose factors you want to use to navigate

ss.phy.summary <- psmelt(ps_phy_glom_ss) %>%
  group_by(site_zone, Phylum) %>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 
head(ss.phy.summary)
ss.phy.summary <- data.frame(ss.phy.summary)

ss.phy.summary[ss.phy.summary$Phylum == "Cyanobacteria" & ss.phy.summary$site_zone == "SM_bay" , ]
#& ss.phy.summary$site_zone == "SM_reef" & ss.phy.summary$site_zone == "SW_bay" & ss.phy.summary$site_zone == "SW_reef" , ]
#and then do this for all of the treatments


#MP 25 june 2023
#looking at differences between samples hosting C and D and their microbial community composition
#for PSA talk!
#cool results - they do correlate BUT also could just be site differences bc the sym correlate so strongly with site soooooo who knows

devtools::install_github("david-barnett/microViz")
library("microViz")

#all data from above:
ps.trim.nd <- readRDS("CW_2020_16S_ps.trim.nd.RDS")
ps.cl.ss <- subset_samples(ps.trim.nd,host_species=="siderea")

#Separate by site
ps.cl.ss.SWB <- subset_samples(ps.cl.ss,site_zone=="SW_bay")
ps.cl.ss.SWR <- subset_samples(ps.cl.ss,site_zone=="SW_reef")
ps.cl.ss.SMB <- subset_samples(ps.cl.ss,site_zone=="SM_bay")
ps.cl.ss.SMR <- subset_samples(ps.cl.ss,site_zone=="SM_reef")

#MAJORITY C
#"10SSSMBMarch2020"|"86SSSMBMarch2020"|"70SSSMBNovember2020"|"39SSSMBNovember2021"
#"77SSSWBMarch2020"

#MAJORITY D
#"3SSSMRMarch2020"|"62SSSMRMarch2020"|"48SSSMRNovember2020"|"15SSSMRNovember2021"|"16SSSMRNovember2021"|"19SSSMRNovember2020"
##"1SSDBMarch2020"|"49SSDBNovember2020"|"96SSDBNovember2020"|"85SSDBNovember2021"

#i'm so so sure there's a better way to do this and I'm annoyed this is how I'm doing it but here goes
#SWB 1 C sample, rest D
ps.cl.ss.swb.cd <- ps.cl.ss.SWB %>%
  ps_mutate(major_clade = case_when(
    full_sample_id == "77SSSWBMarch2020" ~ "Majority C",
    full_sample_id != "77SSSWBMarch2020" ~ "Majority D"))

#SWR 4 D samples, rest C
ps.cl.ss.swr.cd <- ps.cl.ss.SWR %>%
  ps_mutate(major_clade = case_when(
    full_sample_id == "1SSDBMarch2020"|full_sample_id =="49SSDBNovember2020"|full_sample_id =="96SSDBNovember2020"|full_sample_id =="85SSDBNovember2021" ~ "Majority D",
    full_sample_id != "1SSDBMarch2020"|full_sample_id != "49SSDBNovember2020"|full_sample_id != "96SSDBNovember2020"|full_sample_id != "85SSDBNovember2021" ~ "Majority C"))

#SMB 4 C samples, rest D
ps.cl.ss.smb.cd <- ps.cl.ss.SMB %>%
  ps_mutate(major_clade = case_when(
    full_sample_id == "100SSSMBMarch2020"|full_sample_id =="86SSSMBMarch2020"|full_sample_id =="70SSSMBNovember2020"|full_sample_id =="39SSSMBNovember2021" ~ "Majority C",
    full_sample_id != "100SSSMBMarch2020"|full_sample_id !="86SSSMBMarch2020"|full_sample_id !="70SSSMBNovember2020"|full_sample_id !="39SSSMBNovember2021" ~ "Majority D"))
#10 was labeled as 100 for some reason? check labeling, i think it should be 100

#SMR - 3 D samples, rest C
ps.cl.ss.smr.cd <- ps.cl.ss.SMR %>%
  ps_mutate(major_clade = case_when(
    full_sample_id == "3SSSMRMarch2020"|full_sample_id =="62SSSMRMarch2020"|full_sample_id =="15SSSMRNovember2021" ~ "Majority D",
    full_sample_id != "3SSSMRMarch2020"|full_sample_id !="62SSSMRMarch2020"|full_sample_id !="15SSSMRNovember2021"  ~ "Majority C"))
#sample IDs that are not in the dataset (bc trimmed earlier, remember weird bummer time of sequencing low coverage #tragic #disaster):
#"48SSSMRNovember2020"|"16SSSMRNovember2021"|"19SSSMRNovember2020"

#TOTAL 5 Bay C samples, 7 Reef D samples

#merge phyloseq objects back together
ps.cl.ss.cd <- merge_phyloseq(ps.cl.ss.swb.cd,ps.cl.ss.swr.cd,ps.cl.ss.smb.cd,ps.cl.ss.smr.cd)

#might need to just do this later after mutating
ps.cl.ss.cd.rel <- transform_sample_counts(ps.cl.ss.cd, function(x) x / sum(x))
seq.cl.ss.cd.rel <- data.frame(otu_table(ps.cl.ss.cd.rel))
#also non-relative seq table
seq.cl.ss.cd <- data.frame(otu_table(ps.cl.ss.cd))

# Bar plots - ALL
ps_phy_glom_ss_cd <- tax_glom(ps.cl.ss.cd.rel, "Phylum")
ps_gen_glom_ss_cd <- tax_glom(ps.cl.ss.cd.rel, "Genus")

gg.bar.ss.cd <- plot_bar(ps_phy_glom_ss_cd,"sample_full",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(major_clade~site_zone, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Relative Abundance")
gg.bar.ss.cd = gg.bar.ss.cd + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.ss.cd
ggsave(gg.bar.ss.cd,file="bac.phy.barplot.ss.cd.pdf",h=10,w=20)
#full plot - to put in supplementary, with all data across sites and time
gg.bar.ss <- plot_bar(ps_phy_glom_ss,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.ss
#each number corresponds to specific sample here
ggsave(gg.bar.ss,file="bac.phy.barplot.ss.site.season.sampleid.pdf",h=15,w=30)   

#bars grouped by symbiont type
ps.ss.cd <- merge_samples(ps_phy_glom_ss_cd, "major_clade")
ps.ss.cd.rel <- transform_sample_counts(ps.ss.cd, function(x) x / sum(x))
bar.rel.cd.ss <- plot_bar(ps.ss.cd.rel, fill="Phylum")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  xlab("Symbiodiniaceae Genus")+
  ylab("Relative Abundance")+
  #scale_fill_manual(values=c("#FDBF6F","#9ECAE1","#E78AC3","#B2DF8A","#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#CAB2D6"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar.rel.cd.ss
ggsave(bar.rel.cd.ss, file="bar.rel.cd.ss.pdf",w=10,h=8)

#ord plots for SSID to look at beta diversity

ss.ord.all.cd <- plot_ordination(ps.cl.ss.cd,ordinate(ps.cl.ss.cd,"PCoA", "bray"),color="major_clade")+
  stat_ellipse()+
  scale_color_brewer(palette="Set1")
ss.ord.all.cd
ggsave(ss.ord.all.cd, file="ss.ord.all.cd.ps.trim.pdf",w=10,h=8)
#with stats - can only look at differences between C and D, not combined C and D and site bc too few
#but here there also seems to be something special about the reef sids hosting D!!! but still it is only 4 samples

#ok now to do some stats on this shit
#ssid
dist.cd.trim <- vegdist(seq.cl.ss.cd)
samdf.cd.trim <- data.frame(sample_data(ps.cl.ss.cd))
row.names(samdf.cd.trim)==row.names(seq.cl.ss.cd)
bet.cd.ss <- betadisper(dist.cd.trim,samdf.cd.trim$major_clade)
anova(bet.cd.ss) #ns for clade
plot(bet.cd.ss) #not sig, very lots of overlap, nice
permutest(bet.cd.ss, pairwise = TRUE, permutations = 999) #ns
adonis2(dist.cd.trim ~ major_clade, data=samdf.cd.trim, permutations=999) 
#clade *!!!! p=0.011
adonis2(dist.cd.trim ~ major_clade, data=samdf.cd.trim, permutations=999) #p=0.011
adonis2(dist.cd.trim ~ major_clade*site_zone, data=samdf.cd.trim, permutations=999) #yes!! both sig clade p=0.009, site_zone p=0.001 but intxn not sig
adonis2(dist.cd.trim ~ major_clade, strata=samdf.cd.trim$site_zone, data=samdf.cd.trim, permutations=999) #ns
pairwise.adonis2(dist.cd.trim ~ site_zone*major_clade, data=samdf.cd.trim, permutations = 999)
#major clade is significant, unless you put site differences first, so site is driving differences in bacterial community more

#PCA of C vs D
#ssid
ps.cd.ss.pca <- ps.cl.ss.cd %>% # replace this with your phyloseq object
  tax_transform("clr", rank = "Phylum") %>% # change here if you want to plot by a different taxonomic rank
  ord_calc()
ss.cd.pca <- ps.cd.ss.pca  %>% # replace es_hsp_pca with the object you saved above with your transformed phyloseq object
  ord_plot(color = "major_clade", # sample_data variable to color by
           plot_taxa = 1:64, # number of taxa to label on plot, I change
           size=3.5, # size of data points
           tax_vec_length = 1.5, # length of arrows corresponding to taxa vectors
           tax_lab_length = 1.6, # length of taxa label, keep just slightly longer than tax_vec_length
           tax_lab_style = tax_lab_style(type = "text", check_overlap=TRUE, # check_overlap is key, keeps your labels from overlapping
                                         max_angle=90, size=2.8, fontface = "bold.italic")
  ) +
  #theme_biome_utils() +
  scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))+
  ggtitle(label = "Siderastrea siderea PCA") + # obv change to whatever you'd like for the title of your plot
  #theme(plot.title = element_markdown()) +
  stat_ellipse(aes(color=major_clade)) + # sample.type = the sample_data variable I chose to calculate ellipses by
  coord_fixed(ratio = 1, clip = "off", xlim = c(-5, 5), ylim = c(-5, 5)) # I usually just plot a few times and manually change the xlim/ylim coordinates so that everything is nicely contained and my ellipses are within the axes
ss.cd.pca
ggsave(ss.cd.pca, file="ss.phylum.PCA.cd.ps.trim.pdf", height=8, width=8)

#differential abundance using ALDEx2
#paper that looks at all the methods and says ancom2 and aldex2 are the best and aldex2 works better for me so I'm using it YEAH
#BiocManager::install("ALDEx2")
library(ALDEx2)
BiocManager::install("microbiomeMarker")
library(microbiomeMarker)
#ISSUE: ALDEx2 doesn't work anymore giving this set of errors:
#https://github.com/yiluheihei/microbiomeMarker/issues/114
#think I need to use old version of ALDEx2 v1.28.0
#got it at this link https://mghp.osn.xsede.org/bir190004-bucket01/index.html#archive.bioconductor.org/packages/3.15/bioc/src/contrib/Archive/ALDEx2/
remove.packages("ALDEx2")
remove.packages("microbiomeMarker")
install.packages("~/Downloads/archive.bioconductor.org_packages_3.15_bioc_src_contrib_Archive_ALDEx2_ALDEx2_1.28.0.tar.gz", repos = NULL, type="source")
library("ALDEx2") #shows as 1.28.0 in my system library in the packages tab
BiocManager::install("microbiomeMarker")
library("microbiomeMarker")

#aldex2 for species comparisons
aldex2_spp <- run_aldex(ps.less.nd, group="host_species", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(), method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_spp <-marker_table(aldex2_spp)

#run aldex2 for ssid
aldex2_ss <- run_aldex(ps.cl.ss, group="site_zone", taxa_rank = "Genus",
  transform = "identity", norm = "none", norm_para = list(), method = "glm_anova",
  p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_ss <-marker_table(aldex2_ss)
#aldex2 ssid timepoint
aldex2_ss_time <- run_aldex(ps.cl.ss, group="m_y", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_ss_time <-marker_table(aldex2_ss_time)

#aldex2 ssid m_rb
aldex2_ss_m_rb <- run_aldex(ps.cl.ss, group="m_rb", taxa_rank = "Genus",
                            transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                            p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_ss_m_rb <-marker_table(aldex2_ss_m_rb)

#aldex2 ssid m_y_s_z
aldex2_ss_all <- run_aldex(ps.cl.ss, group="m_y_s_z", taxa_rank = "Genus",
                            transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                            p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_ss_all <-marker_table(aldex2_ss_all)

#run aldex2 for ssid reef vs bay
aldex2_ss_bay_reef <- run_aldex(ps.cl.ss, group="reef_bay", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "wilcox.test",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_ss_bay_reef <-marker_table(aldex2_ss_bay_reef)

#aldex2 location for srad
aldex2_sr <- run_aldex(ps.cl.sr, group="site_zone", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_sr <-marker_table(aldex2_sr)

#srad aldex timepoint
aldex2_sr_time <- run_aldex(ps.cl.sr, group="m_y", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_sr_time <-marker_table(aldex2_sr_time)

#srad aldex location:timepoint
aldex2_sr_time_site <- run_aldex(ps.cl.sr, group="m_y_s_z", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
mm_sr_time_site <-marker_table(aldex2_sr_time_site)

#aldex2 for ppor
aldex2_pp <- run_aldex(ps.cl.pp, group="site_zone", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
aldex2_pp_time <- run_aldex(ps.cl.pp, group="m_y", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)
aldex2_pp_time_zone <- run_aldex(ps.cl.pp, group="m_y_s_z", taxa_rank = "Genus",
                       transform = "identity", norm = "none", norm_para = list(),method = "glm_anova",
                       p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "all", paired = FALSE)

#no markers found for porites both at site and timepoint

#interesting! maybe por is more driven by symbionts
#talking with carly scott - how much do corals have "control" over sym vs microbes
#transient vs not - here are porites brooders so 'vertical' sym transmission
#maybe more control/more "STABILITY" whatever that means over microbes too???
#what is the impetus to actually switch sym vs microbe? like reaction energy req kinda thing

saveRDS(aldex2_ss_time, "CW_2020_16S_aldex2_ss_time.RDS")
saveRDS(aldex2_ss, "CW_2020_16S_aldex2_ss.RDS")
saveRDS(aldex2_sr, "CW_2020_16S_aldex2_sr.RDS")
saveRDS(aldex2_ss_bay_reef, "CW_2020_16S_aldex2_ss_bay_reef.RDS")

aldex2_ss <- readRDS("CW_2020_16S_aldex2_ss.RDS")
aldex2_sr <- readRDS("CW_2020_16S_aldex2_sr.RDS")
aldex2_ss_bay_reef <- readRDS("CW_2020_16S_aldex2_ss_bay_reef.RDS")
aldex2_ss_time <- readRDS("CW_2020_16S_aldex2_ss_time.RDS")

#abundance plots
abundance_aldex2_ss <- plot_abundance(aldex2_ss, markers = NULL, group="site_zone") + 
  scale_fill_manual("Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"), labels=c('Santa Martha Bay', 'Santa Martha Reef', 'Spaanse Water Bay','Spannse Water Reef'))+
  theme_classic(base_size = 20)+
  xlab("Abundance")+
  scale_x_continuous(breaks = seq(0, 10000, by = 2000))+
  theme(legend.position="bottom",legend.key.size = unit(1, 'cm'), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color="black"))
abundance_aldex2_ss

ggsave(abundance_aldex2_ss, file = "aldex2_abundance_ss_sites_ps.less.pdf",h=8,w=15)

abund_aldex2_ss_bay_reef <- plot_abundance(aldex2_ss_bay_reef, markers = NULL, group="reef_bay") +
  scale_fill_manual("Habitat", values = c("coral","skyblue3"), labels=c('Bay', 'Reef'))+
  geom_boxplot(inherit.aes = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 20)+
  xlab("Abundance")+
  scale_x_continuous(breaks = seq(0, 10000, by = 2000))+
  theme(legend.position="bottom", legend.key.size = unit(1, 'cm'), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color="black"))
abund_aldex2_ss_bay_reef

ggsave(abund_aldex2_ss_bay_reef, file = "aldex2_abundance_ss_bay_reef_ps.less.pdf",h=8,w=10)

abundance_aldex2_ss_time <- plot_abundance(aldex2_ss_time, markers = NULL, group="m_y") + 
  scale_fill_manual("Timepoint", values = c("goldenrod","cornflowerblue", "skyblue"),labels = c("March 2020", "November 2020", "November 2021"))+
  theme_classic(base_size = 20)+
  xlab("Abundance")+
  scale_x_continuous(breaks = seq(0, 20000, by = 5000))+
  theme(legend.position="bottom",legend.key.size = unit(1, 'cm'), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color="black"))
abundance_aldex2_ss_time

ggsave(abundance_aldex2_ss_time, file = "aldex2_abundance_ss_time_ps.less.pdf",h=8,w=15)

abundance_aldex2_sr <- plot_abundance(aldex2_sr, markers = NULL, group="site_zone")+
  scale_fill_manual("Site", values = c("coral3","lightsalmon1"), labels=c('Santa Martha Bay','Spaanse Water Bay'))+
  theme_classic(base_size = 20)+
  xlab("Abundance")+
  scale_x_continuous(breaks = seq(0, 12000, by = 2000))+
  theme(legend.position="bottom",legend.key.size = unit(1, 'cm'), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color="black"))
abundance_aldex2_sr <- annotate_figure(abundance_aldex2_sr, top = text_grob("Siderastrea radians", color = "black", hjust = 1, face = "italic", size = 40))
ggsave(abundance_aldex2_sr, file = "aldex2_abundance_sr_sites_ps.less.pdf",h=8,w=10)

heat_aldex2_ss <- microbiomeMarker::plot_heatmap(aldex2_ss, transform = "identity", group="site_zone")
heat_aldex2_ss
#had to save this via export to pdf bc it wouldn't work for some reason!!
ggsave(heat_aldex2_ss, file = "heatmap_diff_abundance_site_ss.ps.less.pdf",h=4,w=10)

heat_aldex2_ss_bay_reef <- microbiomeMarker::plot_heatmap(aldex2_ss_bay_reef, transform = "identity", group="reef_bay")
heat_aldex2_ss_bay_reef

heat_aldex2_sr <- microbiomeMarker::plot_heatmap(aldex2_sr,transform="identity",group="site_zone")
heat_aldex2_sr
#also won't save so had to just export
ggsave(heat_aldex2_sr, file="heatmap_diff_abundance_site_sr.ps.less.png",h=4,w=8)

#running ancombc for ssid
run_ancombc(ps.cl.ss, group="site_zone", confounders = character(0), contrast = NULL,
  taxa_rank = "Genus", transform = "identity", norm = "none", p_adjust = "BH",
  prv_cut = 0.1, tol = 1e-05, max_iter = 100, conserve = FALSE, pvalue_cutoff = 0.05)
#running this with no normalization, but could also do CLR normalization... idk what is better
#found no markers so use aldex ig

#emma and maya making pcoa with taxa labels for diazotroph project
#packages
library(dplyr)
install.packages("ggside")
library(ggside)
BiocManager::install("DESeq2")
install.packages("DESeq2")
library(DESeq2)
BiocManager::install("DECIPHER")
library(DECIPHER)
devtools::install_github("RIVM-IIV-Microbiome/biomeViz")
library(biomeViz)
devtools::install_github("RIVM-IIV-Microbiome/biomeUtils", force = TRUE)
library(biomeUtils)
ps.maya <- ps.rare %>%
  microbiome::transform("compositional") %>%
  mutateTaxaTable(FeatureID = taxa_names(ps.rare))
taxa.pcoa <- plotPCoA(x =ps.maya,
                      group_var = "doc.conc",
                      ord_method = "PCoA",
                      dist_method = "bray",
                      seed = 1253,
                      cor_method = "spearman",
                      verbose = TRUE,
                      padj_cutoff = 0.05,
                      padj_method = "fdr",
                      arrows = TRUE,
                      label_col = "grey30",
                      plot_centroids = TRUE,
                      add_side_box = FALSE,
                      axis_plot = c(1:2),
                      point_shape = 21,  # point_shape
                      point_alpha = 1) +
  scale_color_manual(name = "DOC Concentration", values = c("coral3", "lightcyan3","skyblue1","skyblue3","deepskyblue4"),labels = c("Before", "0 ppm", "5 ppm", "10 ppm", "20 ppm")) +
  theme_classic(base_size = 22) +
  scale_fill_manual(name = "DOC Concentration", values = c("coral3", "lightcyan3","skyblue1","skyblue3","deepskyblue4"),labels = c("Before", "0 ppm", "5 ppm", "10 ppm", "20 ppm"))
taxa.pcoa

##not using anymore
### Core
##go through each dataset

sums <- as.data.frame(taxa_sums(ps.all.rel))
colnames(sums)[1] <- 'sum'
low_sums <- subset(sums, sum <0.0001)

low_names <- row.names(low_sums)

no.rare.seq <- seq.all.rel[,!colnames(seq.all.rel) %in% low_names ]
#remake phyloseq object - already relative
ps.no.rare <- phyloseq(otu_table(no.rare.seq, taxa_are_rows=FALSE), 
                   sample_data(samdf.less), 
                   tax_table(taxa2))
# ps.core.all <- core(ps.trim.nd, detection = 0, prevalence = .7)
# core.all.tax <- data.frame(ps.core.all@tax_table) #7 taxa ps.cleanest
# dim(core.all.tax) #9 ps.trim
# ps.core.ss <- core(ps.cl.ss, detection = 0, prevalence = .7)
# core.ss.tax <- data.frame(ps.core.ss@tax_table) #11 taxa 79 samples ps.cleanest
# dim(core.ss.tax) #ps.trim 13 taxa
# ps.core.sr <- core(ps.cl.sr, detection = 0, prevalence = .7)
# core.sr.tax <- data.frame(ps.core.sr@tax_table) #11 taxa 34 samples
# dim(core.sr.tax) #ps.trim 11 taxa
# ps.core.pp <- core(ps.cl.pp, detection = 0, prevalence = .7)
# core.pp.tax <- data.frame(ps.core.pp@tax_table) #5 taxa 15 samples
# dim(core.pp.tax) #4 taxa ps.trim
# ps.core.pa <- core(ps.cl.pa, detection = 0, prevalence = .7)
# core.pa.tax <- data.frame(ps.core.pa@tax_table) #1 taxa 6 samples
# dim(core.pa.tax) #ps.trim 1 taxa
# 
# # calculating core abundances 
# #all together
# #do this with each dataset to see differences
# core.all.ids <- row.names(core.all.tax)
# core.all.tax$id <-row.names(core.all.tax)
# seq.all.rel <- readRDS("seq.all.rel.ps.trim.RDS")
# all.tax <- data.frame(tax_table(ps.all.rel))
# all.tax$id <- rownames(all.tax)
# seq.all.core <- seq.all.rel %>% select(all_of(core.all.ids))
# #core
# core.all.rel <- data.frame(colMeans(seq.all.core))
# core.all.rel <- core.all.rel %>% 
#   dplyr::rename(relative_abundance = colMeans.seq.all.core.) %>% 
#   mutate(id = rownames(core.all.rel)) %>%
#   mutate(species = "all")
# core.all.rel <- merge(core.all.tax,core.all.rel,by="id")
# core.all.rel.ordered <- core.all.rel %>%
#   arrange(desc(relative_abundance))
# #total
# total.all.rel <- data.frame(colMeans(seq.all.rel))
# total.all.rel <- total.all.rel %>% 
#   dplyr::rename(relative_abundance = colMeans.seq.all.rel.) %>% 
#   mutate(id = rownames(total.all.rel)) %>%
#   mutate(species = "all") %>%
#   filter(!relative_abundance==0)
# total.all.rel <- merge(all.tax,total.all.rel,by="id")
# total.all.rel.ordered <- total.all.rel %>%
#   arrange(desc(relative_abundance))
# #now save as tables both core and all
# write.csv(total.all.rel.ordered,"total.all.rel.ordered.ps.trim.csv")
# 
# #ssid
# core.ss.ids <- row.names(core.ss.tax)
# core.ss.tax$id <-row.names(core.ss.tax)
# ps.ss.rel <- transform_sample_counts(ps.cl.ss, function(x) x / sum(x))
# seq.ss.rel <- readRDS("seq.ss.rel.ps.trim.RDS")
# ss.tax <- data.frame(tax_table(ps.ss.rel))
# ss.tax$id <- rownames(ss.tax)
# seq.ss.core <- seq.ss.rel %>% select(all_of(core.ss.ids))
# #core
# core.ss.rel <- data.frame(colMeans(seq.ss.core))
# core.ss.rel <- core.ss.rel %>% 
#   rename(relative_abundance = colMeans.seq.ss.core.) %>% 
#   mutate(id = rownames(core.ss.rel)) %>%
#   mutate(species = "Siderastrea siderea")
# core.ss.rel <- merge(core.ss.tax,core.ss.rel,by="id")
# core.ss.rel.ordered <- core.ss.rel %>%
#   arrange(desc(relative_abundance))
# #total
# total.ss.rel <- data.frame(colMeans(seq.ss.rel))
# total.ss.rel <- total.ss.rel %>% 
#   rename(relative_abundance = colMeans.seq.ss.rel.) %>% 
#   mutate(id = rownames(total.ss.rel)) %>%
#   mutate(species = "Siderastrea siderea")%>%
#   filter(!relative_abundance==0)
# total.ss.rel <- merge(ss.tax,total.ss.rel,by="id")
# total.ss.rel.ordered <- total.ss.rel %>%
#   arrange(desc(relative_abundance))
# #now save as tables both core and ss
# write.csv(total.ss.rel.ordered,"total.ss.rel.ordered.ps.trim.csv")
# 
# #srad
# core.sr.ids <- row.names(core.sr.tax)
# core.sr.tax$id <-row.names(core.sr.tax)
# ps.sr.rel <- transform_sample_counts(ps.cl.sr, function(x) x / sum(x))
# seq.sr.rel <- readRDS("seq.sr.rel.ps.trim.RDS")
# sr.tax <- data.frame(tax_table(ps.sr.rel))
# sr.tax$id <- rownames(sr.tax)
# seq.sr.core <- seq.sr.rel %>% select(all_of(core.sr.ids))
# #core
# core.sr.rel <- data.frame(colMeans(seq.sr.core))
# core.sr.rel <- core.sr.rel %>% 
#   rename(relative_abundance = colMeans.seq.sr.core.) %>% 
#   mutate(id = rownames(core.sr.rel)) %>%
#   mutate(species = "Siderastrea radians")
# core.sr.rel <- merge(core.sr.tax,core.sr.rel,by="id")
# core.sr.rel.ordered <- core.sr.rel %>%
#   arrange(desc(relative_abundance))
# #total
# total.sr.rel <- data.frame(colMeans(seq.sr.rel))
# total.sr.rel <- total.sr.rel %>% 
#   rename(relative_abundance = colMeans.seq.sr.rel.) %>% 
#   mutate(id = rownames(total.sr.rel)) %>%
#   mutate(species = "Siderastrea radians")%>%
#   filter(!relative_abundance==0)
# total.sr.rel <- merge(sr.tax,total.sr.rel,by="id")
# total.sr.rel.ordered <- total.sr.rel %>%
#   arrange(desc(relative_abundance))
# #now save as tables both core and sr
# write.csv(total.sr.rel.ordered,"total.sr.rel.ordered.ps.trim.csv")
# 
# #ppor
# core.pp.ids <- row.names(core.pp.tax)
# core.pp.tax$id <-row.names(core.pp.tax)
# ps.pp.rel <- transform_sample_counts(ps.cl.pp, function(x) x / sum(x))
# seq.pp.rel <- readRDS("seq.pp.rel.ps.trim.RDS")
# pp.tax <- data.frame(tax_table(ps.pp.rel))
# pp.tax$id <- rownames(pp.tax)
# seq.pp.core <- seq.pp.rel %>% select(all_of(core.pp.ids))
# #core
# core.pp.rel <- data.frame(colMeans(seq.pp.core))
# core.pp.rel <- core.pp.rel %>% 
#   rename(relative_abundance = colMeans.seq.pp.core.) %>% 
#   mutate(id = rownames(core.pp.rel)) %>%
#   mutate(species = "Porites porites")
# core.pp.rel <- merge(core.pp.tax,core.pp.rel,by="id")
# core.pp.rel.ordered <- core.pp.rel %>%
#   arrange(desc(relative_abundance))
# #total
# total.pp.rel <- data.frame(colMeans(seq.pp.rel))
# total.pp.rel <- total.pp.rel %>% 
#   rename(relative_abundance = colMeans.seq.pp.rel.) %>% 
#   mutate(id = rownames(total.pp.rel)) %>%
#   mutate(species = "Porites porites") %>%
#   filter(!relative_abundance==0)
# total.pp.rel <- merge(pp.tax,total.pp.rel,by="id")
# total.pp.rel.ordered <- total.pp.rel %>%
#   arrange(desc(relative_abundance))
# #now save as tables both core and pp
# write.csv(total.pp.rel.ordered,"total.pp.rel.ordered.ps.trim.csv")
# 
# #past
# core.pa.ids <- row.names(core.pa.tax)
# core.pa.tax$id <-row.names(core.pa.tax)
# seq.pa.rel <- readRDS("seq.pa.rel.ps.trim.RDS")
# pa.tax <- data.frame(tax_table(ps.pa.rel))
# pa.tax$id <- rownames(pa.tax)
# seq.pa.core <- seq.pa.rel %>% select(all_of(core.pa.ids))
# #core
# core.pa.rel <- data.frame(colMeans(seq.pa.core))
# core.pa.rel <- core.pa.rel %>% 
#   rename(relative_abundance = colMeans.seq.pa.core.) %>% 
#   mutate(id = rownames(core.pa.rel)) %>%
#   mutate(species = "Porites astreoides") 
# core.pa.rel <- merge(core.pa.tax,core.pa.rel,by="id")
# core.pa.rel.ordered <- core.pa.rel %>%
#   arrange(desc(relative_abundance))
# #total
# total.pa.rel <- data.frame(colMeans(seq.pa.rel))
# total.pa.rel <- total.pa.rel %>% 
#   rename(relative_abundance = colMeans.seq.pa.rel.) %>% 
#   mutate(id = rownames(total.pa.rel)) %>%
#   mutate(species = "Porites astreoides")%>%
#   filter(!relative_abundance==0)
# total.pa.rel <- merge(pa.tax,total.pa.rel,by="id")
# total.pa.rel.ordered <- total.pa.rel %>%
#   arrange(desc(relative_abundance))
# #now save as tables both core and pa
# write.csv(total.pa.rel.ordered,"total.pa.rel.ordered.ps.trim.csv")
# 
# #now put them all together!
# core.rel.ordered.all.types <- rbind(core.all.rel.ordered,core.ss.rel.ordered,core.sr.rel.ordered,core.pp.rel.ordered,core.pa.rel.ordered)
# write.csv(core.rel.ordered.all.types,"core.rel.ordered.all.types.ps.trim.csv")
# 
# #time for some graphs!!!
# #all and by species
# #colors:
# #display.brewer.all(colorblindFriendly = TRUE)
# #brewer.pal(12, "Paired")
# "#A6CEE3" "#1F78B4" "#B2DF8A"  "#FB9A99" "#E31A1C"  "#FF7F00"
# "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
# 
# "#FDBF6F"#Aestuariibacter
# "#9ECAE1"#Alteromonas 
# "#33A02C"#Class_Bact
# "#E78AC3"#Photobacterium
# "#B2DF8A"#Porphyro
# "#FFD92F"#Pseudoalteromonas 
# "#1F78B4"#Psychrobacter
# "#8DA0CB"#Ruegeria 
# "#FC8D62"#Streptococcus 
# "#FFFF99"#Synnechococcus
# "#66C2A5"#Thalassotalea
# "#CAB2D6"#Vibrio
# 
# c("#FDBF6F","#9ECAE1","#E78AC3","#B2DF8A","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6","#66C2A5")
# 
# #all species together
# ps_glom <- tax_glom(ps.core.all, "Genus")
# ps.z <- merge_samples(ps_glom, "m_y_s_z")
# ps.rel.z <- transform_sample_counts(ps.z, function(x) x / sum(x))
# bar.core.rel.all <- plot_bar(ps.rel.z, fill="Genus")+
#   geom_bar(stat="identity")+
#   theme_cowplot()+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   #scale_fill_manual(values=c("#9ECAE1","#E78AC3","#FFD92F","#8DA0CB","#FC8D62","#CAB2D6"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# bar.core.rel.all
# ggsave(bar.core.rel.all, file="core.bar.ps.trim.all.pdf",w=8,h=8)
# 
# #ssid
# ps_glom.ss <- tax_glom(ps.core.ss, "Genus")
# ps.ss.z <- merge_samples(ps_glom.ss, "site_zone")
# ps.ss.rel.z <- transform_sample_counts(ps.ss.z, function(x) x / sum(x))
# bar.core.rel.ss <- plot_bar(ps.ss.rel.z, fill="Genus")+
#   geom_bar(stat="identity")+
#   theme_cowplot()+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   scale_fill_manual(values=c("#FDBF6F","#9ECAE1","#E78AC3","#B2DF8A","#FFD92F","#8DA0CB","#FC8D62","#66C2A5","#CAB2D6"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# bar.core.rel.ss
# ggsave(bar.core.rel.ss, file="core.bar.ps.trim.ss.site.pdf",w=10,h=8)
# #srad
# ps_glom.sr <- tax_glom(ps.core.sr, "Species")
# ps.sr.z <- merge_samples(ps_glom.sr, "m_y_s_z")
# ps.sr.rel.z <- transform_sample_counts(ps.sr.z, function(x) x / sum(x))
# bar.core.rel.sr <- plot_bar(ps.sr.rel.z, fill="Genus")+
#   geom_bar(stat="identity")+
#   theme_cowplot()+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   #scale_fill_manual(values=c("#FDBF6F","#9ECAE1","#33A02C","#E78AC3","#8DA0CB","#FC8D62","#FFFF99","#66C2A5","#CAB2D6"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# bar.core.rel.sr
# ggsave(bar.core.rel.sr, file="core.bar.ps.trim.sr.pdf",w=6,h=8)
# #ppor
# ps_glom.pp <- tax_glom(ps.core.pp, "Genus")
# ps.pp.z <- merge_samples(ps_glom.pp, "m_y_s_z")
# ps.pp.rel.z <- transform_sample_counts(ps.pp.z, function(x) x / sum(x))
# bar.core.rel.pp <- plot_bar(ps.pp.rel.z, fill="Genus")+
#   geom_bar(stat="identity")+
#   theme_cowplot()+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   #scale_fill_manual(values=c("#FFD92F","#1F78B4","#FC8D62","#CAB2D6"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# bar.core.rel.pp
# ggsave(bar.core.rel.pp, file="core.bar.ps.trim.pp.pdf",w=6,h=8)
# #past
# ps_glom.pa <- tax_glom(ps.core.pa, "Genus")
# ps.pa.z <- merge_samples(ps_glom.pa, "site_zone")
# ps.pa.rel.z <- transform_sample_counts(ps.pa.z, function(x) x / sum(x))
# bar.core.rel.pa <- plot_bar(ps.pa.rel.z, fill="Genus")+
#   geom_bar(stat="identity")+
#   theme_cowplot()+
#   #scale_fill_manual(values=c("#6A3D9A"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# bar.core.rel.pa
# ggsave(bar.core.rel.pa, file="core.bar.ps.trim.pa.pdf",w=4,h=8)
# 
# #now ordination plots - they are all looking so crazy and yet somehow are not significant
# #ssid ord core
# ss.ord.core <- plot_ordination(ps.core.ss,ordinate(ps.core.ss,"PCoA", "bray"),color="site_zone")+
#   stat_ellipse()+
#   theme_cowplot()+
#   scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB","#FC8D62","#66C2A5"))
# ss.ord.core
# #looks like a hot mess!!
# #grouping by time point and site also doesn't show any patterns/clustering
# ggsave(ss.ord.core, file="ss.ord.core.ps.trim.pdf",w=8,h=8)
# #srad ord core
# sr.ord.core <- plot_ordination(ps.core.sr,ordinate(ps.core.sr,"PCoA", "bray"),color="site_zone")+
#   stat_ellipse()+
#   theme_cowplot()+
#   scale_color_manual(name="Site",values=c("#E78AC3","#FC8D62"))
# sr.ord.core
# #no patterns, grouping by timepoint shows that march 2020 SMB clusters tightly
# ggsave(sr.ord.core, file="sr.ord.core.ps.trim.pdf",w=8,h=8)
# #ppor ord core
# pp.ord.core <- plot_ordination(ps.core.pp,ordinate(ps.core.pp,"PCoA", "bray"),color="site_zone")+
#   stat_ellipse()+
#   theme_cowplot()+
#   scale_color_manual(name="Site",values=c("#E78AC3","#8DA0CB"))
# pp.ord.core
# #no patterns, grouping by timepoint shows that Nov 2021 SMR clusters slightly more
# ggsave(pp.ord.core, file="pp.ord.core.ps.trim.pdf",w=8,h=8)
# #past ord core
# pa.ord.core <- plot_ordination(ps.core.pa,ordinate(ps.core.pa,"PCoA", "bray"),color="site_zone")+
#   stat_ellipse()+
#   theme_cowplot()+
#   scale_color_manual(name="Site",values=c("#FC8D62","#66C2A5"))
# pa.ord.core
# #yeah clusters for SWB but only 6 points not really relevant or interesting
# ggsave(pa.ord.core, file="pa.ord.core.ps.trim.pdf",w=8,h=8)
# 
# 
# ### Core stats
# #all
# dist.core.all <- vegdist(seq.all.core)
# samdf.core.all <- data.frame(sample_data(ps.all.rel))
# row.names(samdf.core.all)==row.names(seq.all.rel)
# set.seed(123)
# bet.core.all <- betadisper(dist.core.all,samdf.core.all$site_zone)
# anova(bet.core.all) #ps.cleanest ns, ps.trim ns
# plot(bet.core.all)
# permutest(bet.core.all, pairwise = TRUE, permutations = 999)
# #ps.cleanest ns, ps.trim ns
# adonis2(dist.core.all ~ site_zone, data=samdf.core.all, permutations=999)
# #ps.cleanest ***, ps.trim ***
# adonis2(dist.core.all ~ site_zone/m_y, data=samdf.core.all, permutations=999)
# #ps.cleanest s_z **, s_z:m_y ***, ps. trim both ***
# pairwise.adonis2(dist.core.all ~ site_zone, data=samdf.core.all)
# #PCOA all very overlapped
# #all sig except SMB vs SMR just site (ps.cleanest, ps.trim)
# 
# #Ssid
# dist.core.ss <- vegdist(seq.ss.core)
# samdf.core.ss <- data.frame(sample_data(ps.ss.rel))
# row.names(samdf.core.ss)==row.names(seq.ss.rel)
# set.seed(123)
# bet.core.ss <- betadisper(dist.core.ss,samdf.core.ss$site_zone)
# anova(bet.core.ss) #ps.cleanest, ps.trim ns 
# plot(bet.core.ss)
# permutest(bet.core.ss, pairwise = TRUE, permutations = 999)
# #ps.cleanest,ps.trim ns
# adonis2(dist.core.ss ~ site_zone, data=samdf.core.ss, permutations=999)
# #ps.cleanest, ps.trim **
# adonis2(dist.core.ss ~ site_zone/m_y, data=samdf.core.ss, permutations=999)
# #ps.cleanest, ps.trim s_z **, s_z:m_y ***
# pairwise.adonis2(dist.core.ss ~ site_zone/m_y, data=samdf.core.ss)
# #all sig except SMB vs SMR just site ps.cleanest
# #all sig except SW_bay_vs_SM_bay ps.trim
# #PCOA all very overlapped
# #betadisper: homogenous dispersion (ns)
# #adonis: significant, and can say different beta diversity bc homo
# 
# #Srad
# dist.core.sr <- vegdist(seq.sr.core)
# samdf.core.sr <- data.frame(sample_data(ps.sr.rel))
# row.names(samdf.core.sr)==row.names(seq.sr.rel)
# set.seed(123)
# bet.core.sr <- betadisper(dist.core.sr,samdf.core.sr$site_zone)
# anova(bet.core.sr) #ps.cleanest ***
# plot(bet.core.sr) #SWB to the side and overlapping with SMB
# permutest(bet.core.sr, pairwise = TRUE, permutations = 999)
# #ps.cleanest,ps.trim SMB vs SWB ***
# adonis2(dist.core.sr ~ site_zone, data=samdf.core.sr, permutations=999)
# #ps.cleanest,ps.trim ***
# adonis2(dist.core.sr ~ site_zone/m_y, data=samdf.core.sr, permutations=999)
# #ps.cleanest s_z ***, s_z:m_y ***, ps.trim both ***
# #betadisper: NOT homogenous dispersion (***) and does not cluster nicely ps.cleanest
# #adonis: significant, but can't interpret anything about beta div bc not homo
# #looking at m_y_s_z -> march 2020 SMB clusters separately
# 
# #ppor
# dist.core.pp <- vegdist(seq.pp.core)
# samdf.core.pp <- data.frame(sample_data(ps.pp.rel))
# row.names(samdf.core.pp)==row.names(seq.pp.rel)
# set.seed(123)
# bet.core.pp <- betadisper(dist.core.pp,samdf.core.pp$site_zone)
# anova(bet.core.pp) #ps.cleanest,ps.trim ns
# plot(bet.core.pp) #very overlapped
# permutest(bet.core.pp, pairwise = TRUE, permutations = 999)
# #ps.cleanest,ps.trim ns
# adonis2(dist.core.pp ~ site_zone, data=samdf.core.pp, permutations=999)
# #ps.cleanest,ps.trim ns
# adonis2(dist.core.pp ~ site_zone/m_y, data=samdf.core.pp, permutations=999)
# #ps.cleanest,ps.trim ns both
# #betadisper: homogenous dispersion (ns)
# #adonis: ns, so no difference in beta diversity
# #looking at m_y_s_z -> march 2020 SMB clusters separately