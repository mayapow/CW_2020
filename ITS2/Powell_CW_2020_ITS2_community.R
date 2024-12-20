#Maya Powell
#3/17/2023
#ITS2 processing
#Curacao samples from March 2020, Nov 2020, Nov 2021

setwd('~/Documents/Castillo Lab/CW_2020/CW_2020_ITS2/ITS2')
#load libraries
library(phyloseq)
library('ggplot2')
library('Rmisc')
library(cowplot)
library("ggpubr")
library("vegan")
library(dplyr)
remotes::install_github("KarstensLab/microshades")
library("microshades")
remotes::install_github("mikemc/speedyseq")
library("speedyseq")
library("microViz")

#Check for low read samples
counts <- read.csv('symportal_profile_counts.csv',header=TRUE,row.names=1,check.names=FALSE)
plot(rowSums(counts)) 
counts$sum <- rowSums(counts)
print(counts$sum == "0")
#4 are zeros
counts.no0 <- subset(counts, sum != 0) 
print(counts.no0$sum == "0")
counts.no0 <- counts.no0 %>% select(-sum)
plot(rowSums(counts.no0)) 

#sample dataframe
samdf<-read.csv("sample_info_curacao_2020.csv")
head(samdf)
samdf <- samdf %>% select(-X)
rownames(samdf) <- samdf$id

#Make all phyloseq objects
#ran this once and then read them in later
#import taxa info
taxa <- read.csv("symportal_taxa.csv",header=TRUE)
rownames(taxa) <- as.factor(taxa$ITS2_type_profile)
mtaxa <- as.matrix(taxa)
# import counts (absolute abundance from its2 type profiles)
mcounts <- as.matrix(counts.no0)
# Construct phyloseq object 
ps <- phyloseq(otu_table(mcounts, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(mtaxa))
ps
#52 taxa and 154 samples
saveRDS(ps,file="ps.its2.RDS")

counts.all <- read.csv("symportal_allcounts.csv",header = T, row.names=1)
plot(rowSums(counts.all)) 
#1 zero, removing
counts.all$sum <- rowSums(counts.all)
print(counts.all$sum == "0")
#1 is zero
counts.all.no0 <- subset(counts.all, sum != 0) 
print(counts.all.no0$sum == "0")
counts.all.no0 <- counts.all.no0 %>% select(-sum)
plot(rowSums(counts.all.no0)) 

# import counts (absolute abundance from its2 post med)
mcounts.all <- as.matrix(counts.all.no0)

#
taxa <- read.csv("symportal_taxa_all.csv",header=TRUE)
rownames(taxa) <- as.factor(taxa$ITS2_type)
mtaxa <- as.matrix(taxa)
# Construct phyloseq object 
ps.all <- phyloseq(otu_table(mcounts.all, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(mtaxa))
ps.all
#528 taxa and 157 samples
saveRDS(ps.all,"ps.all.its2.RDS")
ps.all <- readRDS("ps.all.its2.RDS")

#read in ps objects
ps <- readRDS("ps.its2.RDS")
#add variable to make specific site labels 
ps.all <- ps.all %>% ps_mutate(site_nice =
                      case_when(site == "SMB" ~ "Santa Martha Bay", 
                                site == "SMR" ~ "Santa Martha Reef",
                                site == "SWB" ~ "Spaanse Water Bay", 
                                site == "DB" ~ "Spaanse Water Reef"))

ps.nd = subset_samples(ps.all, id!= "D6" & id!="N10" & id!="N11" & id!="N12" & id!= "N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#save ps.nd
saveRDS(ps.nd, "ps.nd.all.its2.RDS")
ps.nd <- readRDS("ps.nd.its2.RDS")
ps.rel <- transform_sample_counts(ps.nd, function(x) x / sum(x))

#ps.all <- readRDS("ps.all.its2.RDS")

##Bar plots
#relative abundance by sample
plot_bar(ps, "Sample", fill="ITS2_type_profile")

#absolute abundance
plot_bar(ps,"ITS2_type_profile", fill="ITS2_type_profile",facet_grid=host_species~site_zone)

##Type profile by sample##
#all samples
#plot_bar(ps,"sample_full")
ps.rel <- transform_sample_counts(ps.nd, function(x) x / sum(x))
gg.bar <- plot_bar(ps.rel,"sample_full",fill="Majority_ITS2_type")+
  geom_bar(stat="identity")+
  theme_classic()+
  facet_grid(~site_zone, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site")+
  ylab("Relative Abundance")
gg.bar = gg.bar + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar
ggsave(gg.bar,file="sym.barplot.all.pdf",h=8,w=35)    

#now do this for all datasets (relative) split by species
ps.ss.rel <- subset_samples(ps.rel, host_species=="siderea")
ps.sr.rel <- subset_samples(ps.rel, host_species=="radians")
ps.pp.rel <- subset_samples(ps.rel, host_species=="porites")
#ps.pa.rel <- subset_samples(ps.rel, host_species=="astreoides")

mdf.ss <- prep_mdf(ps.ss.rel)

#siderastrea siderea - these are the figures to put in the paper!!!!
#making supp plot of all types
#TOO HARD TO VISUALIZE - USE MICROSHADES
gg.bar.ss <- plot_bar(ps.ss.rel,"sample_full",fill="ITS2_type")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22, legend.position = "bottom")+
  facet_wrap(~site_zone, scales = "free")+
  ylab("Relative Abundance")
gg.bar.ss = gg.bar.ss + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.ss
ggsave(gg.bar.ss,file="sym.barplot.ss.site.pdf",h=10,w=15)
#full plot - to put in supplementary, with all data across sites and time
gg.bar.ss <- plot_bar(ps.ss.rel,"number",fill="Majority_ITS2_sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_nice, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.ss
#each number corresponds to specific sample here
ggsave(gg.bar.ss,file="sym.barplot.ss.site.season.sampleid.pdf",h=15,w=20)    

#siderastrea radians
gg.bar.sr <- plot_bar(ps.sr.rel,"sample_full",fill="Majority_ITS2_sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(~site_nice, scales = "free")+
  ylab("Relative Abundance")
gg.bar.sr = gg.bar.sr + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.sr
ggsave(gg.bar.sr,file="sym.barplot.sr.site.pdf",h=5,w=15)
#full plot - to put in supplementary, with all data across sites and time
gg.bar.sr <- plot_bar(ps.sr.rel,"number",fill="Majority_ITS2_sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_nice, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.sr
#each number corresponds to specific sample here, removed duplicate of D6 = sample ID 91 Nov 2020 SWB
ggsave(gg.bar.sr,file="sym.barplot.sr.site.season.sampleid.pdf",h=15,w=15)    

#porites porites
gg.bar.pp <- plot_bar(ps.pp.rel,"sample_full",fill="Majority_ITS2_sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(~site_nice, scales = "free")+
  ylab("Relative Abundance")
gg.bar.pp = gg.bar.pp + theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.pp
ggsave(gg.bar.pp,file="sym.barplot.pp.site.pdf",h=5,w=15)
#full plot - to put in supplementary, with all data acropp sites and time
gg.bar.pp <- plot_bar(ps.pp.rel,"number",fill="Majority_ITS2_sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_nice, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.pp
#each number corresponds to specific sample here #removed duplicates same as bacterial dataset (Ns)
ggsave(gg.bar.pp,file="sym.barplot.pp.site.season.sampleid.pdf",h=15,w=15)    

##Type profile by site and zone and timepoint
ps.sz <- merge_samples(ps, "site_zone")
ps.rel.sz <- transform_sample_counts(ps.sz, function(x) x / sum(x))
plot_bar(ps.rel.sz, fill="Majority_ITS2_sequence")

ps.my.sz <-merge_samples(ps, "m_y_s_z")
ps.rel.my.sz <- transform_sample_counts(ps.my.sz, function(x) x / sum(x))
plot_bar(ps.rel.my.sz, fill="ITS2_type_profile")

#sidereastrea siderea
ps.ss <- subset_samples(ps, host_species=="siderea")
ps.ss.z <- merge_samples(ps.ss, "site_zone")
ps.ss.z.rel <- transform_sample_counts(ps.ss.z, function(x) x / sum(x))
bar.ss <- plot_bar(ps.ss.z.rel, fill="Majority_ITS2_sequence")+
  theme_classic(base_size=22)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ylab("Relative abundance")+
  xlab("Site")
  #ggtitle("Siderastrea siderea")
bar.ss
ggsave(bar.ss,file="ssid.major.site.pdf",h=10,w=7) 

#siderastrea radians
ps.sr <- subset_samples(ps, host_species=="radians")
ps.sr.z <- merge_samples(ps.sr, "site_zone")
ps.sr.z.rel <- transform_sample_counts(ps.sr.z, function(x) x / sum(x))
bar.sr <- plot_bar(ps.sr.z.rel, fill="Majority_ITS2_sequence")+
  theme_classic(base_size = 22)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+
  ylab("Relative abundance")+
  xlab("Site")
  #ggtitle("Siderastrea radians")
bar.sr
ggsave(bar.sr,file="srad.major.site.pdf",h=10,w=4) 

library(viridis)
#porites porites
ps.pp <- subset_samples(ps, host_species=="porites")
ps.pp.z <- merge_samples(ps.pp, "site_zone")
ps.pp.z.rel <- transform_sample_counts(ps.pp.z, function(x) x / sum(x))
bar.pp <- plot_bar(ps.pp.z.rel, fill="Majority_ITS2_sequence")+
  theme_classic(base_size=22)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+
  ylab("Relative abundance")+
  xlab("Site")
  #ggtitle("Porites porites")+
bar.pp
ggsave(bar.pp,file="ppor.major.site.pdf",h=10,w=4) 

#put these graphs together
gg.panels.spp.site <- ggarrange(bar.ss,bar.sr,bar.pp,nrow=1,ncol=3,labels="AUTO")
gg.panels.spp.site
##ggsave(gg.panels.site,file="bac.div.site.pdf",height=8)


#porites astreoides
ps.pa <- subset_samples(ps, host_species=="astreoides")
ps.pa.z <- merge_samples(ps.pa, "m_y_s_z")
ps.pa.z.rel <- transform_sample_counts(ps.pa.z, function(x) x / sum(x))
bar.pa <- plot_bar(ps.pa.z.rel, fill="Majority_ITS2_sequence")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Relative abundance")+
  xlab("Site & Timepoint")+
  ggtitle("Porites astreoides")
bar.pa
ggsave(bar.pa,file="past.major.pdf",h=10,w=18) 

#tracking specific samples over time
#ssid
ps.ss.35 <- subset_samples(ps.ss, no_sp_lo=="97_siderea_SW_reef")

ps.ss.rel <- transform_sample_counts(ps.ss.35, function(x) x / sum(x))
gg.ss.bar.s <- plot_bar(ps.ss.rel,"full_sample_id",fill="ITS2_type_profile")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Sample ID")+
  ylab("Relative Abundance")
gg.ss.bar.s
ggsave(gg.ss.bar.s,file="gg.ss.bar.97SSSWR.pdf",h=10,w=14) 

#srad
ps.sr.sample <- subset_samples(ps.sr, no_sp_lo=="98_radians_SW_bay")

ps.sr.s.rel <- transform_sample_counts(ps.sr.sample, function(x) x / sum(x))
gg.sr.bar.s <- plot_bar(ps.sr.s.rel,"full_sample_id",fill="ITS2_type_profile")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Sample ID")+
  ylab("Relative Abundance")
gg.sr.bar.s
ggsave(gg.sr.bar.s,file="gg.sr.bar.9SRSMB.pdf",h=10,w=10) 

#ppor
ps.pp.sample <- subset_samples(ps.pp, no_sp_lo=="47_porites_SM_reef")

ps.pp.s.rel <- transform_sample_counts(ps.pp.sample, function(x) x / sum(x))
gg.pp.bar.s <- plot_bar(ps.pp.s.rel,"sample_full",fill="ITS2_type_profile")+
  geom_bar(stat="identity")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Sample ID")+
  ylab("Relative Abundance")
gg.pp.bar.s
ggsave(gg.pp.bar.s,file="gg.pp.bar.47PPSMR.pdf",h=10,w=14) 

###PCAs###

#type profiles by site

#siderastrea siderea
ps.ss.ord <- ordinate(ps.ss,"NMDS",distance="bray")
gg.ss.ord <- plot_ordination(ps.ss, ps.ss.ord, color ="reef_bay", shape="area")+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(linetype=area))+
  theme_cowplot()
gg.ss.ord
ggsave(gg.ss.ord,file="gg.ss.ord.pdf",h=10,w=10) 

#siderastrea radians
ps.sr.ord <- ordinate(ps.sr,"PCoA",distance="bray")
gg.sr.ord <- plot_ordination(ps.sr, ps.sr.ord, color ="reef_bay", shape="area")+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(linetype=area))+
  theme_cowplot()
gg.sr.ord
ggsave(gg.sr.ord,file="gg.sr.ord.pdf",h=10,w=10) 

#porites porites
ps.pp.ord <- ordinate(ps.pp,"PCoA",distance="bray")
gg.pp.ord <- plot_ordination(ps.pp, ps.pp.ord, color ="reef_bay", shape="area")+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(linetype=area))+
  theme_cowplot()
gg.pp.ord
ggsave(gg.pp.ord,file="gg.pp.ord.pdf",h=10,w=10) 


#porites astreoides
ps.pa.ord <- ordinate(ps.pa,"PCoA",distance="bray")
gg.pa.ord <- plot_ordination(ps.pa, ps.pa.ord, color ="reef_bay", shape="area")+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(linetype=area))+
  theme_cowplot()
gg.pa.ord
ggsave(gg.pa.ord,file="gg.pa.ord.pdf",h=10,w=10) 


###Stats###
library(vegan)
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(dplyr)
#BiocManager::install("edgeR")
library(edgeR)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#separate by species and only look at site - not timepoint

#siderastrea siderea
seq.ss.ps <- data.frame(ps.ss@otu_table)
samdf.ss.ps <- data.frame(ps.ss@sam_data)
dist.ss.ps <- vegdist(seq.ss.ps, method = "bray")
bet.ss.ps <- betadisper(dist.ss.ps,samdf.ss.ps$site)
anova(bet.ss.ps) #significant
set.seed(123)
permutest(bet.ss.ps,pairwise=TRUE,permutations=999)
plot(bet.ss.ps)

adonis2(seq.ss.ps ~ site_zone, data=samdf.ss.ps, permutations=999) #p<001***
pairwise.adonis2(seq.ss.ps ~ site_zone, data = samdf.ss.ps, permutations=999)
#               pairs   F.Model         R2 p.value p.adjusted
#1  SW_bay vs SW_reef 18.560955 0.34018751   0.001      0.001
#2   SW_bay vs SM_bay  4.861360 0.11342057   0.003      0.003
#3  SW_bay vs SM_reef 11.402678 0.22183042   0.001      0.001
#4  SW_reef vs SM_bay  6.496329 0.16041772   0.002      0.002
#5 SW_reef vs SM_reef  5.516496 0.13287479   0.001      0.001
#6  SM_bay vs SM_reef  3.704101 0.08881862   0.002      0.002

#siderastrea radians
seq.sr.ps <- data.frame(ps.sr@otu_table)
samdf.sr.ps <- data.frame(ps.sr@sam_data)
dist.sr.ps <- vegdist(seq.sr.ps)
bet.sr.ps <- betadisper(dist.sr.ps,samdf.sr.ps$site)
anova(bet.sr.ps) #ns
set.seed(123)
permutest(bet.sr.ps,pairwise=TRUE,permutations=999)
#       SMB   SWB
#SMB         0.189
#SWB 0.18582
plot(bet.sr.ps)
adonis2(seq.sr.ps ~ site, data=samdf.sr.ps, permutations=999) #ns
pairwise.adonis2(seq.sr.ps, factors=samdf.sr.ps$site, permutations=999) #
#pairs  F.Model         R2 p.value p.adjusted
#1 SMB vs SWB 1.696757 0.04623724   0.128      0.128


#porites porites
seq.pp.ps <- data.frame(ps.pp@otu_table)
samdf.pp.ps <- data.frame(ps.pp@sam_data)
dist.pp.ps <- vegdist(seq.pp.ps)
bet.pp.ps <- betadisper(dist.pp.ps,samdf.pp.ps$site)
anova(bet.pp.ps) #very sig!!
set.seed(123)
permutest(bet.pp.ps,pairwise=TRUE,permutations=999)
#         SMB   SMR
#SMB            0.001
#SMR 8.5674e-05   
plot(bet.pp.ps)#+theme(legend.margin=0.01)
adonis2(seq.pp.ps ~ site, data=samdf.pp.ps, permutations=999) #p<001***
pairwise.adonis(seq.pp.ps, factors=samdf.pp.ps$site, permutations=999)
#       pairs  F.Model        R2 p.value p.adjusted
#1 SMR vs SMB 11.18058 0.2650646   0.001      0.001

###Now look at timepoint as a factor & split them up by site

#siderastrea siderea
ps.ss.SWB <- subset_samples(ps.ss, site_zone=="SW_bay")
ps.ss.SWR <- subset_samples(ps.ss, site_zone=="SW_reef")
ps.ss.SMB <- subset_samples(ps.ss, site_zone=="SM_bay")
ps.ss.SMR <- subset_samples(ps.ss, site_zone=="SM_reef")

seq.ss.site.ps <- data.frame(ps.ss.SWB@otu_table)
samdf.ss.site.ps <- data.frame(ps.ss.SWB@sam_data)
dist.ss.site.ps <- vegdist(seq.ss.site.ps)
bet.ss.site.ps <- betadisper(dist.ss.site.ps,samdf.ss.site.ps$m_y)
anova(bet.ss.site.ps) #SWB (.), SWR (ns), SMB (ns), SMR (ns)
set.seed(123)
permutest(bet.ss.site.ps,pairwise=TRUE,permutations=999)

plot(bet.ss.site.ps)
adonis2(seq.ss.site.ps ~ m_y, data=samdf.ss.site.ps, permutations=999) #SWB (.) SWB (ns) SMB (ns) SMR (ns)
set.seed(123)
pairwise.adonis(seq.ss.site.ps, factors=samdf.ss.site.ps$m_y, permutations=999)

#siderastrea radians
ps.sr.SWB <- subset_samples(ps.sr, site_zone=="SW_bay")
ps.sr.SMB <- subset_samples(ps.sr, site_zone=="SM_bay")

seq.sr.site.ps <- data.frame(ps.sr.SMB@otu_table)
samdf.sr.site.ps <- data.frame(ps.sr.SMB@sam_data)
dist.sr.site.ps <- vegdist(seq.sr.site.ps)
bet.sr.site.ps <- betadisper(dist.sr.site.ps,samdf.sr.site.ps$m_y)
anova(bet.sr.site.ps) #SWB (*), SMB (ns)
set.seed(123)
permutest(bet.sr.site.ps,pairwise=TRUE,permutations=999)

plot(bet.sr.ps)
adonis2(seq.sr.site.ps ~ m_y, data=samdf.sr.site.ps, permutations=999) #SWB (ns) SMB ()
pairwise.adonis(seq.sr.site.ps, factors=samdf.sr.site.ps$m_y, permutations=999) #
#pairs  F.Model         R2 p.value p.adjusted
#1 SMB vs SWB 1.696757 0.04623724   0.128      0.128

#porites porites
ps.pp.SMR <- subset_samples(ps.pp, site_zone=="SM_reef")
ps.pp.SMB <- subset_samples(ps.pp, site_zone=="SM_bay")

seq.pp.site.ps <- data.frame(ps.pp.SMR@otu_table)
samdf.pp.site.ps <- data.frame(ps.pp.SMR@sam_data)
dist.pp.site.ps <- vegdist(seq.pp.site.ps)
bet.pp.site.ps <- betadisper(dist.pp.site.ps,samdf.pp.site.ps$m_y)
anova(bet.pp.site.ps) #SMR ns
set.seed(123)
permutest(bet.pp.site.ps,pairwise=TRUE,permutations=999)
plot(bet.pp.site.ps)
adonis2(seq.pp.site.ps ~ m_y, data=samdf.pp.site.ps, permutations=999) #SMR ns
pairwise.adonis(seq.pp.site.ps, factors=samdf.pp.site.ps$m_y, permutations=999)

#testing out microshades palette
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps.rel <- ps.rel %>% tax_mutate(Majority_ITS2_sequence = gsub('/','', as.character(Majority_ITS2_sequence)))

mdf_prep <- prep_mdf(ps.rel, subgroup_level = "Majority_ITS2_sequence")
# Create a color object for the specified data
color_obj_prep_test <- create_color_dfs(mdf_prep, selected_groups = c("D","C","B","A"), group_level = "Clade", subgroup_level = "Majority_ITS2_sequence", cvd = TRUE)
color_obj_prep_test$mdf$Majority_ITS2_sequence <- as.character(color_obj_prep_test$mdf$Majority_ITS2_sequence)
# Extract
mdf_test <- color_obj_prep_test$mdf
cdf_test <- color_obj_prep_test$cdf

#subset different groups out
mdf_ss <- mdf_test %>% filter(host_species == "siderea")
mdf_sr <- mdf_test %>% filter(host_species == "radians")
mdf_pp <- mdf_test %>% filter(host_species == "porites")

#now make graphs

#ssid
ss_shade <- plot_microshades(mdf_ss, cdf_test, group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_nice, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss_shade
ggsave(ss_shade,file="sym.shade.ss.site.pdf",h=10,w=15)    

#full plot - to put in supplementary, with all data across sites and time
#numbers are individual samples here
ss_shade_time <- plot_microshades(mdf_ss, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_nice, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss_shade_time
ggsave(ss_shade_time,file="sym.shade.ss.site.time.pdf",h=15,w=25)    

#srad
sr_shade <- plot_microshades(mdf_sr, cdf_test, group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_nice, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr_shade
ggsave(sr_shade,file="sym.shade.sr.site.pdf",h=5,w=15)    

#full plot - to put in supplementary, with all data acrosr sites and time
#numbers are individual samples here
sr_shade_time <- plot_microshades(mdf_sr, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_nice, scales = "free", ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr_shade_time
ggsave(sr_shade_time,file="sym.shade.sr.site.time.pdf",h=15,w=15)

#ppor
pp_shade <- plot_microshades(mdf_pp, cdf_test, group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_nice, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Porites spp."))))
pp_shade
ggsave(pp_shade,file="sym.shade.pp.site.pdf",h=5,w=15)

#full plot - to put in supplementary, with all data acropp sites and time
#numbers are individual samples here
pp_shade_time <- plot_microshades(mdf_pp, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
  theme_classic(base_size = 30)+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_nice, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Porites spp."))))
pp_shade_time
ggsave(pp_shade_time,file="sym.shade.pp.site.time.pdf",h=10,w=15)  

#put plots together
sym_pp_sr_plot <- ggarrange(sr_shade,pp_shade, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
sym_all_plot <- ggarrange(ss_shade,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
ggsave(sym_all_plot, file="sym_all_plot_microshades.pdf",width=25,height=12)

#put plots together
sym_all_time <- ggarrange(ss_shade_time,
                          ggarrange(sr_shade_time, pp_shade_time, ncol = 2, 
                                    labels = c("B", "C"), font.label = list(size = 50),
                                    legend = FALSE), 
                          #heights = c(1,2), widths = c(2,1), 
                          nrow = 2,
                          legend = "right", common.legend = TRUE, 
                          labels = c("A"),font.label = list(size = 50))
ggsave(sym_all_time, file="sym_all_time_plot_supp_microshades.pdf",width=25,height=30)

#SUPPLEMENTARY PLOTS OF ALL POST MED SEQUENCES
mdf_prep <- prep_mdf(ps.nd, subgroup_level = "ITS2_type")
# Create a color object for the specified data
color_obj_prep_test <- create_color_dfs(mdf_prep, selected_groups = c("D","C","B","A"), group_level = "Clade", subgroup_level = "ITS2_type", cvd = TRUE)
color_obj_prep_test$mdf$ITS2_type <- as.character(color_obj_prep_test$mdf$ITS2_type)
# Extract
mdf_test <- color_obj_prep_test$mdf
cdf_test <- color_obj_prep_test$cdf

#subset different groups out
mdf_ss <- mdf_test %>% filter(host_species == "siderea")
mdf_sr <- mdf_test %>% filter(host_species == "radians")
mdf_pp <- mdf_test %>% filter(host_species == "porites")

#now make graphs

#ssid
ss_shade <- plot_microshades(mdf_ss, cdf_test, group_label = "Clade - ITS2") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss_shade
ggsave(ss_shade,file="sym.shade.ss.all.site.pdf",h=10,w=15)    

#full plot - to put in supplementary, with all data across sites and time
#numbers are individual samples here
ss_shade_time <- plot_microshades(mdf_ss, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss_shade_time
ggsave(ss_shade_time,file="sym.shade.ss.all.site.time.pdf",h=15,w=25)    

#srad
sr_shade <- plot_microshades(mdf_sr, cdf_test, group_label = "Clade - ITS2 type") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr_shade
ggsave(sr_shade,file="sym.shade.sr.all.site.pdf",h=5,w=15)    

#full plot - to put in supplementary, with all data acrosr sites and time
#numbers are individual samples here
sr_shade_time <- plot_microshades(mdf_sr, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone, scales = "free", ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr_shade_time
ggsave(sr_shade_time,file="sym.shade.sr.all.site.time.pdf",h=15,w=15)    

#ppor
pp_shade <- plot_microshades(mdf_pp, cdf_test, group_label = "Clade - ITS2 type") + 
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone, scales = "free")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Porites sp."))))
pp_shade
ggsave(pp_shade,file="sym.shade.pp.all.site.pdf",h=5,w=15)    

#full plot - to put in supplementary, with all data acropp sites and time
#numbers are individual samples here
pp_shade_time <- plot_microshades(mdf_pp, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(substitute(paste(italic("Porites spp."))))
pp_shade_time
ggsave(pp_shade_time,file="sym.shade.pp.all.site.time.pdf",h=10,w=15)    

#put plots together
sym_pp_sr_plot <- ggarrange(sr_shade,pp_shade, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
sym_all_plot <- ggarrange(ss_shade,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
ggsave(sym_all_plot, file="sym_all_plot_ALL_postmed_microshades.pdf",width=25,height=12)

#put plots together
sym_all_time <- ggarrange(ss_shade_time,
                          ggarrange(sr_shade_time, pp_shade_time, ncol = 2, 
                                    labels = c("B", "C"), font.label = list(size = 50),
                                    legend = FALSE), 
                          #heights = c(1,2), widths = c(2,1), 
                          nrow = 2,
                          legend = "right", common.legend = TRUE, 
                          labels = c("A"),font.label = list(size = 50))
ggsave(sym_all_time, file="sym_all_time_plot_supp_ALL_postmed_microshades.pdf",width=25,height=30)

##########################################################
#### Summary of Dominant ITS2 Majority Types and DIVs ####
##########################################################

setwd('~/Documents/Castillo Lab/CW_2020/CW_2020_ITS2/ITS2')
#all info
ps.all <- readRDS("ps.all.its2.RDS")
ps.all.rel <- transform_sample_counts(ps.all, function(x) x / sum(x))
seq.all <- data.frame(ps.all.rel@otu_table)
#samdf.all <- data.frame(ps.all.rel@sam_data)
samdf.all <- read.csv("sample_info_curacao_2020.csv")
rownames(samdf.all) <- samdf.all$id

#just majority ITS2 types
ps <- readRDS("ps.its2.RDS")
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
seq <- data.frame(ps.rel@otu_table)
seq.test <- seq 
#making the test just to check to make sure it doesn't change the identities of the taxa bc I'm ~paranoid~ hehe
samdf <- data.frame(ps.rel@sam_data)
length(unique(rownames(samdf))) == nrow(samdf)

## First T0
# change column names to majority its2 sequence
sym_taxa = read.csv(file = "symportal_taxa.csv", header = TRUE) 

sym_taxa$Majority_ITS2_sequence #you need to do it this way because you need to add the 1s etc to not have duplicates
colnames(seq) =  c("A4.1","A4.2","A4.3","A4ch","A4.4",         
                   "A4.5","B19.1","B19.2","B19aq","B19.3",        
                   "B19ao","B19ap","B5","B19.4","B19.5",        
                   "B2","C46.1","C3.1","C42ef.C42eg","C3.2",
                   "C46.C3","C46.C1","C47a.1","C46.C42.2","C42eg.C42ef",
                   "C47a.2","C1.1","C1.2","C3.3","C1.3",
                   "C1.4","C47a.3","C45l","C1.5","C3.4",
                   "C45a","C1.6","C3.5","C1.7","C1.8",
                   "C3af","C1.9","C42.2.1","C46.2","C42.2.2",
                   "D1.1","D1.2","D1.3","D1.4","D1.5","D1.6","D1.7")

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.sums <- seq %>%
  mutate(A4_sum = rowSums(select(., starts_with("A4")))) %>%
  #mutate(A4ch_sum = A4ch) %>%
  mutate(B19_sum = rowSums(select(., starts_with("B19")))) %>%
  #mutate(B19aq_sum = B19aq) %>%
  #mutate(B19ao_sum = B19ao) %>%
  #mutate(B19ap_sum = B19ap) %>%
  mutate(B5_sum = B5) %>%
  mutate(B2_sum = B2) %>%
  mutate(C46_sum = rowSums(select(., starts_with("C46.")))) %>%
  mutate(C3_sum = rowSums(select(., starts_with("C3")))) %>%
  #mutate(C42ef.C42eg_sum = C42ef.C42eg) %>%
  #mutate(C46.C3_sum = C46.C3) %>%
  #mutate(C46.C1_sum = C46.C1) %>%
  mutate(C47a_sum = rowSums(select(., starts_with("C47a.")))) %>%
  #mutate(C46.C42.2_sum =C46.C42.2) %>%
  #mutate(C42eg.C42ef_sum = C42eg.C42ef) %>%
  mutate(C1_sum = rowSums(select(., starts_with("C1.")))) %>%
  mutate(C45_sum = rowSums(select(., starts_with("C45")))) %>%
  #mutate(C45l_sum = C45l) %>%
  #mutate(C45a_sum = C45a) %>%
  #mutate(C3af_sum = C3af) %>%
  mutate(C42_sum = rowSums(select(., starts_with("C42")))) %>%
  #mutate(C42.2_sum = rowSums(select(., starts_with("C42.2.")))) %>%
  mutate(D1_sum = rowSums(select(., starts_with("D1.")))) %>%
  rownames_to_column(var = "id") %>%
  select(id, contains("_sum"))

its2.sums.rel <- left_join(seqtab.rel.sums, samdf, by = "id")

#use seqtab rel sum to make graphs
#"A4","B19","B2","C46","C3","C47a","C1","C45","C42","D1"
taxa.maj <- data.frame(colnames(seqtab.rel.sums))
taxa.maj <- data.frame(taxa.maj[-c(1), ])
colnames(taxa.maj)[1] <- 'maj_its2'
taxa.maj$genus <- c("A","B","B","B","C","C","C","C","C","C","D")
rownames(taxa.maj) <- taxa.maj$maj_its2

taxa.maj$maj_its2 = as.factor(taxa.maj$maj_its2)
taxa.maj$genus = as.factor(taxa.maj$genus)

#seqtab rel sums
rownames(seqtab.rel.sums) <- seqtab.rel.sums$id
seqtab.rel.sums.maj <- seqtab.rel.sums %>% select(-id)
#col names of seqtab.rel.sums.maj match rownames of taxa.maj
#rownames of seqtab.rel.sums.maj match rownames of samdf

View(ps.its2.maj@otu_table)
#make ps object of majority type sums
ps.its2.maj.new <- phyloseq(sample_data(samdf.all),
                       otu_table(ps.its2.maj@otu_table,taxa_are_rows=FALSE),
                       tax_table(as.matrix(ps.its2.maj@tax_table)))
ps.its2.maj.new
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 11 taxa and 154 samples ]:
# sample_data() Sample Data:        [ 154 samples by 33 sample variables ]:
# tax_table()   Taxonomy Table:     [ 11 taxa by 2 taxonomic ranks ]:
# taxa are columns

plot_bar(ps.its2.maj, x="id",fill="maj_its2")+
  theme_classic()

its2_colors = c("A4_sum" = "#ffaabb", "C47a_sum" = "#99ddff","C42_sum" = "#225555",
                "B19_sum" =  "#aaaa00","B5_sum" = "#44bb99", "B2_sum" = "#225522",
                "C46_sum" = "#222255","C3_sum" = "#77aadd",
                  "C1_sum" = "#eedd88",
                     "C45_sum" = "#ee8860", "D1_sum" = "#994455")
#"A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1"

ps.its2.maj.nd = subset_samples(ps.its2.maj.new, id!= "D6" & id!="N10" & id!="N11" & id!="N12" & id!= "N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
saveRDS(ps.its2.maj.nd,file = "ps.its2.maj.nd")
ps.its2.maj.nd <- readRDS("ps.its2.maj.nd")
saveRDS(ps.its2.maj.new,file = "ps.its2.maj")
ps.its2.maj <- readRDS("ps.its2.maj")

ps.ss.maj <- subset_samples(ps.its2.maj.nd, host_species=="siderea")
ps.sr.maj <- subset_samples(ps.its2.maj.nd, host_species=="radians")
ps.pp.maj <- subset_samples(ps.its2.maj.nd, host_species=="porites")



###SIDERASTREA SIDEREA###
ss.site <- plot_bar(ps.ss.maj, x="id", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.site
#ggsave(ss.site,file="its2.ss.site.pdf",h=10,w=15)    

#full plot - to put in supplementary, with all data across sites and time
#numbers are individual samples here
ss.site.time <- plot_bar(ps.ss.maj, x = "number", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free") +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.site.time
#ggsave(ss.site.time,file="its2.ss.site.time.pdf",h=15,w=25)    

###Siderastrea radians
sr.site <- plot_bar(ps.sr.maj, x="id", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr.site
#ggsave(sr.site,file="its2.sr.site.pdf",h=5,w=15)   

#full plot - to put in supplementary, with all data acrosr sites and time
#numbers are individual samples here
sr.site.time <- plot_bar(ps.sr.maj, x = "number", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free", nrow = 3, ncol = 2) +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr.site.time
#ggsave(sr.site.time,file="its2.sr.site.time.pdf",h=15,w=15)  

#####Branching Porites sp.
pp.site <- plot_bar(ps.pp.maj, x="id", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle(substitute("Branching Porites sp."))
pp.site
#ggsave(pp.site,file="its2.pp.site.pdf",h=5,w=15)   

#full plot - to put in supplementary, with all data acropp sites and time
#numbers are individual samples here
pp.site.time <- plot_bar(ps.pp.maj, x = "number", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free", nrow = 2, ncol = 2) +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")+
  ggtitle(substitute("Branching Porites sp."))
pp.site.time
#ggsave(pp.site.time,file="its2.pp.site.time.pdf",h=10,w=15)  


###put plots together
sym_pp_sr_plot <- ggarrange(sr.site,pp.site, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
sym_all_plot <- ggarrange(ss.site,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
ggsave(sym_all_plot, file="its2.plot.all.pdf",width=25,height=12)

#put plots together
sym_all_time <- ggarrange(ss.site.time,
                          ggarrange(sr.site.time, pp.site.time, ncol = 2, 
                                    labels = c("B", "C"), font.label = list(size = 50)), 
                          #heights = c(1,2), widths = c(2,1), 
                          nrow = 2,
                          #legend = "right", common.legend = TRUE, 
                          labels = c("A"),font.label = list(size = 50))
ggsave(sym_all_time, file="its2.plot.time.all.pdf",width=25,height=30)

#convert to factors and numeric as needed
its2.sums.rel = its2.sums.rel %>%
  mutate_at(c(2:12), as.numeric)

its2.sums.rel = its2.sums.rel %>%
  mutate_at(c(13:35), as.factor)

# add in dominant and minor distinctions to use in other plots
its2.dom.rel <- its2.sums.rel %>%
  mutate(dominant_type = case_when(A4_sum >= 0.7 ~ "A4",
                                   B19_sum >= 0.7 ~ "B19",
                                   B5_sum >= 0.7 ~ "B5",
                                   B2_sum >= 0.7 ~ "B2",
                                   C46_sum >= 0.7 ~ "C46",
                                   C3_sum >= 0.7 ~ "C3",
                                   C42_sum >= 0.7 ~ "C42",
                                   C47a_sum >= 0.7 ~ "C47a",
                                   C1_sum >= 0.7 ~ "C1",
                                   C45_sum >= 0.7 ~ "C45",
                                   D1_sum >= 0.7 ~ "D1")) %>%
  mutate(minor_type = case_when(A4_sum < 0.7 & A4_sum > 0.0 ~ "A4",
                                B19_sum < 0.7 & B19_sum > 0.0 ~ "B19",
                                B5_sum < 0.7 & B5_sum > 0.0 ~ "B5",
                                B2_sum < 0.7 & B2_sum > 0.0 ~ "B2",
                                C46_sum < 0.7 & C46_sum > 0.0 ~ "C46",
                                C3_sum < 0.7 & C3_sum > 0.0 ~ "C3",
                                C42_sum < 0.7 & C42_sum > 0.0 ~ "C42",
                                C47a_sum < 0.7 & C47a_sum > 0.0 ~ "C47a",
                                C1_sum < 0.7 & C1_sum > 0.0 ~ "C1",
                                C45_sum < 0.7 & C45_sum > 0.0 ~ "C45",
                                D1_sum < 0.7 & D1_sum > 0.0 ~ "D1"))

write.csv(its2.dom.rel, file = "ITS2.dominanttype.CW_2020.csv", row.names = FALSE)

its2.rel.dom <- filter(its2.dom.rel, host_species!="astreoides")
ss.its2.rel.dom <- filter(its2.dom.rel, host_species=="siderea")
sr.its2.rel.dom <- filter(its2.dom.rel, host_species=="radians")
pp.its2.rel.dom <- filter(its2.dom.rel, host_species=="porites")

# summarize proportion of individuals with more than 70% of the sym types
its2.dom.spp <- its2.dom.rel %>%
  group_by(host_species) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            #n_maj_B19 = sum(B19_sum >= 0.7),
            #n_maj_B5 = sum(B5_sum >= 0.7),
            #n_maj_B2 = sum(B2_sum >= 0.7),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            #n_maj_C45 = sum(C45_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            #p_maj_B19 = n_maj_B19/n,
            #p_maj_B5 = n_maj_B5/n,
            #p_maj_B2 = n_maj_B2/n,
            p_maj_C46 = n_maj_C46/n,
            p_maj_C3 = n_maj_C3/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n,
            p_maj_C1 = n_maj_C1/n,
            #p_maj_C45 = n_maj_C45/n,
            p_maj_D1 = n_maj_D1/n)
#based on this - removing samples from the species that do not have any majority of certain sym types
#pa not doing just placeholder for my braaaain
#pp keep = A4, C42, C47a
#sr keep = C46, C1, D1
#ss keep = C3, C1, D1
#all remove = B19, B5, B2, C45

####SSID######

#timepoint
ss.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   m_y               n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>         <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 March 2020       24        6        3       12    0.25    0.125     0.5  
# 2 November 2020    24        6        1       13    0.25    0.0417    0.542
# 3 November 2021    30        7        3       18    0.233   0.1       0.6  

#reef vs bay
ss.its2.rel.dom %>%
  group_by(reef_bay) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   reef_bay     n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>    <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 bay         40        0        3       33      0      0.075    0.825
# 2 reef        38       19        4       10      0.5    0.105    0.263

#site
ss.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   site_zone     n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>     <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 SM_bay       22        0        2       17    0       0.0909    0.773
# 2 SM_reef      21        8        4        6    0.381   0.190     0.286
# 3 SW_bay       18        0        1       16    0       0.0556    0.889
# 4 SW_reef      17       11        0        4    0.647   0         0.235

########SRAD######

#timepoint
sr.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C46 = n_maj_C46/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   m_y               n n_maj_C46 n_maj_C1 n_maj_D1 p_maj_C46 p_maj_C1 p_maj_D1
#   <fct>         <int>     <int>    <int>    <int>     <dbl>    <dbl>    <dbl>
# 1 March 2020       12        12        0        0     1        0       0     
# 2 November 2020    15        11        2        1     0.733    0.133   0.0667
# 3 November 2021    10        10        0        0     1        0       0     

#site
sr.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C46 = n_maj_C46/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   site_zone     n n_maj_C46 n_maj_C1 n_maj_D1 p_maj_C46 p_maj_C1 p_maj_D1
#   <fct>     <int>     <int>    <int>    <int>     <dbl>    <dbl>    <dbl>
# 1 SM_bay       18        16        2        0     0.889    0.111   0     
# 2 SW_bay       19        17        0        1     0.895    0       0.0526

####PPOR######

#timepoint
pp.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n)
#   m_y               n n_maj_A4 n_maj_C42 n_maj_C47a p_maj_A4 p_maj_C42 p_maj_C47a
#   <fct>         <int>    <int>     <int>      <int>    <dbl>     <dbl>      <dbl>
# 1 November 2020    11        6         2          1    0.545     0.182     0.0909
# 2 November 2021    22       10         4          8    0.455     0.182     0.364 

#site
pp.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n)
#   site_zone     n n_maj_A4 n_maj_C42 n_maj_C47a p_maj_A4 p_maj_C42 p_maj_C47a
#   <fct>     <int>    <int>     <int>      <int>    <dbl>     <dbl>      <dbl>
# 1 SM_bay       18       16         2          0    0.889     0.111        0  
# 2 SM_reef      15        0         4          9    0         0.267        0.6


#####STATS#####
# Kruskal-Wallis (non-parametric alternative to one-way anova) to see if there is a difference
# in mean proportion D based on timepoint, reef vs bay, and site

##CORAL SPECIES####
kruskal.test(C47a_sum ~ host_species, data = its2.rel.dom)
dunnTest(C47a_sum ~ host_species, data = its2.rel.dom, method="bonferroni")
# data:  C3_sum by host_species
# Kruskal-Wallis chi-squared = 30.343, df = 2, p-value = 2.577e-07
# Comparison         Z      P.unadj        P.adj
# 1 porites - radians  0.000000 1.000000e+00 1.000000e+00
# 2 porites - siderea -4.367233 1.258306e-05 3.774918e-05
# 3 radians - siderea -4.543210 5.540404e-06 1.662121e-05

#data:  D1_sum by host_species
#Kruskal-Wallis chi-squared = 87.675, df = 2, p-value < 2.2e-16
#         Comparison           Z      P.unadj        P.adj
# 1 porites - radians  0.01162176 9.907274e-01 1.000000e+00
# 2 porites - siderea -7.41652703 1.202314e-13 3.606941e-13
# 3 radians - siderea -7.72931509 1.081267e-14 3.243802e-14

# data:  C1_sum by host_species
# Kruskal-Wallis chi-squared = 8.8846, df = 2, p-value = 0.01177
# Comparison          Z     P.unadj      P.adj
# 1 porites - radians -0.7596603 0.447457678 1.00000000
# 2 porites - siderea -2.7481119 0.005993955 0.01798186
# 3 radians - siderea -1.9476548 0.051456277 0.15436883

# data:  C46_sum by host_species
# Kruskal-Wallis chi-squared = 110, df = 2, p-value < 2.2e-16
# Comparison         Z      P.unadj        P.adj
# 1 porites - radians -8.092023 5.868194e-16 1.760458e-15
# 2 porites - siderea  0.358183 7.202064e-01 1.000000e+00
# 3 radians - siderea 10.078778 6.857379e-24 2.057214e-23

# data:  A4_sum by host_species
# Kruskal-Wallis chi-squared = 53.756, df = 2, p-value = 2.124e-12
# Comparison          Z      P.unadj        P.adj
# 1 porites - radians  6.4502792 1.116443e-10 3.349329e-10
# 2 porites - siderea  6.6946068 2.162523e-11 6.487570e-11
# 3 radians - siderea -0.7725698 4.397771e-01 1.000000e+00

# data:  C42_sum by host_species
# Kruskal-Wallis chi-squared = 33.21, df = 2, p-value = 6.145e-08
# Comparison         Z      P.unadj        P.adj
# 1 porites - radians 4.3334275 1.468056e-05 4.404169e-05
# 2 porites - siderea 5.6528673 1.577932e-08 4.733795e-08
# 3 radians - siderea 0.6828201 4.947205e-01 1.000000e+00

# data:  C47a_sum by host_species
# Kruskal-Wallis chi-squared = 32.669, df = 2, p-value = 8.053e-08
# Comparison          Z      P.unadj        P.adj
# 1 porites - radians  4.8390231 1.304789e-06 3.914367e-06
# 2 porites - siderea  5.3608265 8.284206e-08 2.485262e-07
# 3 radians - siderea -0.2274366 8.200843e-01 1.000000e+00

####SSID####
library(FSA)
#sym types are: D1, C1, and C3
kruskal.test(D1_sum ~ site_zone, data = ss.its2.rel.dom)
# data:  D1_sum by m_y
# Kruskal-Wallis chi-squared = 0.44807, df = 2, p-value = 0.7993
# data:  D1_sum by reef_bay
# Kruskal-Wallis chi-squared = 21.166, df = 1, p-value = 4.212e-06
# data:  D1_sum by site_zone
# Kruskal-Wallis chi-squared = 23.323, df = 3, p-value = 3.458e-05
dunnTest(D1_sum ~ site_zone, data = ss.its2.rel.dom, method="bonferroni")
#          Comparison          Z      P.unadj        P.adj
# 1  SM_bay - SM_reef  2.8629102 4.197695e-03 0.0251861676
# 2   SM_bay - SW_bay -1.4455974 1.482902e-01 0.8897409920
# 3  SM_reef - SW_bay -4.1495122 3.331846e-05 0.0001999107
# 4  SM_bay - SW_reef  2.4418297 1.461304e-02 0.0876782289
# 5 SM_reef - SW_reef -0.2602169 7.946964e-01 1.0000000000
# 6  SW_bay - SW_reef  3.6900069 2.242480e-04 0.0013454879

kruskal.test(C1_sum ~ site_zone, data = ss.its2.rel.dom)
# data:  C1_sum by m_y
# Kruskal-Wallis chi-squared = 0.70146, df = 2, p-value = 0.7042
# data:  C1_sum by reef_bay
# Kruskal-Wallis chi-squared = 0.92886, df = 1, p-value = 0.3352
# data:  C1_sum by site_zone
# Kruskal-Wallis chi-squared = 4.7558, df = 3, p-value = 0.1906
kruskal.test(C3_sum ~ site_zone, data = ss.its2.rel.dom)
# data:  C3_sum by m_y
# Kruskal-Wallis chi-squared = 0.10039, df = 2, p-value = 0.951
# data:  C3_sum by reef_bay
# Kruskal-Wallis chi-squared = 35.042, df = 1, p-value = 3.226e-09
# data:  C3_sum by site_zone
# Kruskal-Wallis chi-squared = 37.907, df = 3, p-value = 2.958e-08
dunnTest(C3_sum ~ site_zone, data = ss.its2.rel.dom, method="bonferroni")
#          Comparison          Z      P.unadj        P.adj
# 1  SM_bay - SM_reef -3.4273855 6.094232e-04 3.656539e-03
# 2   SM_bay - SW_bay  0.3825644 7.020427e-01 1.000000e+00
# 3  SM_reef - SW_bay  3.6338211 2.792546e-04 1.675527e-03
# 4  SM_bay - SW_reef -4.9037320 9.403274e-07 5.641965e-06
# 5 SM_reef - SW_reef -1.6486878 9.921162e-02 5.952697e-01
# 6  SW_bay - SW_reef -5.0417118 4.613857e-07 2.768314e-06

####SRAD####
#sym types are: C46, C1, D1
kruskal.test(D1_sum ~ m_y, data = sr.its2.rel.dom)
# data:  D1_sum by m_y
# Kruskal-Wallis chi-squared = 1.4667, df = 2, p-value = 0.4803
# data:  D1_sum by site_zone
# Kruskal-Wallis chi-squared = 0.94737, df = 1, p-value = 0.3304
kruskal.test(C1_sum ~ site_zone, data = sr.its2.rel.dom)
# data:  C1_sum by m_y
# Kruskal-Wallis chi-squared = 3.0148, df = 2, p-value = 0.2215
# data:  C1_sum by site_zone
# Kruskal-Wallis chi-squared = 2.1698, df = 1, p-value = 0.1407
kruskal.test(C46_sum ~ m_y, data = sr.its2.rel.dom)
# data:  C46_sum by site_zone
# Kruskal-Wallis chi-squared = 1.052, df = 1, p-value = 0.305
# data:  C46_sum by m_y
# Kruskal-Wallis chi-squared = 6.9282, df = 2, p-value = 0.0313
dunnTest(C46_sum ~ m_y, data = sr.its2.rel.dom, method="bonferroni")
#                      Comparison          Z    P.unadj      P.adj
# 1    March 2020 - November 2020  2.4266191 0.01524025 0.04572074
# 2    March 2020 - November 2021  0.3378417 0.73548248 1.00000000
# 3 November 2020 - November 2021 -1.9477617 0.05144349 0.15433046

####PPOR######
#sym types are A4, C42, C47a
kruskal.test(A4_sum ~ site_zone, data = pp.its2.rel.dom)
# data:  A4_sum by m_y
# Kruskal-Wallis chi-squared = 0.56796, df = 1, p-value = 0.4511
# data:  A4_sum by site_zone
# Kruskal-Wallis chi-squared = 20.711, df = 1, p-value = 5.342e-06
kruskal.test(C42_sum ~ m_y, data = pp.its2.rel.dom)
# data:  C42_sum by m_y
# Kruskal-Wallis chi-squared = 0.28407, df = 1, p-value = 0.594
# data:  C42_sum by site_zone
# Kruskal-Wallis chi-squared = 0.015913, df = 1, p-value = 0.8996
kruskal.test(C47a_sum ~ site_zone, data = pp.its2.rel.dom)
# data:  C47a_sum by m_y
# Kruskal-Wallis chi-squared = 1.5364, df = 1, p-value = 0.2152
# data:  C47a_sum by site_zone
# Kruskal-Wallis chi-squared = 16.5, df = 1, p-value = 4.865e-05

