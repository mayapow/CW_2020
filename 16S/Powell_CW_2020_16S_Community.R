#Maya Powell
#Dec 2025 editing for manuscript revisions
#16S COMMUNITY COMPOSITION
#Based on scripts from many people online (see packages for reference if you want lol)
#and also some real live people: Ana Dulskiy (UNC Castillo lab), Nicola Kreifall (BU Davies Lab), Hannah Aichelman (BU Davies Lab) and Steph Smith (UNC Septer Lab)

#load libraries
#install.packages("BiocManager") #use BiocManager to install packages if not available through CRAN
#BiocManager::install("phyloseq") #replace with any packages to install as needed
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
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#library(microbiomeutilities)
#devtools::install_github("david-barnett/microViz")
library("microViz")
library("here")
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)

#####INITIAL DATA PROCESSING####
#commented out because only needed intitially :)
#skip unless you need to remake phyloseq objects or examine another type (e.g. rarefied, trimmed, etc)
#load(here("16S","CW_2020_16S_taxa2.Rdata"))
#load in Phyloseq objects
#ps.cleanest.nd <- readRDS(here("16S","CW_2020_16S_ps.cleanest.nd.RDS"))
#ps.cleanest.nd
#27653 taxa and 134 samples
#ps.cleanest #27653 taxa and 144 samples, RAW DATA untrimmed, unrarefied
#remove duplicates from lane2 - they are very close (not stat sig dif) and lane 2 gets trimmed out earlier
#ps.cleanest.nd = subset_samples(ps.cleanest, id!="N10" & id!="N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.cleanest.nd <- data.frame(sample_data(ps.cleanest.nd))

#ps.less.nd <- readRDS(here("16S","CW_2020_16S_ps.less.nd.RDS"))
#ps.less #27653 taxa and 135 samples, samples witN counts <1,000 removed (9 samples removed)
#ps.less.nd = subset_samples(ps.less.nd, id!="N10" & id!="N11" & id!="N1" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.less.nd <- data.frame(sample_data(ps.less.nd))
#ps.less.nd

#ps.rare.nd <- readRDS(here("16S","CW_2020_16S_ps.rare.nd.RDS"))
#ps.rare #23759 taxa and 98 samples, rarefied to 10,000 (46 samples removed)
#ps.rare.nd = subset_samples(ps.rare, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.rare.nd <- data.frame(sample_data(ps.rare.nd))
#ps.rare.nd
#23759 taxa and 94 samples

#ps.trim.nd <- readRDS(here("16S","CW_2020_16S_ps.trim.nd.RDS"))
#ps.trim #835 taxa and 134 samples, trimmed witN MCMC OTU for low abundance taxa (10 samples removed)
#ps.trim.nd = subset_samples(ps.trim, id!="N10" & id!="N11" & id!="N2" & id!="N3" & id!="N4" & id!="N5" & id!="N6" & id!="N7" & id!="N8" & id!="N9")
#samdf.trim.nd <- data.frame(sample_data(ps.trim.nd))
#ps.trim.nd
#835 taxa and 128 samples

#Rename ASVs to be more informative
#go through and do this with all different dataframes, as needed
#only needed VERY initially - all phyloseq objects have this already!
# tax <- as.data.frame(ps.trim.rare.nd@tax_table@.Data)
# tax.clean <- data.frame(row.names = row.names(tax),
#                         Kingdom = str_replace(tax[,1], "D_0__",""),
#                         Phylum = str_replace(tax[,2], "D_1__",""),
#                         Class = str_replace(tax[,3], "D_2__",""),
#                         Order = str_replace(tax[,4], "D_3__",""),
#                         Family = str_replace(tax[,5], "D_4__",""),
#                         Genus = str_replace(tax[,6], "D_5__",""),
#                         Species = str_replace(tax[,7], "D_6__",""),
#                         stringsAsFactors = FALSE)
# tax.clean[is.na(tax.clean)] <- ""
# for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
# ####### Fill holes in the tax table
# tax.clean[is.na(tax.clean)] <- ""
# for (i in 1:nrow(tax.clean)){
#   if (tax.clean[i,2] == ""){
#     kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
#     tax.clean[i, 2:7] <- kingdom
#   } else if (tax.clean[i,3] == ""){
#     phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
#     tax.clean[i, 3:7] <- phylum
#   } else if (tax.clean[i,4] == ""){
#     class <- paste("Class_", tax.clean[i,3], sep = "")
#     tax.clean[i, 4:7] <- class
#   } else if (tax.clean[i,5] == ""){
#     order <- paste("Order_", tax.clean[i,4], sep = "")
#     tax.clean[i, 5:7] <- order
#   } else if (tax.clean[i,6] == ""){
#     family <- paste("Family_", tax.clean[i,5], sep = "")
#     tax.clean[i, 6:7] <- family
#   } else if (tax.clean[i,7] == ""){
#     tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
#   }
# }
# tax_table(ps.trim.rare.nd) <- as.matrix(tax.clean)
# #do this for each one, then save and subset below, then only need to do step above once

#adding other variables as needed - then resave
# #adding in site_nice column
# ps.less.nd <- ps.less.nd %>% ps_mutate(site_nice =
#                                          case_when(site == "SMB" ~ "Santa Martha Bay", 
#                                                    site == "SMR" ~ "Santa Martha Reef",
#                                                    site == "SWB" ~ "Spaanse Water Bay", 
#                                                    site == "DB" ~ "Spaanse Water Reef"))
# 
# #adding in species and genus column 
# ps.less.nd <- ps.less.nd %>% ps_mutate(coral_species =
#                                          case_when(host_species == "siderea" ~ "Siderastrea siderea", 
#                                                    host_species == "radians" ~ "Siderastrea radians",
#                                                    host_species == "porites" ~ "Porites sp."))
# 
# 
# ps.less.nd <- ps.less.nd %>% ps_mutate(tp_suffix = 
#                                          case_when(
#                                            m_y == "March_2020"    ~ "a",
#                                            m_y == "November_2020" ~ "b",
#                                            m_y == "November_2021" ~ "c")) %>%
#   ps_mutate(number_time = paste0(number,tp_suffix))
# 
# #add column for Habitat:Timepoint
# ps.less.nd@sam_data$m_rb <- paste(ps.less.nd@sam_data$m_y,ps.less.nd@sam_data$reef_bay)
# ps.less.nd <- subset_samples(ps.less.nd,host_species!="astreoides")
# 
# #View(data.frame(ps.less.nd@sam_data))
# 

#resaving and reading in nd (no duplicates) dataframes with cleaned up tax tables
#can always redo from beginning
# saveRDS(ps.cleanest.nd, "CW_2020_16S_ps.cleanest.nd.RDS")
# saveRDS(ps.less.nd, here("16S","CW_2020_16S_ps.less.nd.RDS"))
# saveRDS(ps.rare.nd, "CW_2020_16S_ps.rare.nd.RDS")
# saveRDS(ps.trim.nd, "CW_2020_16S_ps.trim.nd.RDS")
# saveRDS(ps.trim.rare.nd, "CW_2020_16S_ps.trim.rare.nd.RDS")

#coming back to this script - START HERE and just read in phyloseq objects
#these have cleaned up tax tables and no duplicate samples
#and now these have all the environmental data added! MP 4/20/2023
#ps.cleanest.nd <- readRDS(here("16S","CW_2020_16S_ps.cleanest.nd.RDS"))
ps.less.nd <- readRDS(here("16S","CW_2020_16S_ps.less.nd.RDS"))
#ps.rare.nd <- readRDS(here("16S","CW_2020_16S_ps.rare.nd.RDS"))
#ps.trim.nd <- readRDS(here("16S","CW_2020_16S_ps.trim.nd.RDS"))
#ps.trim.rare.nd <- readRDS(here("16S","CW_2020_16S_ps.trim.rare.nd.RDS"))

#split by species for each:
#do this for each type
ps.cl.ss <- subset_samples(ps.less.nd,host_species=="siderea")
ps.cl.sr <- subset_samples(ps.less.nd,host_species=="radians")
ps.cl.pp <- subset_samples(ps.less.nd,host_species=="porites")

#make relative phyloseq objects
ps.all.rel <- transform_sample_counts(ps.less.nd, function(x) x / sum(x))
ps.ss.rel <- transform_sample_counts(ps.cl.ss, function(x) x / sum(x))
ps.sr.rel <- transform_sample_counts(ps.cl.sr, function(x) x / sum(x))
ps.pp.rel <- transform_sample_counts(ps.cl.pp, function(x) x / sum(x))

#make as needed (take a while to make)
seq.all.ps.less <- data.frame(otu_table(ps.less.nd))
seq.ss <- data.frame(otu_table(ps.cl.ss))
seq.sr <- data.frame(otu_table(ps.cl.sr))
seq.pp <- data.frame(otu_table(ps.cl.pp))

#write.csv(seq.all.ps.rare, "seq.all.ps.rare.csv")
#seq.all.ps.rare.t <- as.data.frame(t(seq.all.ps.rare))
#seq.all.ps.rare.t <- seq.all.ps.rare.t %>%
#  mutate(ids=row.names(seq.all.ps.rare.t),
#         .before=A1)
#seq.all.ps.rare.t
#write.table(seq.all.ps.rare.t, file=here("16S","seq.all.ps.rare.t.txt", sep="\t", row.names = FALSE)
#rownames(seq.all.ps.rare.t) <- NULL

#####ORDINATION PLOTS#####
#using PCoA's and bray curtis distances

#sample dataframe as needed
samdf.all <- data.frame(sample_data(ps.all.rel))

#coral species
spp.ord <- plot_ordination(ps.all.rel,ordinate(ps.all.rel,"PCoA", "bray"),color="host_species")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  theme(legend.position = "right")+
  geom_point(size=5)+
  scale_color_manual("Coral Species", values = c("darkkhaki","darksalmon","indianred"), labels = c(substitute(paste(italic("Porites sp."))), substitute(paste(italic("S. radians"))),substitute(paste(italic("S. siderea")))))
spp.ord
#ggsave(spp.ord, file=here("16S","spp.ord.ps.less.pdf"),w=10,h=8)

#species dominant symbiont type
its2_colors = c("A4" = "#ffaabb",
                "C1"="#eedd88",
                "C3"="#77aadd",  
                "C42"="#225555",
                "C45"="#ee8860",
                "C46"="#222255",
                "C47a"="#99ddff",
                "D1"="#994455")
spp.its2.ord <- plot_ordination(ps.all.rel,ordinate(ps.all.rel,"PCoA", "bray"),color="dominant_type")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  theme(legend.position = "right")+
  geom_point(size=5)+
  scale_color_manual("Majority ITS2 Type", values = its2_colors)
spp.its2.ord

#ssid
ss.ord <- plot_ordination(ps.ss.rel,ordinate(ps.ss.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"),labels = c("Santa Martha Bay", "Santa Martha Reef", "Spaanse Water Bay", "Spaanse Water Reef"))
ss.ord
#ggsave(ss.ord, file=here("16S","ss.ord.ps.less.pdf"),w=10,h=8)

#ssid timepoint
ss.time.ord <- plot_ordination(ps.ss.rel,ordinate(ps.ss.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual(name = "Timepoint", values = c("goldenrod","cornflowerblue", "skyblue"),labels = c("March 2020", "November 2020", "November 2021"))
ss.time.ord
#ggsave(ss.time.ord, file=here("16S","ss.ord.time.ps.less.pdf"),w=10,h=8)

#ss dominant symbiont type
dom.ss.ord <- plot_ordination(ps.ss.rel,ordinate(ps.ss.rel,"PCoA", "bray"),color="dominant_type")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  theme(legend.position = "right")+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))+
  scale_color_manual("Majority ITS2 Type", values = its2_colors)
dom.ss.ord

#srad
sr.ord <- plot_ordination(ps.sr.rel,ordinate(ps.sr.rel,"PCoA", "bray"),color="site_zone")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  scale_color_manual(name = "Site", values = c("coral3","lightsalmon1"),labels = c("Santa Martha Bay", "Spaanse Water Bay"))
sr.ord
#ggsave(sr.ord, file=here("16S","sr.ord.ps.less.pdf"),w=10,h=8)

#srad timepoint
sr.time.ord <- plot_ordination(ps.sr.rel,ordinate(ps.sr.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  scale_color_manual(name = "Timepoint", values = c("goldenrod","cornflowerblue", "skyblue"),labels = c("March 2020", "November 2020", "November 2021"))
sr.time.ord
#gsave(sr.time.ord, file=here("16S","sr.ord.time.ps.less.pdf"),w=10,h=8)

#srad dominant symbiont type
dom.sr.ord <- plot_ordination(ps.sr.rel,ordinate(ps.sr.rel,"PCoA", "bray"),color="dominant_type")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  theme(legend.position = "right")+
  geom_point(size=5)+
  ggtitle(substitute(paste(italic("Siderastrea radians"))))+
  scale_color_manual("Majority ITS2 Type", values = its2_colors)
dom.sr.ord

#ppor
pp.ord <- plot_ordination(ps.pp.rel,ordinate(ps.pp.rel,"PCoA", "bray"),color="site_zone")+
  #stat_ellipse(aes(linetype=reef_bay),linewidth=1)+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))+
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise"),labels = c("Santa Martha Bay", "Santa Martha Reef"))
pp.ord
#ggsave(pp.ord, file=here("16S","pp.ord.ps.less.pdf"),w=10,h=8)

#ppor timepoint
pp.time.ord <- plot_ordination(ps.pp.rel,ordinate(ps.pp.rel,"PCoA", "bray"),color="m_y")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  geom_point(size=5)+
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))+
  scale_color_manual(name = "Timepoint", values = c("cornflowerblue", "skyblue"),labels = c("November 2020", "November 2021"))
pp.time.ord
#ggsave(pp.time.ord, file=here("16S","pp.ord.time.ps.less.pdf"),w=10,h=8)

#ppor dominant symbiont type
dom.pp.ord <- plot_ordination(ps.pp.rel,ordinate(ps.pp.rel,"PCoA", "bray"),color="dominant_type")+
  stat_ellipse(linewidth=2)+
  theme_classic(base_size = 22)+
  theme(legend.position = "right")+
  geom_point(size=5)+
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))+
  scale_color_manual("Majority ITS2 Type", values = its2_colors)
dom.pp.ord

pcoa.all <- ggarrange(ss.ord,sr.ord,pp.ord,nrow=1,ncol=3,legend="right",common.legend = TRUE, labels = c("A","B","C"), font.label = list(size = 30))
pcoa.time <- ggarrange(ss.time.ord,sr.time.ord,pp.time.ord,nrow=1,ncol=3,legend="right",common.legend = TRUE, labels = c("D","E","F"), font.label = list(size = 30))
pcoa.sym1 <- ggarrange(dom.ss.ord,dom.sr.ord,dom.pp.ord,nrow=1,ncol=3,legend="none",common.legend = TRUE, labels = c("G","H","I"), font.label = list(size = 30))
legend_col <- get_legend(spp.its2.ord) #use combined ITS2 legend here
pcoa.sym <- ggarrange(pcoa.sym1, legend_col, nrow = 1, ncol = 2, widths = c(1, 0.18))

pcoa.together <- ggarrange(pcoa.all, pcoa.time, pcoa.sym, nrow=3,ncol=1)
ggsave(pcoa.together, file=here("16S","pcoa.site.time.sym.all.pdf"),width=20,height=20)

###BETA DIVERSITY STATS####
#ALL
#seq.all.ps.less <- data.frame(otu_table(ps.less.nd)) #uncomment as needed if not made already above
samdf.all <- data.frame(sample_data(ps.all.rel))

dist.all <- vegdist(seq.all.ps.less)
row.names(samdf.all)==row.names(seq.all.ps.less)
bet.all.spp <- betadisper(dist.all,samdf.all$host_species)

anova(bet.all.spp)
plot(bet.all.spp)

permutest(bet.all.spp, pairwise = TRUE, permutations = 999) 

adonis2(dist.all ~ host_species, data=samdf.all, permutations=999)
pairwise.adonis2(dist.all ~ host_species, data = samdf.all, permutations=999)
#all species different from eachother

#ssid
dist.ss <- vegdist(seq.ss)
samdf.ss <- subset(samdf.all, host_species == "siderea")
row.names(samdf.ss)==row.names(seq.ss)
bet.ss <- betadisper(dist.ss,samdf.ss$dominant_type)
anova(bet.ss) #ps.less ns
#dom = n.s.
plot(bet.ss)
#permutest(bet.ss, pairwise = TRUE, permutations = 999)

adonis2(dist.ss ~ dominant_type,data=samdf.ss, permutations=999) #sig!
#dom = 0.012 *

#go through each iteration of pairwise comparisons bc intxns don't work
pairwise.adonis2(dist.ss ~ dominant_type, data=samdf.ss, permutations=999)

#srad
dist.sr <- vegdist(seq.sr)
samdf.sr <- subset(samdf.all, host_species == "radians")
row.names(samdf.sr)==row.names(seq.sr)
bet.sr <- betadisper(dist.sr,samdf.sr$m_y_s_z)
anova(bet.sr) #ps.lesr ns
plot(bet.sr)

permutest(bet.sr, pairwise = TRUE, permutations = 999)

adonis2(dist.sr ~ area*m_y,data=samdf.sr, permutations=999) #sig!

#go through each iteration of pairwise comparisons bc intxns don't work
pairwise.adonis2(dist.sr ~ m_y_s_z, data=samdf.sr, permutations=999)

#ppor
dist.pp <- vegdist(seq.pp)
samdf.pp <- subset(samdf.all, host_species == "porites")
row.names(samdf.pp)==row.names(seq.pp)

bet.pp <- betadisper(dist.pp,samdf.pp$m_y_s_z)
anova(bet.pp) #ns
plot(bet.pp) #overlap ns

#permutest(bet.pp, pairwise = TRUE, permutations = 999) #ns

adonis2(dist.pp ~ reef_bay*m_y, data=samdf.pp, permutations=999) #ns

pairwise.adonis2(dist.pp ~ m_y, data=samdf.pp, permutations=999)

####Bar plots####
#data - relative abundance dataframes
#ps.all.rel
#ps.ss.rel
#ps.sr.rel
#ps.pp.rel

ps.all.phy <- aggregate_rare(ps.all.rel,
                             level = "Phylum",
                             detection = 0.1,   # 10% RA threshold
                             prevalence = 0,      # keep all genera, just lump low-RA
                             aggregate_name = "Other")
#taxa.phy <- as.data.frame(ps.all.phy@tax_table)
#unique(taxa.phy$Phylum)
all.phy <- psmelt(ps.all.phy) %>%
  mutate(Phylum = as.character(Phylum),
         Phylum = fct_relevel(factor(Phylum), "Other", after = Inf))

phyla_colors <- c(
  "Actinobacteriota"  = "#E69F00",  # orange
  "Bacteroidota"      = "#88CCEE",  # sky blue
  "Chloroflexi"       = "#339911",  # light green
  "Cyanobacteria"     = "#F0E442",  # yellow
  "Desulfobacterota"  = "#0072B2",  # blue
  "Firmicutes"        = "#D55E00",  # vermillion
  "Fusobacteriota"    = "#CC79A7",  # purple/pink
  "Myxococcota"       = "#999933",  # olive
  "Planctomycetota"   = "#332288",  # dark purple
  "Proteobacteria"    = "#115511",  # dark green
  "Verrucomicrobiota" = "#553223",  # brown
  "Other"             = "darkgray"
)

bac.spp.phy <- ggplot(all.phy, aes(x = full_sample_id, y = Abundance, fill = Phylum)) +
  geom_col() +
  facet_wrap(~ coral_species, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phyla_colors, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 15),axis.ticks.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  xlab("Sample ID") +
  ylab("Relative Abundance")
bac.spp.phy

#now at genus level
ps.all.gen <- aggregate_rare(ps.all.rel,
                              level = "Genus",
                              detection = 0.1,   # 10% RA threshold
                              prevalence = 0,      # keep all genera, just lump low-RA
                              aggregate_name = "Other")

all.gen <- psmelt(ps.all.gen) %>%
  mutate(Genus = as.character(Genus),
         Genus = fct_relevel(factor(Genus), "Other", after = Inf))

ss.gen <- filter(all.gen,host_species=="siderea")
sr.gen <- filter(all.gen,host_species=="radians")
pp.gen <- filter(all.gen,host_species=="porites")

#get colors for top 10% of genera
#taxa.gen <- as.data.frame(ps.all.gen@tax_table)
#unique(taxa.gen$Genus)
genus_cols <- c(
  "Algicola"= "#4477AA",
  "Allofrancisella"= "#EE6677",
  "Alteromonas"= "#228833",
  "Amycolatopsis"= "#CCBB44",
  "Aureibacter"= "#66CCEE",
  "BD1-7 clade"= "#AA3377",
  "Candidatus Amoebophilus"= "#552222",
  "Candidatus Fritschea"= "#332288",
  "Cetobacterium"= "#88CCAA",
  "Class_Alphaproteobacteria"= "#CC6677",
  "Class_Bacteroidia"= "#117733",
  "Cutibacterium"= "#DDCC77",
  "Desulfobacter"= "#77AADD",
  "Desulfocella"= "#44AA99",
  "Dyella"= "#999933",
  "Endozoicomonas"= "#AA4499",
  "Epibacterium"= "#6699CC",
  "Family_Cyclobacteriaceae"= "#66AA55",
  "Family_Desulfobacteraceae"= "#BB5566",
  "Family_Flavobacteriaceae"= "#DDAA33",
  "Family_Rhodobacteraceae"= "#33BBEE",
  "Family_Rickettsiaceae"= "#225555",
  "Fodinicurvata"= "#885588",
  "Fulvivirga"= "#77CCDD",
  "Hormoscilla SI04-45"= "#5566AA",
  "Marinomonas"= "#44BB99",
  "Order_Burkholderiales"= "#AA8833",
  "Other"= "darkgray",
  "Photobacterium"= "#CC4444",
  "Phycisphaera"= "#66C2A5",
  "Phyllobacterium"= "#8DA0CB",
  "Phylum_Proteobacteria"= "#1B9E77",
  "Propionigenium"= "#E6AB02",
  "Prosthecochloris"= "#A6761D",
  "Pseudoalteromonas"= "#7570B3",
  "Psychrosphaera" = "#4DAF4A",
  "Ruegeria"= "#E7298A",
  "Shewanella"= "#377EB8",
  "Shimia"= "#984EA3",
  "Streptococcus"= "#D95F02",
  "Thalassotalea"= "#A6CEE3",
  "Vibrio"= "#234565")

#bars grouped by species
bac.spp.gen <- ggplot(all.gen, aes(x = full_sample_id, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ coral_species, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text = element_text(face = "italic", size = 30)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")
bac.spp.gen

#save species comparison graphs - to add with ordination and diversity plots
#read in alpha diversity plots
all_spp_alpha_plot <- readRDS(file = here("16S","all_spp_alpha_plot.rds")) #labeled with B,C,D,E

alpha.beta <- ggarrange(all_spp_alpha_plot,spp.ord,nrow=1,ncol=2,common.legend = F, 
                        labels = c("A","E"), font.label = list(size = 50),
                        widths = c(3,1))
its2.bar <- readRDS(file = here("ITS2","its2.coralspecies.rds"))
bac.spp <- ggarrange(bac.spp.phy,bac.spp.gen,its2.bar,nrow=3,ncol=1,common.legend = F, 
                     labels = c("F", "G", "H"), font.label = list(size = 50))
all.bac.spp <- ggarrange(alpha.beta,bac.spp,nrow=2,ncol=1,common.legend = F,
                        heights = c(1,6))
ggsave(all.bac.spp, file=here("16S","bac.coralspecies.pdf"),width=30,height=40)


#siderastrea siderea site
bac.ss <- ggplot(ss.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste(italic("Siderastrea siderea"))))
bac.ss

bac.ss.time <- ggplot(ss.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(m_y~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste(italic("Siderastrea siderea"))))
bac.ss.time

#siderastrea radians site
bac.sr <- ggplot(sr.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste(italic("Siderastrea radians"))))
bac.sr

bac.sr.time <- ggplot(sr.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(m_y~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste(italic("Siderastrea radians"))))
bac.sr.time

#branching porites sp
bac.pp <- ggplot(pp.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
bac.pp

bac.pp.time <- ggplot(pp.gen, aes(x = number_time, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(m_y~ site_nice, scales = "free_x") +
  scale_fill_manual(values = genus_cols, drop = FALSE) +
  theme_classic(base_size = 30) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 20)) +
  xlab("Sample ID") +
  ylab("Relative Abundance")+
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
bac.pp.time

###put plots together site
bac.sr.pp <- ggarrange(bac.sr,bac.pp, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 50))
bac.genus.site <- ggarrange(bac.ss,bac.sr.pp, nrow=1, ncol = 2, legend = "bottom", common.legend = TRUE, labels = c("A"),font.label = list(size = 50))
bac.genus.site
ggsave(bac.genus.site, file=here("16S","bac.genus.site.pdf"),width=25,height=20)

#put plots together site timepoint
bac.time.sr.pp <- ggarrange(bac.sr.time, bac.pp.time, ncol = 2, 
                            labels = c("B", "C"), font.label = list(size = 50),
                            legend = "none")
bac.genus.site.time <- ggarrange(bac.ss.time,bac.time.sr.pp,
                          nrow = 2, heights = c(1.5,1),
                          labels = c("A"),font.label = list(size = 50),
                          legend = "bottom", common.legend = T)
ggsave(bac.genus.site.time, file=here("16S","bac.genus.site.time.pdf"),width=25,height=30)

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
ggsave(gg.bar.ss,file=here("16S","bac.phy.barplot.ss.site.pdf"),h=10,w=35)
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
ggsave(gg.bar.ss,file=here("16S","bac.phy.barplot.ss.site.season.sampleid.pdf"),h=15,w=30)    

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
# ggsave(gg.bar.sr,file=here("16S","bac.phy.barplot.sr.site.pdf"),h=5,w=15)
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
ggsave(gg.bar.sr,file=here("16S","bac.phy.barplot.sr.site.pdf"),h=10,w=30)

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
ggsave(gg.bar.sr,file=here("16S","bac.phy.barplot.sr.site.season.sampleid.pdf"),h=15,w=15)    

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
# ggsave(gg.bar.pp,file=here("16S","bac.phy.barplot.pp.site.pdf"),h=5,w=15)
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
ggsave(gg.bar.pp,file=here("16S","bac.phy.barplot.pp.site.pdf"),h=10,w=27)

twotogether <- ggarrange(gg.bar.sr,gg.bar.pp,nrow = 1, legend = "none", labels = c("B","C"),font.label = list(size = 50))
composition_all_phylum_ps.less <- ggarrange(gg.bar.ss,twotogether,nrow = 1, common.legend = TRUE, legend = "bottom",labels = c("A"), font.label = list(size = 50))
ggsave(composition_all_phylum_ps.less, file=here("16S","composition_all_phylum_ps.less.pdf"),width=30,height=15)

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
ggsave(gg.bar.pp,file=here("16S","bac.phy.barplot.pp.site.season.sampleid.pdf"),h=15,w=15)    

gg.bar.pp <- plot_bar(ps.cl.pp,"number",fill="Phylum")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(m_y~site_zone, scales = "free", ncol = 2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site & Timepoint")+
  ylab("Relative Abundance")
gg.bar.pp

####RDA of Symbiont type####
library(microViz)

#Bray distance & sample data
#use distance matricies and sample data from above
#dist.ss
#dist.sr
#dist.pp
#samdf.ss
#samdf.sr
#samdf.pp

#Does dominant_type explain bacteria after controlling for site?
ss_mod_symb <- capscale(dist.ss ~ dominant_type + Condition(site_nice), data = samdf.ss)
anova.cca(ss_mod_symb, permutations = 999) #ns
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist.ss ~ dominant_type + Condition(site_nice), data = samdf.ss)
# Df SumOfSqs      F Pr(>F)
# Model     2   0.8323 1.0648  0.305
# Residual 68  26.5776   

#Does site explain bacteria after controlling for dominant_type?
ss_mod_site <- capscale(dist.ss ~ site_nice + Condition(dominant_type), data = samdf.ss)
anova.cca(ss_mod_site, permutations = 999) #sig

# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist.ss ~ site_nice + Condition(dominant_type), data = samdf.ss)
# Df SumOfSqs      F Pr(>F)    
# Model     3   1.9987 1.7046  0.001 ***
#   Residual 68  26.5776                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ss.rda.partial <- ordinate(ps.cl.ss, method = "CAP", distance = "bray",
                           formula = ~ dominant_type + Condition(site_nice))
its2.rda <- plot_ordination(ps.cl.ss, ss.rda.partial, color = "dominant_type", shape = "site_nice") +
  geom_point(size = 3) +
  stat_ellipse(type = "t") +
  theme_classic()
its2.rda

#siderastrea radians
#Does dominant_type explain bacteria after controlling for site?
sr_mod_symb <- capscale(dist.sr ~ dominant_type + Condition(site_nice), data = samdf.sr)
anova.cca(sr_mod_symb, permutations = 999) #ns
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist.sr ~ dominant_type + Condition(site_nice), data = samdf.sr)
# Df SumOfSqs      F Pr(>F)
# Model     3   1.2389 1.0833  0.253
# Residual 27  10.2924     

#Does site explain bacteria after controlling for dominant_type?
sr_mod_site <- capscale(dist.sr ~ site_nice + Condition(dominant_type), data = samdf.sr)
anova.cca(sr_mod_site, permutations = 999) #sig
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist.sr ~ site_nice + Condition(dominant_type), data = samdf.sr)
# Df SumOfSqs      F Pr(>F)    
# Model     1   0.8893 2.3329  0.001 ***
#   Residual 27  10.2924 

sr.rda.partial <- ordinate(ps.cl.sr, method = "CAP", distance = "bray",
                           formula = ~ dominant_type + Condition(site_nice))
its2.rda <- plot_ordination(ps.cl.sr, sr.rda.partial, color = "dominant_type", shape = "site_nice") +
  geom_point(size = 3) +
  stat_ellipse(type = "t") +
  theme_classic()
its2.rda

#porites sp
#Does dominant_type explain bacteria after controlling for site?
pp_mod_symb <- capscale(dist.pp ~ dominant_type + Condition(site_nice), data = samdf.pp)
anova.cca(pp_mod_symb, permutations = 999) #not possible

# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist.pp ~ dominant_type + Condition(site_nice), data = samdf.pp)
# Df SumOfSqs      F Pr(>F)
# Model     2   0.9570 1.1367      1
# Residual 11   4.6305   

#Does site explain bacteria after controlling for dominant_type?
pp_mod_site <- capscale(dist.pp ~ site_nice + Condition(dominant_type), data = samdf.pp)
anova.cca(pp_mod_site, permutations = 999) #not possible

# No constrained component
# 
# Model: capscale(formula = dist.pp ~ site_nice + Condition(dominant_type), data = samdf.pp)
# Df SumOfSqs  F Pr(>F)
# Model     0   0.0000  0       
# Residual 14   4.6305  

pp.rda.partial <- ordinate(ps.cl.pp, method = "CAP", distance = "bray",
                           formula = ~ dominant_type + Condition(site_nice))
its2.rda <- plot_ordination(ps.cl.pp, pp.rda.partial, color = "dominant_type", shape = "site_nice") +
  geom_point(size = 3) +
  stat_ellipse(type = "t") +
  theme_classic()
its2.rda

#generating CAP plots of environmental variables
library(ggcorrplot)
library(ggcorrplot2)
#first briefly look at correlations between environmental variables
env_vars <- c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site","mean_sal_site")
              #"range_temp_site","range_pH_site","range_do_site","range_par_site","range_sal_site")

ct_ss <- corr.test(samdf.ss %>% select(site_nice, all_of(env_vars)) %>%
                     transmute(
                       `Mean Temp`  = mean_temp_site,
                       `Mean pH`    = mean_pH_site,
                       `Mean DO`    = mean_do_site,
                       `Mean PAR`   = mean_par_site,
                       `Mean Salinity`   = mean_sal_site,
                       #`Range Temp` = range_temp_site,
                       #`Range pH`   = range_pH_site,
                       #`Range DO`   = range_do_site,
                       #`Range PAR`  = range_par_site,
                       #`Range Salinity`  = range_sal_site
                     ), adjust = "none")
#extract
corr_ss  <- ct_ss$r
p_ss     <- ct_ss$p
#plot
p_corr_ss <- ggcorrplot.mixed(corr_ss, 
                              upper = "ellipse", 
                              lower = "number", 
                              p.mat = p_ss, 
                              insig = "label_sig", 
                              sig.lvl = c(0.05, 0.01, 0.001))

p_corr_ss

#now test using dbRDA

#siderastrea siderea
null_ss <- capscale(dist.ss ~ 1, data = samdf.ss)
full_ss <- capscale(dist.ss ~ mean_temp_site + mean_pH_site + mean_do_site + mean_par_site + mean_sal_site + range_temp_site + range_pH_site + range_do_site + range_par_site + range_sal_site, data = samdf.ss)
#Some constraints or conditions were aliased because they were redundant. This can happen if terms are constant or linearly
#dependent (collinear): ‘mean_par_site’, ‘mean_sal_site’, ‘range_temp_site’, ‘range_pH_site’, ‘range_do_site’,
#‘range_par_site’, ‘range_sal_site’
alias(full_ss)
#keep mean temp, pH, and do

mod_ss <- capscale(dist.ss ~ mean_temp_site + mean_pH_site + mean_do_site, data = samdf.ss)

#permutation tests for the selected model
anova.cca(mod_ss, permutations = 999)
# Model: capscale(formula = dist.ss ~ mean_temp_site + mean_pH_site + mean_do_site, data = samdf.ss)
#          Df SumOfSqs      F Pr(>F)    
# Model     3   2.3237 1.9781  0.001 ***
anova.cca(mod_ss, by = "term", permutations = 999)
# Model: capscale(formula = dist.ss ~ mean_temp_site + mean_pH_site + mean_do_site, data = samdf.ss)
#                Df SumOfSqs      F Pr(>F)    
# mean_temp_site  1   0.6119 1.5628  0.010 ** 
# mean_pH_site    1   0.6818 1.7412  0.006 ** 
# mean_do_site    1   1.0299 2.6302  0.001 ***

#siderastrea radians
null_sr <- capscale(dist.sr ~ 1, data = samdf.sr)
full_sr <- capscale(dist.sr ~ mean_temp_site + mean_pH_site + mean_do_site + mean_par_site + mean_sal_site + range_temp_site + range_pH_site + range_do_site + range_par_site + range_sal_site, data = samdf.sr)
#Some constraints or conditions were aliased because they were redundant. This can happen if terms are constant or linearly
#dependent (collinear): ‘mean_pH_site’, ‘mean_do_site’, ‘mean_par_site’, ‘mean_sal_site’, ‘range_temp_site’, ‘range_pH_site’,
#range_do_site’, ‘range_par_site’, ‘range_sal_site’
alias(full_sr)
#keep only mean temp
#ALL CO-LINEAR BC JUST 2 SITES

mod_sr <- capscale(dist.sr ~ mean_temp_site, data = samdf.sr)

#permutation tests for the selected model
anova.cca(mod_sr, permutations = 999)
# Model: capscale(formula = dist.sr ~ mean_temp_site, data = samdf.sr)
#          Df SumOfSqs      F Pr(>F)    
# Model     1   0.8802 2.2901  0.001 ***
# Residual 30  11.5312  

#branching porites sp
null_pp <- capscale(dist.pp ~ 1, data = samdf.pp)
full_pp <- capscale(dist.pp ~ mean_temp_site + mean_pH_site + mean_do_site + mean_par_site + mean_sal_site + range_temp_site + range_pH_site + range_do_site + range_par_site + range_sal_site, data = samdf.pp)
#Some constraints or conditions were aliased because they were redundant. This can happen if terms are constant or linearly
#dependent (collinear): ‘mean_pH_site’, ‘mean_do_site’, ‘mean_par_site’, ‘mean_sal_site’, ‘range_temp_site’, ‘range_pH_site’,
#range_do_site’, ‘range_par_site’, ‘range_sal_site’
alias(full_pp)
#keep only mean temp
#ALL CO-LINEAR BC JUST 2 SITES

mod_pp <- capscale(dist.pp ~ mean_temp_site, data = samdf.pp)

#permutation tests for the selected model
anova.cca(mod_pp, permutations = 999)
# Model: capscale(formula = dist.pp ~ mean_temp_site, data = samdf.pp)
#          Df SumOfSqs      F Pr(>F)  
# Model     1   0.5541 1.2891  0.056 .
# Residual 13   5.5875        
#and not even different what the heck

#siderastrea siderea
bray.ss <- ps.cl.ss %>% tax_transform("identity") %>% dist_calc("bray")
ss.perm <- bray.ss %>% dist_permanova(variables = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site"), seed = 321)
cap.site.ss <- bray.ss %>%
  ord_calc(constraints = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_nice", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(type = "text", size = 8, aspect_ratio = 1, colour = "black", justify = "side")) +
  stat_ellipse(aes(colour = site_nice), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise","lightsalmon1","deepskyblue4"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme_classic(base_size = 22)
cap.site.ss
#ggsave(cap.site.ss, file=here("16S","CAP.env.ss.ps.less.pdf"))


#siderastrea radians
bray.sr <- ps.cl.sr %>% tax_transform("identity") %>% dist_calc("bray")
sr.perm <- bray.sr %>% dist_permanova(variables = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site"), seed = 321)
cap.site.sr <- bray.sr %>%
  ord_calc(constraints = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_nice", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(type = "text", size = 8, aspect_ratio = 1, colour = "black", justify = "side")) +
  stat_ellipse(aes(colour = site_nice), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(name = "Site", values = c("coral3","lightsalmon1"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme_classic(base_size = 22)
cap.site.sr

#porites sp
bray.pp <- ps.cl.pp %>% tax_transform("identity") %>% dist_calc("bray")
pp.perm <- bray.pp %>% dist_permanova(variables = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site"), seed = 321)
cap.site.pp <- bray.pp %>%
  ord_calc(constraints = c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site")) %>%
  ord_plot(
    colour = "site_nice", size = 2.5,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(1.5, colour = "black"),
    constraint_lab_style = constraint_lab_style(type = "text", size = 8, aspect_ratio = 1, colour = "black", justify = "side")) +
  stat_ellipse(aes(colour = site_nice), linewidth = 0.2) + # linewidth not size since ggplot 3.4.0
  scale_color_manual(name = "Site", values = c("coral3", "paleturquoise"))+
  coord_fixed(ratio = 1, clip = "off") +
  theme_classic(base_size = 22)
cap.site.pp

#stats
env_vars <- c("mean_temp_site","mean_pH_site","mean_do_site","mean_par_site", "mean_sal_site")

cor(sam_sr2[, env_vars], use = "pairwise.complete.obs")

ss_mod_symb <- capscale(dist.ss ~ dominant_type + Condition(site_nice), data = samdf.ss)
anova.cca(ss_mod_symb, permutations = 999) #ns

# dbRDA: constrain by dominant_type
ss_dom_rda <- capscale(dist.ss ~ dominant_type, data = samdf.ss)

# Overall test (does dominant_type explain variation?)
anova_overall <- anova.cca(ss_dom_rda, permutations = 999)
anova_terms   <- anova.cca(ss_dom_rda, by = "terms", permutations = 999)
anova_axes    <- anova.cca(ss_dom_rda, by = "axis", permutations = 999)

anova_overall
anova_terms
anova_axes

# Effect size: proportion constrained
RsquareAdj(ss_dom_rda)

ggplot(ss_dom_rda) +
  geom_point(aes(color = dominant_type), size = 3)
  #stat_ellipse(aes(color = dominant_type), linetype = 2) +
  #labs(title = "dbRDA: dominant_type") +
  #theme_minimal()

# ####Microshades plots - not using anymore#####
# #remotes::install_github("KarstensLab/microshades")
# library(microshades)
# #remotes::install_github("mikemc/speedyseq")
# library(speedyseq)
# 
# taxa.less <- as.data.frame(ps.less.nd@tax_table)
# 
# #all
# #porites porites
# ps_phy_glom_all <- tax_glom(ps.less.nd, "Phylum")
# # Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
# mdf_prep_all <- prep_mdf(ps.less.nd, subgroup_level = "Family")
# 
# # Create a color object for the specified data
# color_obj_all <- create_color_dfs(mdf_prep_all, selected_groups = c("Proteobacteria","Bacteroidota","Desulfobacterota","Planctomycetota","Firmicutes"), group_level = "Phylum", subgroup_level = "Family", cvd = TRUE)
# #can only do up to 5 based on the palette number
# #using top 5 based on colors here
# sums <- as.data.frame(taxa_sums(ps_phy_glom_all))
# View(data.frame(ps.less.nd@tax_table))
# #MAKE SURE TO SORT BY MOST ABUNDANT
# #sq2 Proteobacteria
# #sq16 Bacteroidota
# #sq22 Desulfobacterota
# #sq132 Planctomycetota
# #sq35 Firmicutes
# #sq47 Cyanobacteria
# 
# #most abundant taxa for all overall determined by above are:
# 
# # Extract
# mdf_all <- color_obj_all$mdf
# cdf_all <- color_obj_all$cdf
# 
# #plot all - species comparisons
# ms_bar_plot_all <- plot_microshades(mdf_all, cdf_all, group_label = "Phylum - Family")
# ms_bar_plot_final <- ms_bar_plot_all + scale_y_continuous(labels = scales::percent, expand = expansion(0))+
#   theme_classic(base_size=30)+
#   theme(legend.position = "bottom", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   facet_wrap(host_species~.,scales = "free")+
#   #labeller=labeller(site_zone = c("Santa Martha Bay", "Santa Martha Reef", "Spaanse Water Bay", "Spaanse Water Reef")))+
#   xlab("Coral Species")+
#   ylab("Relative Abundance")
#   #ggtitle(substitute(paste(italic("All Species"))))
# ms_bar_plot_final
# ggsave(ms_bar_plot_final, file = "ms.all.spp.ps.less.pdf"), w=30, h=10)
# 
# #subset different groups out
# mdf_ss <- mdf_all %>% filter(host_species == "siderea")
# mdf_sr <- mdf_all %>% filter(host_species == "radians")
# mdf_pp <- mdf_all %>% filter(host_species == "porites")
# 
# #porites porites
# 
# #plot pp
# ms_bar_plot_pp <- plot_microshades(mdf_pp, cdf_all, group_label = "Phylum - Family")
# ms_bar_plot_pp_final <- ms_bar_plot_pp + scale_y_continuous(labels = scales::percent, expand = expansion(0))+
#   theme_classic(base_size=22)+
#   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   facet_wrap(.~m_y,scales = "free")+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   ggtitle(substitute(paste(italic("Porites sp."))))
# ms_bar_plot_pp_final
# ggsave(ms_bar_plot_pp_final, file = "ms.pp.ps.less.pdf"), w=30, h=10)
# ggsave(ms_bar_plot_pp_final, file = "ms.pp.time.ps.less.pdf"), w=30, h=10)
# 
# #siderastrea siderea
# #plot ss
# ms_bar_plot_ss <- plot_microshades(mdf_ss, cdf_all, group_label = "Phylum - Family")
# ms_bar_plot_ss_final <- ms_bar_plot_ss + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
#   theme_classic(base_size=22)+
#   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   facet_wrap(.~site_nice,scales = "free")+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   ggtitle(substitute(paste(italic("Siderastrea siderea"))))
# ms_bar_plot_ss_final
# ggsave(ms_bar_plot_ss_final, file = here("16S","ms.ss.ps.less.pdf"), w=30, h=15)
# ggsave(ms_bar_plot_ss_final, file = "ms.ss.time.ps.less.pdf"), w=30, h=10)
# 
# #siderastrea radians
# #plot sr
# ms_bar_plot_sr <- plot_microshades(mdf_sr, cdf_all, group_label = "Phylum - Family")
# ms_bar_plot_sr_final <- ms_bar_plot_sr + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
#   theme_classic(base_size=22)+
#   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   facet_wrap(.~m_y,scales = "free")+
#   xlab("Site")+
#   ylab("Relative Abundance")+
#   ggtitle(substitute(paste(italic("Siderastrea radians"))))
# ms_bar_plot_sr_final
# ggsave(ms_bar_plot_sr_final, file = "ms.sr.ps.less.pdf"), w=30, h=10)
# ggsave(ms_bar_plot_sr_final, file = "ms.sr.time.ps.less.pdf"), w=30, h=10)
# 
# #saving plots
# ggsave(ms_bar_plot_sr_final, file = "microshades_bar_sr.ps.less.png",h=5,w=20)
# ggsave(ms_bar_plot_ss_final, file = "microshades_bar_ss.ps.less.png",h=10,w=20)
# ggsave(ms_bar_plot_pp_final, file = "microshades_bar_pp.ps.less.png",h=5,w=20)
# 
# #putting plots together
# bac_pp_sr_plot <- ggarrange(ms_bar_plot_sr_final,ms_bar_plot_pp_final, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 50))
# bac_all_plot <- ggarrange(ms_bar_plot_ss_final,bac_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 50))
# ggsave(bac_all_plot, file=here("16S","microshades_all_bar.ps.less.pdf"),width=30,height=12)
# 
# #putting time plots together
# bac_all_time_plot <- ggarrange(ms_bar_plot_ss_final,ms_bar_plot_sr_final,ms_bar_plot_pp_final, nrow=3, ncol = 1, legend = "none", common.legend = TRUE, labels = c("A","B","C"),font.label = list(size = 50))
# ggsave(bac_all_time_plot, file=here("16S","microshades_all_bar_time.ps.less.pdf"),width=15,height=20)
# 
# #just ss and pp for emes poster april 2nd 2024
# ms_ss_pp <- ggarrange(ms_bar_plot_ss_final,ms_bar_plot_pp_final,nrow = 2, ncol=1, heights= c(2,1), common.legend = TRUE, legend = "right", labels = c("A","B"),font.label = list(size = 50))
# ggsave(ms_ss_pp, file=here("16S","microshades_ss_pp_bar_ps.less.png",width=20,height=15)
