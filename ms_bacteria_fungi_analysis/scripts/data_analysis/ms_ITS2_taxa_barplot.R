######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Taxonomic composition of samples
######################

# In this script we will:
# 1. Transform count to relative abundace data
# 2. Generate barplots for top-10 community composition for breastmilk and stool samples separately

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
# helper functions
source("ms_bacteria_fungi_analysis/scripts/functions/function_plot_taxa.R")

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_ITS2_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# prepare group.ages data
sample_data(ps)$group.ages.short <- factor(sample_data(ps)$group.ages.short, levels = as.character(cols$group.ages.short))

# 1. Transform to relative abundance data ####
ps.rab <- transform_sample_counts(ps, function(x) x/sum(x))

# check your new OTU table
otu.tab <- data.frame(otu_table(ps.rab))

# each sample should add up to 1
check <- colSums(otu.tab)
unique(check) #:) all is good

# 2. Prepare data for plotting ####
# subset sample types
stool <- subset_samples(ps.rab, sample.type == "stool")
bm <- subset_samples(ps.rab, sample.type == "breastmilk")

# Genus-level data ####

# stool: individual samples
genus.s <- get_taxbar_data(stool, num.taxa = 15, rank = "Genus")

# breastmilk: by individual samples
genus.b <- get_taxbar_data(bm, num.taxa = 20, rank = "Genus")

# stool, by age stratification
genus.s.g <- get_taxbar_data.grp(stool, num.taxa = 25, rank = "Genus", grp = "group.ages")

### Plot ####
cols.10 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a', "#525252",  "#969696")

cols.15 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', "#ae017e", "#cc4c02", "#fe9929", "#525252","#969696")


# plot taxa ####

# Function for individual barplots
plot.ind <- function(data, cols, title){
  
  ggplot(data, aes(x = MBI.sample.ID, y = Abundance, fill = Taxa)) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    facet_grid(~group.ages.short, scale = "free", space = "free") + 
    scale_fill_manual(title, values = c(cols, "grey")) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.ticks.x = element_blank()) +
    #ggtitle(title) +
    labs(x = "", y = "Relative Abundance")
  
}

# Function for barplots by metadata group
plot.groups <- function(data, cols, title){
  
  ggplot(data, aes(x = group.ages.short, y = rel.abund, fill = Taxa)) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    #facet_grid(~group.ages, scale = "free", space = "free") + 
    scale_fill_manual(title, values = c(cols, "grey")) +
    theme(axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          strip.text.x = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.ticks.x = element_blank()) +
    #ggtitle(title) +
    labs(x = "", y = "Relative Abundance") +
    coord_flip()
  
}

# Genus-stool
# add annotation
group.annot <- cols[,c("group.ages", "group.ages.short")]
genus.s.g <- join(genus.s.g, group.annot, by = "group.ages")
genus.s.g$group.ages.short <- factor(genus.s.g$group.ages.short, levels = as.character(cols$group.ages.short))

p.gen.sg <- plot.groups(genus.s.g,
                        c(cols.15[1:15], cols.10, cols.15[16:17]),
                        "Genus")

p.gen.sg

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/stool_ITS2_taxabars.pdf", device = "pdf", dpi = 300, width = 10, height = 4.5)

# Breastmilk ind:  Genus ####
p.bm <- plot.ind(genus.b,
                 c(cols.15[1:15], cols.10[1:5], cols.15[16:17]),
                 "Genus")
p.bm

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/bm_ITS2_taxabars.pdf", device = "pdf", dpi = 300, width = 10, height = 4.5)
