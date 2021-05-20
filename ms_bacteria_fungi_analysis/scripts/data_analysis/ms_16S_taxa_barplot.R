######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S: Taxonomic composition of samples
######################

# In this script we will:
# 1. Transform count to relative abundace data
# 2. Generate barplots for top-10 community composition for breastmilk and stool samples separately
# 3. Generate boxplots to visualize temporal trends for associated taxa over age

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
# helper functions
source("ms_bacteria_fungi_analysis/scripts/functions/function_plot_taxa.R")

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_16S_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# prepare group.ages data
sample_data(ps)$group.ages <- factor(sample_data(ps)$group.ages, levels = as.character(cols$group.ages))

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

# Get taxonomy data ####
# By individual sample
# stool
genus.s <- get_taxbar_data(stool, num.taxa = 15, rank = "Genus")

# breastmilk
genus.bm <- get_taxbar_data(bm, num.taxa = 20, rank = "Genus")

# By group annotation
# stool
genus.s.g <- get_taxbar_data.grp(stool, num.taxa = 25, rank = "Genus", grp = "group.ages")

# Taxa Barplots ####

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
          strip.text.x = element_text(size = 10, face = "bold"),
          strip.background = element_blank(),
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

# annotate by group.ages.short
group.annot <- cols[,c("group.ages", "group.ages.short")]
genus.s.g <- join(genus.s.g, group.annot, by = "group.ages")
genus.s.g$group.ages.short <- factor(genus.s.g$group.ages.short, levels = as.character(group.annot$group.ages.short))

# plot
p.gen.sg <- plot.groups(genus.s.g,
                        c(cols.15[1:15], cols.10, cols.15[16:17]),
                        "Genus")
p.gen.sg

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/stool_16S_taxabars.pdf", device = "pdf", dpi = 300, width = 10, height = 4.5)

# plot bacteroides and prevotella only
bac.prev <- filter(genus.s, Taxa %in% c("g__Bacteroides", "g__Prevotella"))
bac.prev$group.ages.short <- factor(bac.prev$group.ages.short, levels = c(as.character(cols$group.ages.short)))

p.bac.prev <- plot.ind(bac.prev, cols.15[c(2,4)], "Genus")
p.bac.prev

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/bar_prev_bac.pdf", device = "pdf", dpi = 300, width = 12, height = 2.5)

# Breastmilk taxa ####

p.bm <- plot.ind(genus.bm,
                 c(cols.15[1:15], cols.10[1:5], cols.15[16:17]),
                 "Genus")
p.bm

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/bm_16S_taxabars.pdf", device = "pdf", dpi = 300, width = 10, height = 4.5)

# Stool: Boxplots of shifting bacterial taxa ####
# Genera:  Boxplots ####

# create column to differentiate mothers and infants
dat <- plot_bar(stool)
dat <- dat$data

# Aggregate the data by taxonomic rank
dat.ag <- aggregate(Abundance ~ stdyid + Genus + moth.inf + group.ages.short, data = dat, sum)

dat.ag$group.ages.short <- factor(dat.ag$group.ages.short, levels = as.character(cols$group.ages.short))

# taxa boxplots ####

grp.cols <- as.character(cols$col)
names(grp.cols) <- cols$group.ages.short

# Function to generate boxplots of relative abundance based on names of genera selected
plot.box <- function(dat, genera){
  
  dat <- filter(dat, Genus %in% genera)
  
  p <- ggplot(dat, aes(x = group.ages.short, y = Abundance, fill = group.ages.short)) +
    geom_boxplot(outlier.size = NA) +
    scale_fill_manual(values = grp.cols) +
    geom_point(shape = 21, position = position_jitter(width = 0.2)) +
    theme_bw() +
    facet_grid(~moth.inf, space = "free", scales = "free") +
    theme(axis.text.y = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0, size = 12, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.line = element_line(size = 0.8),
          legend.text = element_text(size = 12),
          legend.position = "none",
          strip.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold", size = 14),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white")) +
    ggtitle(genera) +
    labs(x = "", y = "Abundance") 
  
  p
  
}

# these genera were selected based on DESeq2 results, as taxa that differ between mothers by post-partum time

genera <- c("g__Prevotella", "g__Escherichia", "g__Faecalibacterium", "g__Blautia") 

# plot
gen.list <- llply(as.list(genera), function(i){
  p <- plot.box(dat.ag, i)
  p
})

names(gen.list) <- genera

# print
gen.list$g__Prevotella
#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/prev_select.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

gen.list$g__Escherichia
#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/ecoli_select.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

gen.list$g__Faecalibacterium
#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/faec_select.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

gen.list$g__Blautia
#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/blaut_select.pdf", device = "pdf", dpi = 300, width = 5, height = 3.5)

# run stats on select ####
# wilcoxon test for difference in relative abundance between mothers by post-partum time.  Results will be used to annontate the "Mother" panels of the figures above.

wilcox.res <- ldply(as.list(genera), function(i){
  
  dat.05 <- filter(dat.ag, group.ages.short == "Mothers_0-5 days", Genus == i)
  dat.25 <- filter(dat.ag, group.ages.short == "Mothers_26-35 days", Genus == i)
  test <- wilcox.test(dat.05$Abundance, dat.25$Abundance)
  p.value <- test$p.value
  genus <- i
  res <- data.frame(genus, p.value)
  res
  
})

wilcox.res





