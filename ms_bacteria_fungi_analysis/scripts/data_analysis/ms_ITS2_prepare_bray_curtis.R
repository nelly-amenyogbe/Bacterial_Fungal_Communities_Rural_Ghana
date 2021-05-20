######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Beta Diversity: Bray-Curtis distance between all samples
######################

# In this script, we create a Bray-Curtis distance matrix, and transform this into a data frame of every pairwise combination.  Metadata is added for the first and second pair. This data will be exported for the next analysis.

## Load packages
library(phyloseq)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(vegan)

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_ITS2_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# 1.  Bray-Curtis distnances between all samples ####
# get BC data
otu.table <- data.frame(t(otu_table(ps)))
dist <- vegdist(otu.table, method="bray")

# create pairwise data frame
mat <- as.matrix(dist) # 131*131 = 17161

mat <- as.data.frame(mat)
mat$id <- rownames(mat)
m.melt <- melt(mat, id.vars = "id")

rm <- which(m.melt$id == m.melt$variable) # remove sample-sample pairs. 131 flagged
m.melt <- m.melt[-c(rm),] # 17030 (17161-131)

# add metadata
meta <- data.frame(sample_data(ps))

colnames(m.melt)[1] <- "MBI.sample.ID"
m.melt <- join(m.melt, meta, by = "MBI.sample.ID")

# add metadata for pair
colnames(meta)
meta.pair <- meta[,c("MBI.sample.ID", "sample.name", "stdyid", "sample.type", "partna", "age.in.days", "pair.id", "moth.inf", "group.ages", "group.ages.short")]

colnames(meta.pair) <- c("MBI.sample.ID", "pair.sample.name", "pair.subject", "pair.sample.type", "pair.part", "pair.age.d", "pair.pid", "pair.mf", "pair.ga", "pair.gas")

colnames(meta.pair)[1] <- "variable"

m.melt <- join(m.melt, meta.pair, by = "variable") # this dataset has all pairwise combinations, with metadata for the first and second pair.

check <- filter(m.melt, MBI.sample.ID == "S0014-0001")
check2 <- filter(check, pair.sample.type == "breastmilk") # one child stool is paired with 23 breastmilks.  

# save data 
#write.csv(m.melt, "ms_bacteria_fungi_analysis/Rdata/R_export/ITS2_bray_curtis_dist.csv")
