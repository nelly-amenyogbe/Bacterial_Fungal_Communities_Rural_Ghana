######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S: Beta Diversity: Bray-Curtis distance between all samples
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
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_16S_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# 1.  Bray-Curtis distnances between all samples ####
# get BC data
otu.table <- data.frame(t(otu_table(ps)))
dist <- vegdist(otu.table, method="bray")

# create pairwise data frame
mat <- as.matrix(dist) # 187*187 = 34969

mat <- as.data.frame(mat)
mat$id <- rownames(mat)
m.melt <- melt(mat, id.vars = "id")

rm <- which(m.melt$id == m.melt$variable) # remove sample-sample pairs. 187 flagged
m.melt <- m.melt[-c(rm),] # 34782 (34969-34782)

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
check2 <- filter(check, pair.sample.type == "breastmilk") # one child stool is paired with 27 breastmilks.  

# plot similarity to newborn stool

# to newborn stool

nb.dat <- filter(m.melt, pair.ga == "DOL 0-5", sample.type == "stool", pair.sample.type == "stool")

# set colours
p.cols <- as.character(cols$col)
names(p.cols) <- cols$group.ages.short

ggplot(filter(nb.dat, group.ages %in% c("Mothers_DOL 0-5", "Mothers_DOL 26-35")), aes(x = group.ages.short, y = 1-value, fill = group.ages.short)) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = p.cols) +
  labs(y = "Similarity", x = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14)) +
  ggtitle("Similarity to newborn stool")

#ggsave("ms_bacteria_fungi_analysis/figures/beta_diversity/prev_16s_nb_st_distance.pdf", device = "pdf", dpi = 300, width = 5, height = 2.5)

# test for significance of mothers to newborns

mot.df <- filter(nb.dat, group.ages %in% c("Mothers_DOL 0-5", "Mothers_DOL 26-35"))

moms.05 <- filter(mot.df, group.ages == "Mothers_DOL 0-5")
moms.26 <- filter(mot.df, group.ages == "Mothers_DOL 26-35")

wilcox.test(moms.05$value, moms.26$value) # < 2.2e-16

# save data 
#write.csv(m.melt, "ms_bacteria_fungi_analysis/Rdata/R_export/16S_bray_curtis_dist.csv")
