######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S: plot sample depth and remove samples with low sample counts
######################

# In this script we plot the sequencing depth of all samples, and then remove samples with depth below 1000 reads, together with OTUs with counts less than 3 across the dataset

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/prev_16S_physeq_filtered.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# generate sample summary ####
# filter to OTUs present in counts over 3
ps.lcr <- prune_taxa(taxa_sums(ps) > 3, ps) # 4492 OTUs total

# Summary depth #### 
depth <- colSums(otu_table(ps.lcr))
sample_data(ps.lcr)$depth <- depth

# plot depth #### 
dat <- data.frame(sample_data(ps.lcr))

# set factor levels
dat$group.ages.short <- factor(dat$group.ages.short, levels = as.character(cols$group.ages.short))

dat$remove <- ifelse(dat$depth < 1000, "Y", "N")
dat$sample.type <- gsub("breastmilk", "bm", dat$sample.type)

check <- filter(dat, remove == "Y") # 5 samples

ggplot(dat, aes(x = group.ages.short, y = depth/1000)) +
  geom_boxplot(outlier.size = NA) +
  geom_point(position = position_jitter(width = 0.35), aes(color = remove)) +
  scale_color_manual(values = c("black", "blue")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        legend.position = "none") +
  labs(y = "Number of Reads (thousands)") +
  geom_hline(yintercept = 1, color = "blue") +
  facet_grid(sample.type~., space = "free", scales = "free") +
  coord_flip()

#ggsave("ms_bacteria_fungi_analysis/figures/sequencing_depth/prev_16S_depth.pdf", device = "pdf", dpi = 300, width = 5, height = 5.3)

# generate final PS for analysis ####
ps.final <- ps.lcr
ps.final <- subset_samples(ps.final, depth > 1000)

# Samples PB106 and PB190 were recruited outside of their age brackets.  Thus, they will be omitted from further analysis.  
subj.rm <- c("PB106", "PB190")

subjects <- unique(as.character(dat$stdyid))
subjects.keep <- subjects[-which(subjects %in% subj.rm)]

ps.final <- subset_samples(ps.final, stdyid %in% subjects.keep)

# save final physeq ####
#saveRDS(ps.final, "ms_bacteria_fungi_analysis/Rdata/R_export/ps_16S_for_analysis.rds")

