######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Beta Diversity
######################

# Ordinate all samples and use PERMANOVA to determine variance explained by age

# load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)
library(vegan)

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_ITS2_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# Ordination:  NMDS on Bray-Curtis Distance
ord.nmds <- ordinate(ps, method = "NMDS", distance = "bray")

p.ord <- plot_ordination(ps, ord.nmds)
ord.data <- p.ord$data
ord.labs <- list(x = p.ord$labels$x, y = p.ord$labels$y)

ord.data$group.ages.short <- factor(ord.data$group.ages.short, levels = as.character(cols$group.ages.short))

ggplot(ord.data, aes(x = NMDS1, y = NMDS2, fill = group.ages.short, shape = sample.type)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_fill_manual(values = as.character(cols$col)) +
  scale_shape_manual(values = c(21, 24)) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.position = "right",
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  labs(x = ord.labs$x, y = ord.labs$y)

#ggsave("ms_bacteria_fungi_analysis/figures/beta_diversity/ord_ITS2_nmds.pdf", device = "pdf", dpi = 300, width = 6.8, height = 4)

# adonis for age ####
# subset data for only child stool samples
ps.stool <- subset_samples(ps, sample.type == "stool")
ps.stool <- subset_samples(ps, moth.inf == "Child")

stool.meta <- data.frame(sample_data(ps.stool))

otu.tab <- data.frame(t(otu_table(ps.stool)))

ad.test <- adonis(otu.tab ~ group.ages, data = stool.meta, method = "bray")
ad.test # R2 = 0.120, p = 0.003.

ad.test <- adonis(otu.tab ~ age.in.days, data = stool.meta, method = "bray")
ad.test # R2 = 0.017, p = 0.044

# END ####

# SessionInfo ####
sessionInfo()

#R version 3.3.2 (2016-10-31)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: macOS  10.15.4

#locale:
#  [1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] vegan_2.4-5     lattice_0.20-38 permute_0.9-5   ggplot2_3.1.0   phyloseq_1.19.1
#[6] dplyr_0.8.0.1   plyr_1.8.4