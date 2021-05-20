######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Alpha Diversity
######################

# Determine Shannon Diversity and Observed Richness.  Perform wilcoxon test to compare each age stratification to maternal richness

# load packages
library(phyloseq)
library(ggplot2)
library(plyr)
library(dplyr)
library(vegan)

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_ITS2_for_analysis.rds")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# rarefy data ####
depth <- sample_sums(ps)
range(depth) # rarefy to 1212 reads per sample

ps.rar <- rarefy_even_depth(ps, rngseed = 711)

# calculate diversity

p.rich <- plot_richness(ps.rar, measures = c("Observed", "Shannon")) 

# get richness data
rich.data <- p.rich$data
rich.data$variable <- gsub("Observed", "Observed Richness", rich.data$variable)
rich.data$variable <- gsub("Shannon", "Shannon Diversity", rich.data$variable)

# richness statistics ####
#wilcoxon test compared to Mothers_DOL 26-35

groups <- as.character(cols$group.ages.short)
groups <- groups[-which(groups == "Mothers_26-35 days")] # take away the comparator group

# perform wilcox test for all age stratifications
wilcox.res.obs <- ldply(as.list(groups), function(i){
  
  dat.i <- filter(rich.data, variable == "Observed Richness", group.ages.short == i, sample.type == "stool")
  dat.mot <- filter(rich.data, variable == "Observed Richness", group.ages.short == "Mothers_26-35 days", sample.type == "stool")
  test <- wilcox.test(dat.i$value, dat.mot$value)
  p.value <- test$p.value
  group <- i
  variable <- "Observed Richness"
  
  res <- data.frame(group, variable, p.value)
  res
  
})

wilcox.res.obs$p.adj <- p.adjust(wilcox.res.obs$p.value, method = "bonferroni")

# shannon diversity
wilcox.res.sh <- ldply(as.list(groups), function(i){
  
  dat.i <- filter(rich.data, variable == "Shannon Diversity", group.ages.short == i, sample.type == "stool")
  dat.mot <- filter(rich.data, variable == "Shannon Diversity", group.ages.short == "Mothers_26-35 days", sample.type == "stool")
  test <- wilcox.test(dat.i$value, dat.mot$value)
  p.value <- test$p.value
  group <- i
  variable <- "Shannon Diversity"
  
  res <- data.frame(group, variable, p.value)
  res
  
})

wilcox.res.sh$p.adj <- p.adjust(wilcox.res.sh$p.value, method = "bonferroni")

wilcox.res <- rbind(wilcox.res.obs, wilcox.res.sh)

# annotate
wilcox.res$sig.annot <- ifelse(wilcox.res$p.adj <= 0.0001, "****",
                               ifelse(wilcox.res$p.adj <= 0.001, "***",
                                      ifelse(wilcox.res$p.adj <= 0.01, "**",
                                             ifelse(wilcox.res$p.adj <= 0.05, "*",
                                                    ifelse(wilcox.res$p.adj > 0.05, "ns", wilcox.res$p.adj)))))

# Plot ####
# get signif aesthetics
rich.sum <- ddply(rich.data, c("variable"), summarise, max = max(value))

colnames(wilcox.res)[1] <- "group.ages.short"

wilcox.res <- join(wilcox.res, rich.sum, by = "variable")

# set factor levels
rich.data$group.ages.short <- factor(rich.data$group.ages.short, levels = as.character(cols$group.ages.short))

p.rich.stool <- ggplot(filter(rich.data, sample.type == "stool"), aes(x = group.ages.short, y = value, fill = group.ages.short)) +
  geom_boxplot(outlier.size = NA) +
  scale_fill_manual(values = as.character(cols$col)) +
  geom_point(shape = 21, position = position_jitter(width = 0.3)) +
  theme_classic() +
  facet_wrap(~variable, scales = "free") +
  theme(axis.text.x = element_blank(),
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
  ggtitle("Fungi") +
  labs(x = "", y = "Alpha Diversity") + 
  geom_text(data = wilcox.res, aes(x = group.ages.short, y = max+0.2, label = sig.annot), inherit.aes = FALSE)

p.rich.stool

#ggsave("ms_bacteria_fungi_analysis/figures/alpha_diversity/stool_ITS2_alphadiv.pdf", device = "pdf", dpi = 300, width = 9.4, height = 3)

# richness statistics: bm ####

#wilcoxon test compared to Mothers_DOL 26-35: observed

dat.05 <- filter(rich.data, variable == "Observed Richness", group.ages == "Mothers_DOL 0-5", sample.type == "breastmilk")

dat.25 <- filter(rich.data, variable == "Observed Richness", group.ages == "Mothers_DOL 26-35", sample.type == "breastmilk")

test <- wilcox.test(dat.05$value, dat.25$value)
p.value <- test$p.value
p.value # 0.18
variable <- "Observed Richness"
res <- data.frame(variable, p.value)
res

#wilcoxon test compared to Mothers_DOL 26-35: Shannon

dat.05 <- filter(rich.data, variable == "Shannon Diversity", group.ages == "Mothers_DOL 0-5", sample.type == "breastmilk")

dat.25 <- filter(rich.data, variable == "Shannon Diversity", group.ages == "Mothers_DOL 26-35", sample.type == "breastmilk")

test <- wilcox.test(dat.05$value, dat.25$value)
p.value <- test$p.value
p.value # 0.40
variable <- "Shannon Diversity"
res <- data.frame(variable, p.value)
res

# No significant findings

# END ####