######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S and ITS3: Beta Diversity: Bray-Curtis distance between all samples: cross-kingdom comparison
######################

# In this script, we compare the bray-curtis distance between all age stratifications with the oldest children for bacteria and fungi

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
# helper function
source("ms_bacteria_fungi_analysis/scripts/functions/function_sig_annot.R")

# load data
s6 <- read.csv("ms_bacteria_fungi_analysis/Rdata/R_export/16S_bray_curtis_dist.csv")
it <- read.csv("ms_bacteria_fungi_analysis/Rdata/R_export/ITS2_bray_curtis_dist.csv")
cols <- read.csv("ms_bacteria_fungi_analysis/Rdata/group_colours.csv")

# combine data
s6$community <- "Bacteria"
it$community <- "Fungi"

# get matching column names
s6 <- s6[,c("MBI.sample.ID", "variable", "value", "Sample_Plate", "Sample_Well", "sample.name", "sample.type", "stdyid", "partna", "pair.id", "age.in.days", "moth.inf", "group.ages", "group.ages.short", "participant.sex", "col", "depth", "pair.sample.name", "pair.subject", "pair.sample.type", "pair.part", "pair.age.d", "pair.pid", "pair.mf", "pair.ga", "pair.gas", "community")]

it <- it[,c("MBI.sample.ID", "variable", "value", "Sample_Plate", "Sample_Well", "sample.name", "sample.type", "stdyid", "partna", "pair.id", "age.in.days", "moth.inf", "group.ages", "group.ages.short", "participant.sex", "col", "depth", "pair.sample.name", "pair.subject", "pair.sample.type", "pair.part", "pair.age.d", "pair.pid", "pair.mf", "pair.ga", "pair.gas", "community")]

dat <- rbind(s6, it)

# prepare data to plot ####
dat$group.ages.short <- factor(dat$group.ages.short, levels = c(as.character(cols$group.ages.short)))

# run child stats ####
#get df
st.dat <- filter(dat, pair.ga == "MOL 57-63", sample.type != "breastmilk", pair.ga != "breastmilk", moth.inf == "Child", pair.mf == "Child") # data of only child-child pairs where the first pair is the oldest age stratification

group.ages.short <- cols$group.ages.short
group.ages.test <- group.ages.short[-which(group.ages.short %in% c("5 years", "Mothers_0-5 days", "Mothers_26-35 days"))]

# bacteria res
st.res.b <- ldply(as.list(group.ages.test), function(i){
  
  dat.test <- filter(st.dat, group.ages.short == i, community == "Bacteria")
  dat.mot <- filter(st.dat, group.ages.short == "5 years", community == "Bacteria")
  
  test <- wilcox.test(dat.test$value, dat.mot$value)
  p.value <- test$p.value
  group.ages.short <- i
  res <- data.frame(group.ages.short, p.value)
  res
  
  
})

st.res.b$p.adj <- p.adjust(st.res.b$p.value, method = "bonferroni")
st.res.b$community <- "Bacteria"

# fungi res
st.res.f <- ldply(as.list(group.ages.test), function(i){
  
  dat.test <- filter(st.dat, group.ages.short == i, community == "Fungi")
  dat.mot <- filter(st.dat, group.ages.short == "5 years", community == "Fungi")
  
  test <- wilcox.test(dat.test$value, dat.mot$value)
  p.value <- test$p.value
  group.ages.short <- i
  res <- data.frame(group.ages.short, p.value)
  res
  
  
})

st.res.f$p.adj <- p.adjust(st.res.f$p.value, method = "bonferroni")
st.res.f$community <- "Fungi"

st.res <- rbind(st.res.b, st.res.f)

# annotate pv
st.res <- add.sig.annot(st.res)

# get min values
st.sum <- ddply(st.dat, "community", summarise, max = max(value))
st.res <- join(st.res, st.sum, by = "community")
st.res <- filter(st.res, group.ages.short != "5 years")

# plot child ####
p.cols <- as.character(cols$col)
names(p.cols) <- cols$group.ages.short

ggplot(filter(st.dat, group.ages!= "MOL 57-63"), aes(x = group.ages.short, y = 1-value, fill = group.ages.short)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_point() +
  theme_classic() +
  facet_grid(~community) +
  coord_flip() +
  scale_fill_manual(values = p.cols) +
  labs(y = "Similarity to 5-year-old children", x = "") +
  theme(axis.text.x = element_text(size = 12, angle = 40, hjust = 1.0, vjust = 1.0),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text =  element_text(size = 14),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 14)) +
  geom_text(data = st.res, aes(x = group.ages.short, y = 0.95, label = sig.annot), inherit.aes = FALSE) +
  ylim(c(0,1))

#ggsave("ms_bacteria_fungi_analysis/figures/beta_diversity/its_s16_child_dist.pdf", device = "pdf", dpi = 300, width = 5.8, height = 4)
