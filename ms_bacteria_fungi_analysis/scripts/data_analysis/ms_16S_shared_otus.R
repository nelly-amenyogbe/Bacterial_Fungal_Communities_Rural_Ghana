######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S: Shared OTUs between mother-infant pairs
######################

# In this script, we identify number of shared OTUs between mother-infant related and unrelated pairs, to determine if related pairs share more OTUs in common.  We plot the most shared OTUS.

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(reshape2)
library(ggplot2)
# helper functions
source("ms_bacteria_fungi_analysis/scripts/functions/function_sig_annot.R")

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_16S_for_analysis.rds")

# Make relative abundance ###
ps.rab <- transform_sample_counts(ps, function(x) x/sum(x))

# get taxonomy
tax <- data.frame(tax_table(ps))

# subset to pairs
sample_data(ps.rab)$pair.yn <- ifelse(is.na(sample_data(ps.rab)$pair.id), "N", "Y")

ps.p <- subset_samples(ps.rab, pair.yn == "Y")
ps.p <- subset_samples(ps.p, group.ages != "DOL 13-17")

# Get shared OTUs 
meta <- data.frame(sample_data(ps.p))

# melt data to have sample x otu
otu.tab <- data.frame(t(otu_table(ps.p)))
otu.tab$MBI.sample.ID <- rownames(otu.tab)

otu.m <- melt(otu.tab, id.vars = "MBI.sample.ID")
otu.m <- filter(otu.m, value > 0)
otu.m <- join(otu.m, meta, by = "MBI.sample.ID")

# get shared

inf.samples <- meta %>% 
  filter(moth.inf == "Child") %>%
  select(stdyid) %>%
  c() %>%
  as.vector()

bm.samples <- meta %>% 
  filter(moth.inf == "Mother", sample.type == "breastmilk") %>%
  select(stdyid) %>%
  c()

ms.samples <- meta %>% 
  filter(moth.inf == "Mother", sample.type == "stool") %>%
  select(stdyid) %>%
  c()

pairs <- meta[,c("stdyid", "pair.id", "partna")]

p.inf <- filter(pairs, as.character(stdyid) %in% inf.samples$stdyid) # 27 infant samples
# choose only those where the partner is also present
p.bm <- filter(p.inf, partna %in% bm.samples$stdyid) # only 20 pairs present
bm.pairs <- unique(as.character(p.bm$pair.id))

p.ms <- filter(p.inf, partna %in% ms.samples$stdyid) # only 23 pairs present
ms.pairs <- unique(as.character(p.ms$pair.id))

# get taxa shared with mom ####
otu.f <- filter(otu.m, value > 1e-4)
otu.f$age.sample <- paste(otu.f$moth.inf, otu.f$sample.type)

# breastmilk
otu.bm <- filter(otu.f, age.sample != "Mother stool", pair.id %in% bm.pairs)

otu.bm.sum <- ddply(otu.bm, c("pair.id", "variable"), summarise, shared = length(variable))
otu.bm.sum$shared <- ifelse(otu.bm.sum$shared == 2, "Y", "N")
otu.bm.sum$pair.otu <- paste(otu.bm.sum$pair.id, otu.bm.sum$variable)

otu.bm$pair.otu <- paste(otu.bm$pair.id, otu.bm$variable)
otu.bm <- join(otu.bm, otu.bm.sum[,c("pair.otu", "shared")], by = "pair.otu")

# add taxonomy
colnames(otu.bm)[2] <-"OTU.num"
otu.bm <- join(otu.bm, tax, by = "OTU.num")

# summarize shared taxa:BM ####
bm.shared <- filter(otu.bm, shared == "Y", moth.inf == "Child")

# how many taxa?
length(unique(bm.shared$Genus)) # 37 Genera
table(as.character(bm.shared$Genus)) # strep, staph, Lacto, Rothia...

# stool shared ####
otu.st <- filter(otu.f, age.sample != "Mother breastmilk", pair.id %in% ms.pairs)

otu.st.sum <- ddply(otu.st, c("pair.id", "variable"), summarise, shared = length(variable))
otu.st.sum$shared <- ifelse(otu.st.sum$shared == 2, "Y", "N")
otu.st.sum$pair.otu <- paste(otu.st.sum$pair.id, otu.st.sum$variable)

otu.st$pair.otu <- paste(otu.st$pair.id, otu.st$variable)
otu.st <- join(otu.st, otu.st.sum[,c("pair.otu", "shared")], by = "pair.otu")

# add taxonomy
colnames(otu.st)[2] <-"OTU.num"
otu.st <- join(otu.st, tax, by = "OTU.num")

# summarize shared taxa:BM ####
st.shared <- filter(otu.st, shared == "Y")

# how many taxa?
length(unique(st.shared$Genus)) # 75 Genera
table(as.character(st.shared$Genus)) # strep, staph, Lacto, Rothia...

# calculate number of shared taxa compared to other pairs ####
# BM and inf stool ####
bm.samples <- as.character(bm.samples$stdyid)
inf.samples <- as.character(inf.samples$stdyid)

# summarize number of shared OTUs between all pairs
bm.shared.sum <- ldply(as.list(inf.samples), function(i){
  
  dat <- filter(otu.bm, stdyid == i)
  
  shared <- ldply(as.list(bm.samples), function(j){
    
    mom.dat <- filter(otu.bm, stdyid == j)
    shared <- length(which(dat$OTU.num %in% mom.dat$OTU.num))
    mom.id <- j
    shared.dat <- data.frame(mom.id, shared)
    shared.dat
    
  })
  
  shared$stdyid <- i
  shared <- shared[,c(3,1,2)]
  shared
  
})

# add metadata to result
bm.shared.sum <- join(bm.shared.sum, meta, by = "stdyid")

# Identify which are related pairs
bm.shared.sum$is.pair <- ifelse(as.character(bm.shared.sum$mom.id) == as.character(bm.shared.sum$partna), "Y", "N")

# get data for stats: BM ####
# shared
bm.inf.s <- filter(bm.shared.sum, moth.inf == "Child", is.pair == "Y")

bm.inf.ns <- filter(bm.shared.sum, moth.inf == "Child", is.pair == "N") # this is where we will need to pull a random samples

set.seed(777) # for reproducibility

select <- sample(c(1:574), size = 20, replace = FALSE)

bm.inf.sample <- bm.inf.ns[c(select),] # this is a random sample of the non-shared pairs

bm.inf.test <- rbind(bm.inf.s, bm.inf.sample) # this is the final dataset to use

# wilcox.test for groups: BM

test.groups <- unique(as.character(bm.inf.test$group.ages))

bm.res <- ldply(as.list(test.groups), function(i){
  
  dat.n <- filter(bm.inf.test, is.pair == "N", group.ages == i)
  dat.y <- filter(bm.inf.test, is.pair == "Y", group.ages == i)
  test <- wilcox.test(dat.n$shared, dat.y$shared)
  p.value <- test$p.value
  group.ages <- i
  res <- data.frame(group.ages, p.value)
  res
  
})

bm.res$p.adj <- 2*bm.res$p.value # adjust manually for two age comparisons
bm.res <- add.sig.annot(bm.res)

# mom stool and inf stool ####
ms.samples <- as.character(ms.samples$stdyid)

ms.shared.sum <- ldply(as.list(inf.samples), function(i){
  
  dat <- filter(otu.st, stdyid == i)
  
  shared <- ldply(as.list(ms.samples), function(j){
    
    mom.dat <- filter(otu.st, stdyid == j)
    shared <- length(which(dat$OTU.num %in% mom.dat$OTU.num))
    mom.id <- j
    shared.dat <- data.frame(mom.id, shared)
    shared.dat
    
  })
  
  shared$stdyid <- i
  shared <- shared[,c(3,1,2)]
  shared
  
})

# add metadata
ms.shared.sum <- join(ms.shared.sum, meta, by = "stdyid")
ms.shared.sum$is.pair <- ifelse(as.character(ms.shared.sum$mom.id) == as.character(ms.shared.sum$partna), "Y", "N")

# get data for stats: stool ####
# shared
ms.inf.s <- filter(ms.shared.sum, moth.inf == "Child", is.pair == "Y") #23

ms.inf.ns <- filter(ms.shared.sum, moth.inf == "Child", is.pair == "N") # this is where we will need to pull a random samples (length = 625)

set.seed(777) # for reproducibility

select <- sample(c(1:625), size = 23, replace = FALSE)

ms.inf.sample <- ms.inf.ns[c(select),] # this is a random sample of the non-shared pairs

ms.inf.test <- rbind(ms.inf.s, ms.inf.sample) # this is the final dataset to use

# wilcox.test for groups: stool

test.groups <- unique(as.character(ms.inf.test$group.ages))

ms.res <- ldply(as.list(test.groups), function(i){
  
  dat.n <- filter(ms.inf.test, is.pair == "N", group.ages == i)
  dat.y <- filter(ms.inf.test, is.pair == "Y", group.ages == i)
  test <- wilcox.test(dat.n$shared, dat.y$shared)
  p.value <- test$p.value
  group.ages <- i
  res <- data.frame(group.ages, p.value)
  res
  
})

ms.res$p.adj <- 2*ms.res$p.value
ms.res <- add.sig.annot(ms.res)

# plot shared taxa ####

# breastmilk
ggplot(bm.inf.test, aes(x = group.ages, y = shared, fill = is.pair)) +
  geom_boxplot(outlier.size = NA) +
  #facet_grid(~group.ages) +
  geom_point(shape = 21, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  theme_classic() +
  scale_fill_manual("Related", values = c("#525252", "#bdbdbd")) +
  labs(x = "", y = "Number of OTUs shared") +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Bacteria: Breastmilk") +
  geom_text(data = bm.res, size = 8, aes(x = group.ages, y = 17, label = sig.annot), inherit.aes = FALSE)

#ggsave("ms_bacteria_fungi_analysis/figures/shared_otus/shared_16s_bm_otus.pdf", device = "pdf", dpi = 300, width = 3, height = 4)

# stool
ggplot(ms.inf.test, aes(x = group.ages, y = shared, fill = is.pair)) +
  geom_boxplot(outlier.size = NA) +
  #facet_grid(~group.ages) +
  geom_point(shape = 21, position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2)) +
  theme_classic() +
  scale_fill_manual("Related", values = c("#525252", "#bdbdbd")) +
  labs(x = "", y = "Number of OTUs shared") +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Bacteria: Stool") +
  geom_text(data = ms.res, size = 6, aes(x = group.ages, y = 160, label = sig.annot), inherit.aes = FALSE)

#ggsave("ms_bacteria_fungi_analysis/figures/shared_otus/shared_16s_ms_otus.pdf", device = "pdf", dpi = 300, width = 3, height = 4)

# get OTUs most often shared ####

# most shared BM ####
# summarize most shared OTUs
bm.most.shared <- data.frame(table(bm.shared$OTU.num))
colnames(bm.most.shared)[1] <- "OTU.num"
bm.most.shared <- join(bm.most.shared, tax, by = "OTU.num")
bm.most.shared <- bm.most.shared[order(bm.most.shared$Freq, decreasing = TRUE),]

bm.ms <- filter(bm.most.shared, Freq > 3) # 10 OTUs

# prepare data for plotting
bm.ms.p <- filter(otu.bm, shared == "Y", OTU.num %in% bm.ms$OTU.num)
bm.ms.p$OTU.num <- factor(bm.ms.p$OTU.num, levels = as.character(bm.ms$OTU.num)) # set factor levels for OTUs

# set aesthetics
bm.ms.p$age.bin <- ifelse(bm.ms.p$group.ages == "Mothers_DOL 0-5", "DOL 0-5",
                          ifelse(bm.ms.p$group.ages == "Mothers_DOL 26-35", "DOL 26-35", as.character(bm.ms.p$group.ages)))

# plot 

# shorten long bacterial names
bm.ms.p$genus.p <- gsub("f__Enterobacteriaceae_unclassified", "f__Enterobacteriaceae", bm.ms.p$Genus)

ggplot(bm.ms.p, aes(x = moth.inf, y = log(value), group = pair.id, fill = age.bin)) +
  geom_point(shape = 21, size = 2) +
  geom_line() +
  scale_fill_manual(values = c("#525252", "#bdbdbd")) +
  facet_wrap(~OTU.num + genus.p, nrow = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Bacteria: Breastmilk") +
  labs(x = "", y = "Log Abundance")

#ggsave("ms_bacteria_fungi_analysis/figures/shared_otus/bm_16s_shared_linegraph.pdf", device = "pdf", dpi = 300, width = 8, height = 5.3)

# most shared stool ####
# summarize most shared OTUs
st.most.shared <- data.frame(table(st.shared$OTU.num))
colnames(st.most.shared)[1] <- "OTU.num"
st.most.shared <- join(st.most.shared, tax, by = "OTU.num")
st.most.shared <- st.most.shared[order(st.most.shared$Freq, decreasing = TRUE),]

st.ms <- filter(st.most.shared, Freq > 8) # 12 OTUs

# prepare data for plotting
st.ms.p <- filter(otu.st, shared == "Y", OTU.num %in% st.ms$OTU.num)
st.ms.p$OTU.num <- factor(st.ms.p$OTU.num, levels = as.character(st.ms$OTU.num)) # set factor levels for OTUs

# set aesthetics
st.ms.p$age.bin <- ifelse(st.ms.p$group.ages == "Mothers_DOL 0-5", "DOL 0-5",
                          ifelse(st.ms.p$group.ages == "Mothers_DOL 26-35", "DOL 26-35", as.character(st.ms.p$group.ages)))

# plot 
# shorten taxa names
unique(st.ms.p$Genus)

st.ms.p$genus.p <- gsub("f__Enterobacteriaceae_unclassified", "f__Enterobacteriaceae", st.ms.p$Genus)
st.ms.p$genus.p <- gsub("f__Coriobacteriaceae_unclassified", "f__Coriobacteriaceae", st.ms.p$genus.p)


ggplot(st.ms.p, aes(x = moth.inf, y = log(value), group = pair.id, fill = age.bin)) +
  geom_point(shape = 21, size = 2) +
  geom_line() +
  scale_fill_manual(values = c("#525252", "#bdbdbd")) +
  facet_wrap(~OTU.num + Genus, nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 0.8),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Bacteria: Stool") +
  labs(x = "", y = "Log Abundance")

#ggsave("ms_bacteria_fungi_analysis/figures/shared_otus/prev_16s_st_shared_linegraph.pdf", device = "pdf", dpi = 300, width = 8, height = 5.3)


