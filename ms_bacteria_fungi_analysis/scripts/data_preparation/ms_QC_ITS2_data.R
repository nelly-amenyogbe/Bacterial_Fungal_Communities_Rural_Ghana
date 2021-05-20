#################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Data QC
################

# In this script, we prepare the data for biological analysis by:
# 1. Flagging contaimating OTUs and removing them from analysis
# 2. Trimming the blank samples from the dataset
# Breastmilk and stool samples are processed separately

# load packages
library(plyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(reshape2)

# helper function
source("ms_bacteria_fungi_analysis/scripts/functions/OTU_filter_function.R")

# load data
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_ITS2/physeq_ITS2_raw.rds")

sample_data(ps)$Sample_Plate <- as.character(sample_data(ps)$Sample_Plate)

# subset data to blanks vs samples ####

# create blanks only phyloseq
# stool blanks
ps.bl.s <- subset_samples(ps, Sample_Plate != "2")
ps.bl.s <- subset_samples(ps.bl.s, sample.type == "blank")

# bm blanls
ps.bl.bm <- subset_samples(ps, sample.type %in% c("blank", "bm.blank"))
ps.bl.bm <- subset_samples(ps.bl.bm, Sample_Plate == "2")

# samples only phyloseq
# create stools only subject
ps.st <- subset_samples(ps, sample.type == "stool")
ps.bm <- subset_samples(ps, sample.type == "breastmilk")

# Contaminants in Stool ####

ps.blanks <- ps.bl.s
ps.samples <- ps.st

# trim OTU table to taxa present in at least 50% of blanks
b.otutab <- t(data.frame(otu_table(ps.blanks)))
b.otu.f <- remove.low.counts(as.data.frame(b.otutab), 1, 50) # find OTUs present in at least half of blanks
b.otu.f <- b.otu.f$df 
dim(b.otu.f) # 6 possible contaminants

# For each OTU, calculate geomean and sd 

otu.df <- b.otu.f
otu.df$sample.id <- rownames(otu.df)
otu.m <- melt(otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# summarise the data to contain the gm and stdev of each OTU
otu.summary <- summarise(group_by(otu.m, OTU.num),
                         geomean.bl = geo_mean(1+count),
                         stdev.bl = sd(count))

otu.summary$cutoff <- otu.summary$geomean.bl + otu.summary$stdev.bl # the cutoff is the geomean of sample counts in the blanks plus one standard deviation.

# Summarize contaminating OTU counts in the samples
s.bl <- subset_taxa(ps.samples, OTU.num %in% as.character(otu.m$OTU.num))

s.otu.df <- t(data.frame(otu_table(s.bl)))
s.otu.df <- as.data.frame(s.otu.df)
s.otu.df$sample.id <- rownames(s.otu.df)
s.otu.m <- melt(s.otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# add cutoffs to the OTU table
s.otu.m <- join(s.otu.m, otu.summary[,c("OTU.num", "cutoff")], by = "OTU.num")

# summarize the number of samples per OTU below cutoff
s.otu.summary <- summarise(group_by(s.otu.m, OTU.num),
                           geomean.sample = geo_mean(1+count),
                           stdev.sample = sd(count),
                           prev.sample = length(which(count > 0)) / length(count),
                           prev.below.cutoff  = length(which(count <= cutoff)) / length(count))

# Flag OTUs as contaminats.  These are OTUs where 90% of samples have counts that are below the average count plus one standard deviation in all the blanks
s.otu.summary$flag.contam <- ifelse(s.otu.summary$prev.below.cutoff >= 0.9, TRUE, FALSE)

sum.all <- join(otu.summary[,c("OTU.num","geomean.bl", "stdev.bl", "cutoff")], s.otu.summary, by = "OTU.num")

sum.all <- sum.all[order(sum.all$prev.below.cutoff, decreasing = TRUE),]

# add taxa
tax <- data.frame(tax_table(ps))
sum.all <- join(sum.all, tax[,c("OTU.num", "Genus")], by = "OTU.num")

# final stool file 
stool.contams <- sum.all

rm(ps.blanks, ps.samples)
# Breastmilk contams ####

ps.blanks <- ps.bl.bm
ps.samples <- ps.st

# trim OTU table to taxa present in at least 50% of blanks
b.otutab <- t(data.frame(otu_table(ps.blanks)))
b.otu.f <- remove.low.counts(as.data.frame(b.otutab), 1, 50) # find OTUs present in at least half of blanks
b.otu.f <- b.otu.f$df 
dim(b.otu.f) # 3 possible OTUs

# For each OTU, calculate geomean and sd 
otu.df <- b.otu.f
otu.df$sample.id <- rownames(otu.df)
otu.m <- melt(otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")
plot(otu.m$count) # one extreme outlier of over 3000.  Remove, as this will sqew the data
otu.m <- filter(otu.m, count < 1000)

# summarise the data to contain the gm and stdev of each OTU
otu.summary <- summarise(group_by(otu.m, OTU.num),
                         geomean.bl = geo_mean(1+count),
                         stdev.bl = sd(count))

otu.summary$cutoff <- otu.summary$geomean.bl + otu.summary$stdev.bl

# Summary of sample counts ####
s.bl <- subset_taxa(ps.samples, OTU.num %in% as.character(otu.m$OTU.num))

s.otu.df <- t(data.frame(otu_table(s.bl)))
s.otu.df <- as.data.frame(s.otu.df)
s.otu.df$sample.id <- rownames(s.otu.df)
s.otu.m <- melt(s.otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# add cutoffs to the OTU table
s.otu.m <- join(s.otu.m, otu.summary[,c("OTU.num", "cutoff")], by = "OTU.num")

s.otu.summary <- summarise(group_by(s.otu.m, OTU.num),
                           geomean.sample = geo_mean(1+count),
                           stdev.sample = sd(count),
                           prev.sample = length(which(count > 0)) / length(count),
                           prev.below.cutoff  = length(which(count <= cutoff)) / length(count))

# Flag OTUs as contaminats.  These are OTUs where 90% of samples have counts that are below the average count plus one standard deviation in all the blanks
s.otu.summary$flag.contam <- ifelse(s.otu.summary$prev.below.cutoff >= 0.9, TRUE, FALSE)

sum.all <- join(otu.summary[,c("OTU.num","geomean.bl", "stdev.bl", "cutoff")], s.otu.summary, by = "OTU.num")

sum.all <- sum.all[order(sum.all$prev.below.cutoff, decreasing = TRUE),]

# add taxa
sum.all <- join(sum.all, tax[,c("OTU.num", "Genus")], by = "OTU.num")

bm.contams <- sum.all

# check for stool contams in breastmilk ####
st.contam.otus <- stool.contams$OTU.num
st.true <- filter(stool.contams, flag.contam == "TRUE")

ps.st.contams <- subset_taxa(ps, OTU.num %in% st.contam.otus)

p.st <- plot_bar(ps.st.contams, x = "MBI.sample.ID", fill = "OTU.num") + facet_wrap(~sample.type, scales = "free")
p.st

p.st.dat <- p.st$data
p.st.dat$flag.contam <- ifelse(p.st.dat$OTU.num %in% st.true$OTU.num, "TRUE", "FALSE")

# plot distribution of stool contaminants across sample types
p.st.contams <- ggplot(p.st.dat, aes(x = MBI.sample.ID, y = Abundance, fill = OTU.num)) +
  geom_bar(stat = "identity") +
  facet_wrap(sample.type~flag.contam, scales = "free")

p.st.contams # TRUE stool contams appear very sporadically in bm samples

# get overlap between stool and bm contams ####
bm.true <- filter(bm.contams, flag.contam == "TRUE")
st.true <- filter(stool.contams, flag.contam == "TRUE")

contams.overlap <- Reduce(intersect, list(st.true$OTU.num, bm.true$OTU.num))
contams.overlap # no overlap.  Conclude:  since the bm contams do not appear in stool, remove from dataset.  Stool contaminants appear rarely in breastmilk samples, and at low counts under ~10.  Thus, it makes sense to remove all bm and stool TRUE contaminants from the dataset entirely.  For the "FALSE" contaminants, will still remove their mean counts from the samples.  

# remove contaminating OTUs 
true.contams <- c(as.character(st.true$OTU.num), as.character(bm.true$OTU.num)) # 6 OTUs will be removed.  

all.otus <- as.character(rownames(otu_table(ps)))
otus.keep <- all.otus[-which(all.otus %in% true.contams)]

# adjust physeq counts ####
# for the remaineder of samples, remove the mean + SD counts of false contaminants from sample counts in the dataset.

adjust.physeq.counts <- function(physeq, cutoff.dataframe){
  
  cutoff.df <- cutoff.dataframe
  
  # Adjust stool sample counts ####
  # create new OTU table with counts subtracted
  orig <- data.frame(otu_table(physeq))
  
  # separate out OTUs to filter from clean OTUs
  orig.short <- orig[rownames(orig) %in% cutoff.df$OTU.num,]
  orig.rest <- orig[-which(rownames(orig) %in% cutoff.df$OTU.num),]
  
  for(i in 1:nrow(orig.short)){
    orig.short[i,] <- orig.short[i,] - cutoff.df[i,]$cutoff
  }
  
  # change negatives to zero
  orig.short <- replace(orig.short, orig.short < 0, 0)
  orig.short <- round(orig.short)
  
  # recreate full table
  orig.recom <- rbind(orig.rest, orig.short)
  orig.recom <- orig.recom[order(rownames(orig.recom)),]
  
  # Replace values in phyloseq otu table 
  for(i in 1:ncol(otu_table(physeq))){
    otu_table(physeq)[,i] <- orig.recom[,i]
  }
  
  return(physeq)
}

# adjust breastmilk counts ####

# Filter BM samples ####
bm.cutoff.df <- filter(bm.contams, flag.contam == "FALSE")
bm.cutoff.df <- bm.cutoff.df[,c("OTU.num", "cutoff")]

bm.samples <- c(as.character(sample_data(ps.bm)$MBI.sample.ID), as.character(sample_data(ps.bl.bm)$MBI.sample.ID))
bm.ps <- subset_samples(ps, MBI.sample.ID %in% bm.samples)

bm.filter <- adjust.physeq.counts(bm.ps, bm.cutoff.df)
bm.filter <- subset_taxa(bm.filter, OTU.num %in% otus.keep)

# check depth
depth <- colSums(otu_table(ps.bm))
length(which(depth < 1000)) # 4 samples under 1000 reads to remove after comparison with 16S data

bm.final <-  subset_samples(bm.filter, sample.type == "breastmilk")

# adjust stool counts ####
st.cutoff.df <- filter(stool.contams, flag.contam == "FALSE")
st.cutoff.df <- st.cutoff.df[,c("OTU.num", "cutoff")]

st.samples <- c(as.character(sample_data(ps.st)$MBI.sample.ID), as.character(sample_data(ps.bl.s)$MBI.sample.ID))
st.ps <- subset_samples(ps, MBI.sample.ID %in% st.samples)

# Filter OTUs by subtracting OTU counts for contams that cannot be removed, and remove 4 breastmilk the contaminating OTUs from the dataset 
st.filter <- adjust.physeq.counts(st.ps, st.cutoff.df)
st.filter <- subset_taxa(st.filter, OTU.num %in% otus.keep)

st.final <- subset_samples(st.filter, sample.type == "stool")

# combine final ps ####
physeq.final <- merge_phyloseq(st.final, bm.final)

# save data ####
#saveRDS(physeq.final, "ms_bacteria_fungi_analysis/Rdata/R_export/prev_ITS2_physeq_filtered.rds")
