#################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S: Data QC
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
ps <- readRDS("ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_16S/physeq_16S_raw.rds")

# Contaminants in Stool ####

# create blanks only object
bl <- subset_samples(ps, sample.type == "blank")
bl <- subset_samples(bl, Sample_Plate != 3) # no stools on plate 3
bl.s <- subset_samples(bl, blank.type %in% c("n.pcr.blank", "s.ext.blank"))

# create stools only subject
st <- subset_samples(ps, sample.type == "stool")

# find stool contams ####

ps.blanks <- bl.s
ps.samples <- st

# trim OTU table to taxa present in at least 50% of blanks
b.otutab <- t(data.frame(otu_table(ps.blanks)))
b.otu.f <- remove.low.counts(as.data.frame(b.otutab), 1, 50)
b.otu.f <- b.otu.f$df 
dim(b.otu.f) # 58 OTUs present in at least half of blanks

# For each OTU, calculate geomean and sd 
otu.df <- b.otu.f
otu.df$sample.id <- rownames(otu.df)
otu.m <- melt(otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# summarise the data to contain the gm and stdev of each OTU
otu.summary <- summarise(group_by(otu.m, OTU.num), geomean.bl = geo_mean(1+count), stdev.bl = sd(count))
otu.summary$cutoff <- otu.summary$geomean.bl + otu.summary$stdev.bl

# Summary of sample counts ####
s.bl <- subset_taxa(ps.samples, OTU.num %in% as.character(otu.m$OTU.num))

s.otu.df <- t(data.frame(otu_table(s.bl)))
s.otu.df <- as.data.frame(s.otu.df)
s.otu.df$sample.id <- rownames(s.otu.df)
s.otu.m <- melt(s.otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# add cutoffs to the OTU table
s.otu.m <- join(s.otu.m, otu.summary[,c("OTU.num", "cutoff")], by = "OTU.num")

# Summarize contaminating OTU data
s.otu.summary <- summarise(group_by(s.otu.m, OTU.num),
                           geomean.sample = geo_mean(1+count),
                           stdev.sample = sd(count),
                           prev.sample = length(which(count > 0)) / length(count),
                           prev.below.cutoff  = length(which(count <= cutoff)) / length(count))

# Flag OTUs as contaminats.  These are OTUs where 90% of samples have counts that are below the average count plus one standard deviation in all the blanks
s.otu.summary$flag.contam <- ifelse(s.otu.summary$prev.below.cutoff >= 0.9, TRUE, FALSE)

sum.all <- join(otu.summary[,c("OTU.num","geomean.bl", "stdev.bl", "cutoff")], s.otu.summary, by = "OTU.num")

sum.all <- sum.all[order(sum.all$prev.below.cutoff, decreasing = TRUE),]

stool.contams <- sum.all

rm(ps.blanks, ps.samples)

# Contaminants in breastmilk ####
# create bm blank phyloseq
bm.blanks <- subset_samples(ps, Sample_Plate == 3)
bm.blanks <- subset_samples(bm.blanks, blank.type %in% c("bm.ext.blank", "n.pcr.blank")) 

# create bm samples phyloseq
bm.samples <- subset_samples(ps, sample.type == "breastmilk")

# get contaminating OTUs in breastmilks

ps.blanks <- bm.blanks
ps.samples <- bm.samples

# trim OTU table to taxa present in at least 50% of blanks
b.otutab <- t(data.frame(otu_table(ps.blanks)))
b.otu.f <- remove.low.counts(as.data.frame(b.otutab), 1, 50)
b.otu.f <- b.otu.f$df 
dim(b.otu.f) # 8 OTUs present in at least half of blanks

# For each blank OTU, calculate geomean and sd 
otu.df <- b.otu.f
otu.df$sample.id <- rownames(otu.df)
otu.m <- melt(otu.df, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")

# summarise the data to contain the gm and stdev of each OTU
otu.summary <- summarise(group_by(otu.m, OTU.num),
                         geomean.bl = geo_mean(1+count),
                         stdev.bl = sd(count))

otu.summary$cutoff <- otu.summary$geomean.bl + otu.summary$stdev.bl # the cutoff to flag a contaminant is geometric mean across all the blanks plus one standard deviation

# Summary of sample counts
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

s.otu.summary$flag.contam <- ifelse(s.otu.summary$prev.below.cutoff >= 0.9, TRUE, FALSE)

# join metrics for blanks and samples
sum.all <- join(otu.summary[,c("OTU.num","geomean.bl", "stdev.bl", "cutoff")], s.otu.summary, by = "OTU.num")

sum.all <- sum.all[order(sum.all$prev.below.cutoff, decreasing = TRUE),]

bm.contams <- sum.all

# Find contaminant OTUs to remove from all data ####
otu.contam.stools <- filter(stool.contams, flag.contam == TRUE)
otu.contam.bm <- filter(bm.contams, flag.contam == TRUE)

otu.contams <- Reduce(intersect, list(otu.contam.stools$OTU.num, otu.contam.bm$OTU.num))

# add taxonomy
tax <- data.frame(tax_table(ps))
tax.contams <- filter(tax, OTU.num %in% otu.contams)

# Two OTUs:  Otu00088_Halomonas and Otu00137_Shewanella

# Check the prevalence of other stool contaminants in the breastmilk data

bm.inv <- subset_taxa(bm.samples, OTU.num %in% otu.contam.stools$OTU.num) # 19 taxa in breastmilk are flagged as contaminants in the stool data

# Summarize the counts of these OTUs across the breastmilk samples
otab <- t(data.frame(otu_table(bm.inv)))
otab <- data.frame(otab)
otab$sample.id <- rownames(otab)
otab.m <- melt(otab, id.vars = "sample.id", variable.name = "OTU.num", value.name = "count")
otab.sum <- summarise(group_by(otab.m, OTU.num),
                      prevalence = length(which(count > 0)) / length(count)) # prevalence is the fraction of samples in which the OTU is detected
otab.sum

otab.contam <- filter(otab.sum, prevalence < 0.1) # these are OTUs present in fewer than 10% of breastmilk samples

otab.contam # 14 OTUs present in 0 to 7% of breastmilk samples

all.otus.contam <- c(as.character(otab.contam$OTU.num), otu.contams) #14 total OTUs to remove

# Filter Phyloseqs ####
# Because the breastmilk and stool data will be analyzed together, unless an OTU is a contaminant in both sample types, it cannot be removed from the dataset completely. Thus, for OTUs that demonstrate low levels of contamination among breastmilk samples, we will subtract the average count of those OTUs from the counts in the individual samples.

# make a list of contaminant OTUs
otus.all <- as.character(rownames(otu_table(ps)))
otu.keep <- otus.all[-which(otus.all %in% all.otus.contam)]

# Adjust sample counts function ####

# This function take a phyloseq object, and removed OTU counts based on their levels in contaminata
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

# breastmilk
bm.filter <- subset_taxa(bm.samples, OTU.num %in% otu.keep)

bm.depth <- colSums(otu_table(bm.samples))
range(bm.depth) # no samples below 1000 reads

# Remove background contaminants 
bm.cutoff <- bm.contams[,c("OTU.num", "cutoff")]
bm.cutoff <- bm.cutoff[-which(bm.cutoff$OTU.num %in% all.otus.contam),]
bm.cutoff <- bm.cutoff[order(bm.cutoff$OTU.num),]

bm.physeq.final <- adjust.physeq.counts(bm.filter, cutoff.dataframe = bm.cutoff)

# stool
stool.filter <- subset_taxa(st, OTU.num %in% otu.keep)
stool.depth <- colSums(otu_table(stool.filter))
length(which(stool.depth < 1000)) # 4 samples with low sequencing depth

# Samples with low counts will be removed at a later stage, after comparison to fungal data

# Remove counts of contaminant OTUs
st.cutoff <- stool.contams[,c("OTU.num", "cutoff")]
st.cutoff <- st.cutoff[-which(st.cutoff$OTU.num %in% all.otus.contam),]
st.cutoff <- st.cutoff[order(st.cutoff$OTU.num),]

st.physeq.final <- adjust.physeq.counts(stool.filter, cutoff.dataframe = st.cutoff)

# combine both subsets ####
physeq.final <- merge_phyloseq(st.physeq.final, bm.physeq.final)

# remove extra metadata
df <- data.frame(sample_data(physeq.final))
colnames(df)
# remove columns "X" and "blank.type"
sample_data(physeq.final) <- sample_data(physeq.final)[,-which(colnames( sample_data(physeq.final)) %in% c("X", "blank.type"))]

check <- data.frame(sample_data(physeq.final))

# save relevant data ####
#saveRDS(physeq.final, "ms_bacteria_fungi_analysis/Rdata/R_export/prev_16S_physeq_filtered.rds")


