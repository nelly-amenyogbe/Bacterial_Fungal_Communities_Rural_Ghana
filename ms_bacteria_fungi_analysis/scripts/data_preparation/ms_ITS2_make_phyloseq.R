#################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# ITS2: Create phyloseq object
################

# In this script, we take the output OTU table (.shared) and taxonomy file (.taxonomy) from MOTHUR processing of the fastq files, and combine these with sample metadata to produce a phyloseq object.  This data-friendly format will ease the process of performing any analyses, including QA, on these data.

# load packages
library(plyr)
library(dplyr)
library(phyloseq)

# Import MOTHUR

physeq <- import_mothur(mothur_shared_file = "ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_ITS2/its2.stability.trim.contigs.unique.good.precluster.pick.pick.agc.unique_list.shared", mothur_constaxonomy_file = "ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_ITS2/its2.stability.trim.contigs.unique.good.precluster.pick.pick.agc.unique_list.0.03.cons.taxonomy")

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Get sample data and metadata
sample.data <- read.csv("ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_ITS2/its2_sample_names.csv")
meta <- read.csv("ms_bacteria_fungi_analysis/Rdata/sample_metadata.csv")
meta$X <- NULL

# Add metadata to sample data
colnames(sample.data)
colnames(sample.data)[which(colnames(sample.data) == "subject.id")] <- "stdyid"
sample.data <- join(sample.data, meta, by = "stdyid")

# change column name for sample to match 16S data
colnames(sample.data)[which(colnames(sample.data) == "Sample_Name")] <- "MBI.sample.ID"

# get list of samples in OTU data
samples <- colnames(otu_table(physeq)) # 227 samples

length(which(samples %in% sample.data$MBI.sample.ID)) # all accounted for

# Add sample data to phyloseq object ####
rownames(sample.data) <- sample.data$MBI.sample.ID

# check that samples are in the same order as they are in the sequencing data
length(which(rownames(sample.data) ==  samples)) # not in the same order

sample.dat.2 <- sample.data[match(samples, sample.data$MBI.sample.ID),] #re-order

# check
length(which(rownames(sample.dat.2) == samples)) # 227

# now add to pyloseq object
sampledat <- sample_data(sample.dat.2)
physeq <- merge_phyloseq(physeq, sampledat)
physeq

# remove blanks from plate 1, since there were no biological samples from this study on that plate.
physeq <- subset_samples(physeq, Sample_Plate != 1)

# create OTU label in tax tab
tax <- data.frame(tax_table(physeq))
tax$OTU.num <- rownames(tax)
tax <- as.matrix(tax)
taxtab <- tax_table(tax)

physeq <- merge_phyloseq(physeq, taxtab)

# save relevant file ####
#saveRDS(physeq, "ms_bacteria_fungi_analysis/Rdata/raw_data/MOTHUR_ITS2/physeq_ITS2_raw.rds")
