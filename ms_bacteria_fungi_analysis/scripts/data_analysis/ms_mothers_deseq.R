######################
# Nelly Amenyogbe
# Ghana early-life bacteria-fungi Manusript analysis
# 16S and ITS2: DESeq2 for maternal microbiomes
######################

# In this script, we use DESeq2 to identify differentially abundant OTUs between mothers one week (Mother_0-5 days) and one month (Mothers_26-35 days) post partum for stool and breastmilk.

# Load packages
library(phyloseq)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(DESeq2)
# helper functions
source("ms_bacteria_fungi_analysis/scripts/functions/functions_otu_transformations.R")
source("ms_bacteria_fungi_analysis/scripts/functions/matrix_missing_values_cleanup.R")

# load data
s16 <-  readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_16S_for_analysis.rds")
its <-  readRDS("ms_bacteria_fungi_analysis/Rdata/R_export/ps_ITS2_for_analysis.rds")

# subset mothers bm and stool ####
its.bm <- subset_samples(its, sample.type == "breastmilk")
s16.bm <- subset_samples(s16, sample.type == "breastmilk")

# subset mothers stool
moms <- c("Mothers_DOL 0-5", "Mothers_DOL 26-35")

its.stool <- subset_samples(its, sample.type == "stool")
its.stool <- subset_samples(its.stool, group.ages %in% moms)

s16.stool <-  subset_samples(s16, sample.type == "stool")
s16.stool <- subset_samples(s16.stool, group.ages %in% moms)

# make list of physeqs ####
physeqs <- list(its.bm, s16.bm, its.stool, s16.stool)
names(physeqs) <- c("its.bm", "s16.bm", "its.stool", "s16.stool")

# Filter OTU table to only include abundant OTUs ####
# Create matrix of OTU data

otu.tabs <- llply(physeqs, function(i){
  
  tab <- t(data.frame(otu_table(i)))
  tab
  
})

# retain only OTUs present in at least 5% of samples
otu.f <- llply(otu.tabs, function(i){
  
  otu.f <- remove.low.counts(as.data.frame(i), 3, 5)
  otu.f <- otu.f$df
  otu.f
  
})

llply(otu.f, function(i){dim(i)})

# filter physeq objects to include relevant OTUs only
its.bm.f <- subset_taxa(its.bm, OTU.num %in% colnames(otu.f$its.bm))
its.st.f <- subset_taxa(its.stool, OTU.num %in% colnames(otu.f$its.stool))
s16.bm.f <- subset_taxa(s16.bm, OTU.num %in% colnames(otu.f$s16.bm))
s16.st.f <- subset_taxa(s16.stool, OTU.num %in% colnames(otu.f$s16.stool))

ps.f <- list(its.bm.f, its.st.f, s16.bm.f, s16.st.f)
names(ps.f) <- c("its.bm", "its.st", "s16.bm", "s16.st")

# Run DESeq2 ####
# Create DESeq objects
ds.obj <- llply(ps.f, function(i){
  
  design <- ~group.ages
  
  dds <- phyloseq_to_deseq2(i, design)
  
  # because for most data every row has a zero, need to manually specify geomeans.  
  # source: https://support.bioconductor.org/p/62246/#62250
  cts <- counts(dds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
  dds
  
})

# run test
ds.tests <- llply(ds.obj, function(i){
  
  ds.wald <- DESeq(i, test = "Wald", fitType = "parametric")
  ds.wald
  
})

# Get deseq results ####
taxtab.its <- data.frame(tax_table(its))
taxtab.s16 <- data.frame(tax_table(s16))

ds.res.its <- llply(ds.tests[1:2], function(i){
  
  res <- data.frame(results(i))
  res$OTU.num <- rownames(res)
  res <- join(res, taxtab.its, by = "OTU.num")
  res
  
})

ds.res.s16 <- llply(ds.tests[3:4], function(i){
  
  res <- data.frame(results(i))
  res$OTU.num <- rownames(res)
  res <- join(res, taxtab.s16, by = "OTU.num")
  res
  
})

ds.res <- c(ds.res.its, ds.res.s16)
names(ds.res)

# get significant results only
sigtabs <- llply(ds.res, function(i){
  
  sigtab <- filter(i, padj < 0.1)
  sigtab
  
})

llply(sigtabs, function(i){dim(i)})

# plot results ####

# Function to plot significant results, using supplied sigtab, and character for desired title
ds.genus.logchange <- function(sigtab, title){
  
  # Arrange data
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels = names(x))
  
  # Plot
  p <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, fill = Phylum)) +
    geom_point(size = 6, alpha = 0.8, color = "black", shape = 21) +
    theme_bw() +
    #scale_color_manual("Phylum", values = phy.cols) +
    theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(face = "bold"),
          axis.line = element_line(size = 0.8),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    ggtitle(title) +
    labs(x = "", y = "Log2 Fold-Change") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#08519c")
  
  p
}


# Mothers 16S stool
ds.genus.logchange(sigtabs$s16.st, "Mothers: 16S Stool")

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/deseq_mot_16s_stool.pdf", device = "pdf", dpi = 300, width = 8, height = 8)

# Mothers 16S BM
ds.genus.logchange(sigtabs$s16.bm, "Mothers: 16S Breastmilk")

#ggsave("ms_bacteria_fungi_analysis/figures/taxonomy/deseq_mot_16s_bm.pdf", device = "pdf", dpi = 300, width = 7, height = 2.5)

# Mothers ITS2 BM
ds.genus.logchange(sigtabs$its.bm, "ITS2 Breastmilk") # Only two OTUs.  Do not save.

# Mothers ITS2 stool
ds.genus.logchange(sigtabs$its.st, "ITS2 Stool") # None

# END ####

# SessionInfo ####
sessionInfo()

#R version 3.3.2 (2016-10-31)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: macOS  10.15.4

#locale:
#  [1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

#attached base packages:
#  [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
#[9] base     

#other attached packages:
#  [1] DESeq2_1.14.1              SummarizedExperiment_1.4.0 Biobase_2.34.0            
#[4] GenomicRanges_1.26.4       GenomeInfoDb_1.10.3        IRanges_2.8.2             
#[7] S4Vectors_0.12.2           BiocGenerics_0.20.0        reshape2_1.4.3            
#[10] ggplot2_3.1.0              dplyr_0.8.0.1              plyr_1.8.4                
#[13] phyloseq_1.19.1 
