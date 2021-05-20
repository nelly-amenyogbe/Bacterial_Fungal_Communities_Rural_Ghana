#### Filtering Function ####
# Author:  Pedro Dimitrui
# 22-September-2017

simpletrim2 = function(physeq, minobs, nreads) {
  Ji = nsamples(physeq)
  # Force orientation to be sample-by-OTU
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  # `prevalence` is the fraction of total samples in which an OTU is observed
  # at least `minobs` times.
  prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
    return(sum(x > minobs))
  }, minobs)/(Ji)
  # Will only keep OTUs that appear in more than X% of samples and have total
  # reads greater than `nreads`
  keepOTUs = prevalence > 0.01 & taxa_sums(physeq) >= nreads
  return(prune_taxa(keepOTUs, physeq))
}


#### Remove low counts ####

remove.low.counts <- function(a, abundance.cutoff, percent.cutoff){
  
  beads.low <- c()
  a$remove <- NA
  beads <- colnames(a)
  for(i in 1:ncol(a)){
    beads.low[i] = 100*length(which(a[,i] >= abundance.cutoff )) / length(a[,i])
    #beads.low[i] <- 100*length(which(a[,i] <= percent.cutoff))/length(a[,i])
  }
  names(beads.low) <- beads
  beads.rm <- names(which(beads.low <= percent.cutoff))
  beads.rm <- c(beads.rm, "remove")
  beads.rm.pos <- which(colnames(a) %in% beads.rm)
  a <- a[,-c(beads.rm.pos)]
  data <- list(a, ncol(a))
  names(data) <- c("df", "otus.keep")
  return(data)
}

# Geometric Mean ####
geo_mean <- function(x){
  exp(mean(log(x)))
} 


