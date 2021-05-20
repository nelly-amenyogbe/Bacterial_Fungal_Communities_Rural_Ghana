# Nelly Amenyogbe
# 17 June 2016
# Useful functions for transforming OTU data for use with sPLS-DA

#### Remove low counts ####
# Function to remove low-count functions
# This function removes any OTUs with counts less than the specified percentage across all samples
# source:  http://mixomics.org/mixmc/prefiltering/

low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

#### Remove NZV ####
# Function to remove near zero variance features
# returns a list of:
   # 1.  New data frame with nzv features removed
   # 2.  A metrics table with the names and positions of all nzv features
   # 3.  The number of features removed

remove.nzv <- function(mat){
  
  df <- data.frame(mat)
  df$rm <- 0
  df.nzv <- nearZeroVar(df) # returns near zero variance features with freq.cutoff = 95/5, and uniqueCut = 10 
  metrics <- df.nzv$Metrics
  nzv.pos <- df.nzv$Position # get column positions of nzv features
  no.removed <- length(nzv.pos)
  
  df <- df[-c(nzv.pos)]
  
 
  rt <- list(df, metrics, no.removed - 1)
  names(rt) <- c("df.no.nzv", "metrics", "no.removed")
  return(rt) # return data frame with nzv feaures removed
  
}


#### TSS transformation ####
# This function normalizes count data by dividing the count of each OTU by the total number of counts in each sample
# source:  http://mixomics.org/mixmc/normalisation/

TSS.transform <- function(count.matrix){
  TSS.divide = function(x){
    x/sum(x)
  }
  
  # function is applied to each row
  data.TSS <- t(apply(count.matrix, 1, TSS.divide))
  
}

#### Centered Log Ratio (CLR) ####
# This function CLR-transforms the supplied data matrix
# Source:  Embedded within mixOmicsv6 source code, under the splsda function... splsda(logratio = "CLR)

# NOTE:  X must be of class "matrix".  A matrix - looking data frame will not do!!

clr <- function (x)
{
  min.value = min(x[which(x != 0)]) * 0.01
  if (dim(x)[2] == 1) {
    res <- list(x.clr = x, gm = rep(1, dim(x)[1]))
  }
  else {
    geometricmean <- function(x) {
      exp(mean(log(x + min.value)))
    }
    gm <- apply(x, 1, geometricmean)
    x.clr <- log((x + min.value)/(gm))
    res <- x.clr
  }
  #class(res) <- "clr"
  return(res)
}

clr.0 <- function (x)
{
  #min.value = min(x[which(x != 0)]) * 0.01
  min.value = 1
  if (dim(x)[2] == 1) {
    res <- list(x.clr = x, gm = rep(1, dim(x)[1]))
  }
  else {
    geometricmean <- function(x) {
      exp(mean(log(x + min.value)))
    }
    gm <- apply(x, 1, geometricmean)
    x.clr <- log((x + min.value)/(gm))
    res <- x.clr
  }
  #class(res) <- "clr"
  return(res)
}


# CSS Transform ####
# require(metagenomeseq)
# for this package, the input is an OTU matrix that's sample columns by OTU rows

CSS.transform <- function(otu.tab){
  
  data.mgSeq <- newMRexperiment(otu.tab, featureData = NULL, libSize = NULL, normFactors = NULL)
  
  # Create CSS data
  p = cumNormStat(data.mgSeq)
  data.cumnorm <- cumNorm(data.mgSeq, p = p)
  data.CSS <- t(MRcounts(data.cumnorm, norm = TRUE, log = TRUE))
  return(data.CSS)
  
}


