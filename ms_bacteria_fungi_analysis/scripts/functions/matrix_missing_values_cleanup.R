#### Nelly Amenyogbe
#### 07 December 2016
#### Functions to remove columns/rows with missing values and replace NAs with column means

#### Column remove List ####
beads.rm <- function(mat, percent.cutoff){
  beads.na <- c()
  beads <- colnames(mat)
  for(i in 1:ncol(mat)){
    beads.na[i] <- 100*length(which(is.na(mat[,i])))/length(mat[,i])
  }
  names(beads.na) <- beads
  beads.rm <- names(which(beads.na > percent.cutoff))
  beads.rm.pos <- which(colnames(mat) %in% beads.rm)
  df <- data.frame(beads.rm, beads.rm.pos)
  df
}  # creates a list of beads that were removed

#### Remove columns ####
remove.beads <- function(a, percent.cutoff){
  beads.na <- c()
  a$remove <- NA
  beads <- colnames(a)
  for(i in 1:ncol(a)){
    beads.na[i] <- 100*length(which(is.na(a[,i])))/length(a[,i])
  }
  names(beads.na) <- beads
  beads.rm <- names(which(beads.na > percent.cutoff))
  beads.rm.pos <- which(colnames(a) %in% beads.rm)
  a <- a[,-c(beads.rm.pos)]
  a
} # removes beads

#### Remove columns with threshold low values

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


#### List rows to remove ####
which.sub.rm <- function(a, percent.cutoff){
  subj <- rownames(a)
  sub.na <- c()
  
  for(i in 1:nrow(a)){
    sub.na[i] <- 100*length(which(is.na(a[i,])))/length(a[i,])
  }
  names(sub.na) <- subj
  subj.rm <- names(which(sub.na > percent.cutoff))
  subj.rm.pos <- which(rownames(a) %in% subj.rm)
  df <- data.frame(subj.rm, subj.rm.pos)
  df
}

#### Remove rows ####

subj.rm <- function(a, percent.cutoff){
  sub.na <- c()
  a[nrow(a) + 1,] <- NA
  subj <- rownames(a)
  for(i in 1:nrow(a)){
    sub.na[i] <- 100*length(which(is.na(a[i,])))/length(a[i,])
  }
  names(sub.na) <- subj
  subj.rm <- names(which(sub.na > 30))
  subj.rm.pos <- which(rownames(a) %in% subj.rm)
  a <- a[-c(subj.rm.pos),]
  a
}

#### Replace NAs ####

na.replace <- function(df){
  f1 <- function(vec) { 
    m <- mean(vec, na.rm = TRUE) 
    vec[is.na(vec)] <- m 
    return(vec) 
  } 
  df.narm <- apply(df, 2, f1)
  colnames(df.narm) <- colnames(df)
  rownames(df.narm) <- rownames(df)
  df.narm <- as.matrix(df.narm)
} # replaces NA with column mean