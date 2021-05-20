#################
# Nelly Amenyogbe
# 27-March-2019
# Function:  add significance annotation 
################

# Input data frame with p.adj column, and will return annotated significance column useful for graphing.

add.sig.annot <- function(dat){
  
  dat$sig.annot <- ifelse(dat$p.adj <= 0.0001, "****",
                             ifelse(dat$p.adj <= 0.001, "***",
                                    ifelse(dat$p.adj <= 0.01, "**",
                                           ifelse(dat$p.adj <= 0.05, "*",
                                                  ifelse(dat$p.adj > 0.05, "ns", dat$p.adj)))))
  
  return(dat)

}