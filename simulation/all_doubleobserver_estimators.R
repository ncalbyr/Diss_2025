# Lincoln-Petersen Method (N = (M * S) / R)
# N = abundance estimate
# M = Marked (seen) individuals on occasion 1
# S = Sample size on occasion 2
# R = marked animals (from occasion 1) that are RE-captured on occasion 2

M <- nrow(detection_0[detection_0$obs==1 & detection_0$detect==1, ])  # first occasion marks
S <- nrow(detection_0[detection_0$obs==2, ]) # total number of occasion 2 samples (regardless of detection)
R <- detection_0$detect[detection_0$obs==1]==1 & detection_0$detect[detection_0$obs==2]==1  
R <- length(R[R==TRUE])  # caught by both occasions

lp_0 <- M * S / R

# Random movement
M1 <- nrow(detection_1[detection_1$obs==1 & detection_1$detect==1, ])  # first occasion marks
S1 <- nrow(detection_1[detection_1$obs==2, ]) # total number of occasion 2 samples (regardless of detection)
R1 <- detection_1$detect[detection_1$obs==1]==1 & detection_1$detect[detection_1$obs==2]==1  
R1 <- length(R1[R1==TRUE])  # caught by both occasions

lp_1 <- M1 * S1 / R1

# Compare estimates
lp_0
lp_1

##### OR... #####

chapman.mr <- function(df, mismatch){
  if (mismatch==TRUE){df <- sim.mismatch(df)}
  
  S1 <- nrow(df[df$obs==1 & df$detect==1, ])  # first occasion
  S2 <- nrow(df[df$obs==2 & df$detect==1, ])  # second occasion
  B <- df$detect[df$obs==1]==1 & df$detect[df$obs==2]==1  
  B <- length(B[B==TRUE])  # caught by both occasions
  N.hat <- (S1+1)*(S2+1)/(B+1)-1  # abundance estimate
  var.N <- (S1+1)*(S2+1)*(S1-B)*(S2-B)/(((B+1)^2)*(B+2))
  d <- exp(1.96*sqrt(log(1+(var.N/(N.hat^2)))))
  lcl <- N.hat/d; ucl <- N.hat*d
  return(c(N.hat, lcl, ucl))
}

chap_0 <- chapman.mr(df = detection_0, mismatch = FALSE)
chap_1 <- chapman.mr(df = detection_1, mismatch = FALSE)

# Compare estimates
chap_0
chap_1

### MRDS
library(mrds)
fit2.mrds2 <- function(df, mismatch){
  if(mismatch){
    df <- sim.mismatch(df)
    names(df) <- c( "x","y","object","observer","detected")
  }else{ #change names to be in format for ddf
    names(df) <- c("object","observer", "x","y","detected")
  }
  df$distance <- abs(df$x)
  df$distance[df$distance > 0.03] <- 0.03
  df$Region.Label = rep(1,dim(df)[1])
  df$Sample.Label = rep(1,dim(df)[1])
  try({
    model <- ddf(method = "io", dsmodel =~cds(key ="hr"),
                 mrmodel =~glm(link = "logit", formula = ~distance),
                 data = df, meta.data = list(width = 0.03), control = list(refit = T, nrefit = 5, debug = T))
    ests <- dht(model, region.table = data.frame(Region.Label = 1, Area = 100*0.03*2),
                sample.table = data.frame(Region.Label = 1, Sample.Label = 1,Effort = 100))
    
    N <- ests$individuals$N$Estimate
    lci <- ests$individuals$N$lcl
    uci <- ests$individuals$N$ucl
    return(c(N, lci, uci))
  })
}

# Build objects
mark_dist20 <- fit2.mrds2(df = detection_0, mismatch = FALSE)
mark_dist21 <- fit2.mrds2(df= detection_1, mismatch = FALSE)

# Compare estimates
mark_dist20
mark_dist21
