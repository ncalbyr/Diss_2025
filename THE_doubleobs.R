library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(spatstat)
library(ggplot2)
library(truncnorm)
library(grid)
library(gridExtra)
library(circular)
library(fields)
##### 1) SIMULATE 210 ANIMALS
# (from perpendicular distribution of "primate.dat")
# This object uses pi.norm density distribution for pi.x
Fit.n.ip0$fit$par

hist(primate.dat$x)

# Parameters for "pi.x" distribution from "Fit.n.ip0" in file "primate.fit.R"
beta1 <- Fit.n.ip0$fit$par[1]
beta2 <- Fit.n.ip0$fit$par[2]

lphi1 <- Fit.n.ip0$fit$par[3]
lphi2 <- Fit.n.ip0$fit$par[4]

# Generate x-values from the known perpendicular distribution
x_set <- rtruncnorm(n = 210, # number of items to generate
                    a = 0, # lower bound
                    b = 0.05, # upper bound
                    mean = lphi1, sd = exp(lphi2) # parameters for normal distribution
)
openGraph(h=4,w=9) 
par(mfrow=c(1,2))
hist(x_set)
plot(x_set)

##### 2) FIND FUNCTION (WITH HAZARD SIMULATEOR) TO PRODUCE Y's
# Generate forward distances (y) ONLY for detected individuals
detected_xy <- detect2DLT(x = x_set, hr = "ip0", b = c(beta1,beta2),
                          ystart = 0.05, ny = 1000,
                          getIDs = T)
head(detected_xy)
hist(detected_xy$x)
length(detected_xy$x)
##### 4) MOVE X&Y's TO X2&Y2's (Add them back into the fitted object?)

moving <- move.data(df = detected_xy,
                    move = 0,
                    keep_angle = F)
# make sure the NEW version of move.data is loaded in
unmoved <- move.data(df = detected_xy,
                     move=2,
                     keep_angle = F)
##### 5) USE DETECTION FUNCTION ON NEW LOCATIONS
# Build function
detect.data <- function(sigma) {
  
  # Create long-format data
  df <- data.frame(
    id = rep(1:nrow(moving), each = 2),
    obs = rep(1:2, times = nrow(moving)),
    x = c(moving$x, moving$x2),
    y = c(moving$y, moving$y2),
    detect = NA
  )
  
  # b. Simulate detection (only needed for observation 2)
  ys <- seq(0, 0.05, length.out=100)  
  # Generate detection probability with this function from LT2D package
  obs2.probs <- p.approx(ys, df$x[df$obs==2], ip0, b=c(beta1, beta2), what = "px")
  
  df$detect[df$obs==2] <- rbinom(n, 1, obs2.probs)  # second observer detection
  df$detect[df$obs==1] <- 1 # Because all obs. 1's are already known detections
  
  df$detect[abs(df$x)>0.05] <- 0 # Anything beyond observable range goes undetected
  # return dataset
  return(df)
}


for_dobs <- detect.data(sigma = 0.2)
head(for_dobs)

##### 6) LOOP SEVERAL TIMES TO PRODUCE DISTRIBUTIONS OF MEAN ESTIMATES

# Double-Observer Estimates of Abundance
chapman.mr <- function(df, mismatch){
  if (mismatch==TRUE){df <- mismatch(df)}
  
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

chap_F <- chapman.mr(df = for_dobs, mismatch = FALSE)
chap_T <- chapman.mr(df = for_dobs, mismatch = TRUE)

# Compare estimates
chap_F
chap_T

nrow(for_dobs[for_dobs$obs==2 & for_dobs$detect==1, ])
nrow(for_dobs[for_dobs$obs==2, ])
nrow(for_dobs)