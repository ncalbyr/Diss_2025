library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
###
library(ggplot2)
library(grid)
library(gridExtra)
library(circular)
library(ggplot2)
library(spatstat)
library(truncnorm)
library(fields)
library(dplyr)

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
png(filename = "detected_xy.png",
    width = 500, height = 500)
grid.table(head(detected_xy))
dev.off()
##### 4) MOVE X&Y's TO X2&Y2's (Add them back into the fitted object?)
##### 5) USE DETECTION FUNCTION ON NEW LOCATIONS
# Build function

moving <- move.data(df = detected_xy, move = 0, keep_angle = F)

detect.data <- function(sigma) {
  
  n <- nrow(df)
  
  # Create long-format data
  df <- data.frame( # this isn't alternating x and x2's, it's putting them all in one long line (all x's and THEN all x2's)
    id = rep(1:nrow(moving), each = 2),
    obs = rep(1:2, times = nrow(moving)),
    x = as.vector(rbind(moving$x, moving$x2)),
    y = as.vector(rbind(moving$y, moving$y2)),
    detect = NA
  )
  
  # b. Simulate detection (only needed for observation 2)
  ys <- seq(0, 0.05, length.out=100)  
  # Generate detection probability with this function from LT2D package
  obs2.probs <- p.approx(ys, df$x[df$obs==2], ip0, b=c(beta1, beta2), what = "px")
  
  df$detect[df$obs==2] <- rbinom(length(obs2.probs), 1, obs2.probs)  # second observer detection
  df$detect[df$obs==1] <- 1 # Because all obs. 1's are already known detections
  
  df$detect[abs(df$x)>0.05] <- 0 # Anything beyond observable range goes undetected
  # return dataset
  return(df)
}


for_dobs <- detect.data(sigma = 0.2)
head(for_dobs)

# Make an image of this table for the report
png(filename = "two_occasion_detected.png",
    width = 480, height = 480)
grid.table(head(for_dobs))
dev.off()

##### 6) LOOP SEVERAL TIMES TO PRODUCE DISTRIBUTIONS OF MEAN ESTIMATES

# Double-Observer Estimates of Abundance
## (sourced from LT2Dcal/inst/report/report.simulation.R)
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
  return(c(N.hat, lcl, ucl, S1, S2, B))
}

chap_F <- chapman.mr(df = for_dobs, mismatch = FALSE)
chap_T <- chapman.mr(df = for_dobs, mismatch = TRUE)

# Compare estimates
chap_F
chap_T

nrow(for_dobs[for_dobs$obs==2 & for_dobs$detect==1, ])
nrow(for_dobs[for_dobs$obs==2, ])
nrow(for_dobs)

# PLOT MOVEMENT
# Visualize movement with truncation lines
movement_data <- for_dobs %>%
  filter(obs %in% c(1, 2)) %>%
  arrange(id, obs) %>%
  group_by(id) %>%
  summarize(
    x1 = x[obs == 1],
    y1 = y[obs == 1],
    x2 = x[obs == 2],
    y2 = y[obs == 2],
    detect1 = detect[obs == 1],
    detect2 = detect[obs == 2]
  ) %>%
  ungroup()

ggplot(data = movement_data,
       aes(x = x1, y = y1)) +
  geom_segment(aes(xend = x2, yend = y2),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "blue") +
  geom_point(aes(x = x1, y = y1), color = "green", size = 2) +
  geom_point(aes(x = x2, y = y2), color = "red", size = 2) +
  labs(title = "Movement from Observation 1 to 2",
       x = "X coordinate", y = "Y coordinate") +
  theme_minimal()

