# Load the packages, establish the GitHub connection
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(spatstat)
library(ggplot2)
library(truncnorm)
library(grid)
library(gridExtra)
library(circular)
##### 1) SIMULATE 210 ANIMALS
# (from perpendicular distribution of "primate.dat")
# This object uses pi.norm density distribution for pi.x
Fit.n.ip0$fit$par

hist(primate.dat$x)

# Parameters for "pi.x" distribution from "Fit.n.ip0"
beta1 <- Fit.n.ip0$fit$par[1]
beta2 <- Fit.n.ip0$fit$par[2]

lphi1 <- Fit.n.ip0$fit$par[3]
lphi2 <- Fit.n.ip0$fit$par[4]

# Generate x-values from the known perpendicular distribution
x_set <- rtruncnorm(n = 210, # number of items to generate
                    a = 0, # lower bound
                    b = 50, # upper bound
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

# Save head of table as PNG
png("detected_xy.png", width = 600, height = 250)
grid.table(head(detected_xy, 6), rows = NULL)
dev.off()

##### 3) PLUG X&Y's INTO LT2D MODEL FOR ESTIMATES
# Build data frame skeleton
n = dim(detected_xy)[1]
new.df = data.frame(stratum=rep(1,n), transect=rep(1,n), object=1:n, 
                        size=rep(1,n), area=rep(100,n), L=rep(100,n),
                        x=detected_xy$x, y=detected_xy$y)

# Try fitting:
Fit.new = LT2D.fit(new.df,b=c(beta1,beta2),hr="ip0",ystart=0.05,
                   pi.x="pi.norm",logphi=c(lphi1,lphi2),w=0.03,hessian=TRUE,
                   control=list(trace=5,maxit=1000))
Fit.new$fit$par

# Plot GoF and print Go0F ststistics
openGraph(h=4,w=9)
par(mfrow=c(1,2))
gof.ip0 = gof.LT2D(fit.n.ip0,plot=TRUE)
gof.ip0

# Plot fits
openGraph(h=4,w=9)
par(mfrow=c(1,2))
plot(Fit.new,smooth.fy=TRUE)

# Look at abundance estimates
Fit.new$ests
names(Fit.new)
names(Fit.new$fit)
Fit.new
Fit.new$fit$par
Fit.new$fit$pi.x
summary(Fit.new$fit)

##### 4) MOVE X&Y's TO X2&Y2's (Add them back into the fitted object?)

moving <- move.data(df = Fit.new$fit$dat,
                       move = 0,
                       keep_angle = F)

##### 5) USE DETECTION FUNCTION ON NEW LOCATIONS
# Test "ip0" probability
det.fun <- ip0(y = moving$y2,
            x = moving$x2,
            b = c(log(beta1),log(beta2)))
max(det.fun)
min(det.fun)
# Test rbinom()
rbinom(n = nrow(moving), size = 1,
       prob = ip0(y = moving$y2,
                  x = moving$x2,
                  b = c(beta1,beta2)))

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
  df$detect <- rbinom(nrow(df), 1, prob = exp(-(df$x^2) / (2 * sigma^2)))

  return(df)
}


detect.data(sigma = 0.2)

