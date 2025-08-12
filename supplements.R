# Example distributions or plots
# This file is for any extra supplementary figures or calculations
## used to further enhance the report
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
help(ip0)
inverse_hazard <- ip0(y = primate.dat$y,
                  x = primate.dat$x,
                  b = c(beta1,beta2))
curve(expr = ip0(y = primate.dat$y,
                 x = primate.dat$x,
                 b = c(beta1,beta2)),
      n = 100)
plot(inverse_hazard)
curve(inverse_hazard)
lines(inverse_hazard)
hist(inverse_hazard)
##### CODE CHUNK FROM "david-borchers/LT2Dcal/inst/2DLTfunctions mixture_simpson.R" #####
b_test=c(beta1, beta2) # set parameters
#hh=ip0(y=primate.dat$y, x=primate.dat$x, b=b_test) # pick hazard function
yy=seq(0,0.05,length=100);xx=rep(0,100) # set truncation according to how models were fit
hh=ip0(yy,xx,b=b)
plot(yy,hh,type="l", xlab = "Distance (km)",main = "Inverse Power Hazard Detection Function")
abline(v = 0.005, col = "red", lty = 5, lwd = 2)
abline(v = 0.01, col= "red", lty = 5,lwd = 2)
abline(v = 0.03, col = "blue", lwd = 2)
###############################################################################
# Test pi.norm (from LT2D package)
plot(seq(0,0.03,length=100),
     pi.norm(x = seq(0,0.03,length=100),
             logphi = c(0.02, -4.4), w = 0.03),
     xlab = "Perpendicular Distance (km)",
     ylab = "Frequency",
     main = "Perpendicular Distance Distribution (pi(x))"
     )
# So it doesn't need to use "pi.chnorm", the parameters do the work
# "logphi" parameters represent MEAN and DISPERSION (standard deviation)!!!!
# They are DESCRIBING the "logphi" distribution which is a PDF
library(circular)
plot(circular(0))

sum(primate$x <= 0) # There are no negative distances here
sum(primate$x > 0)
plot(rwrappedcauchy(n = sum(primate$x > 0),
               mu = circular(0),
               rho = 0.9))
moveright <- rwrappedcauchy(n = sum(primate$x > 0),
               mu = circular(0),
               rho = 0.9)

primate$angle[primate$x > 0] <- as.numeric(moveright)
# Deciding distance
plot((max(abs(primate$x)) - abs(primate$x)) * runif(n, 0, 0.1) * 2.5
)
