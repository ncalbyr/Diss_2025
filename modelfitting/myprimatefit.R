library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')

#####################################
#####################################

data(primate.dat)

x=primate.dat$x

y=primate.dat$y



# Look at the data

openGraph(h=3,w=12)

par(mfrow=c(1,3))

pdlab="Perpendicular distance"

fdlab="Distance along transect"

plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="")

hist(y,breaks=seq(0,max(na.omit(y)),length=16),xlab=fdlab,main="")

hist(x,breaks=seq(0,max(na.omit(x)),length=12),xlab=pdlab,main="")



# Fit the model

# Normal bump with ip0 hazard function (SELECTED MODEL in paper):
ystart = 0.35
b=c(5, 8)

logphi=c(0.02, -4.4)

w=0.03 # this is what is used in the paper

fit.n.ip0=fityx(y[x<=w],x[x<=w],b=b,hr="ip0",ystart=ystart,
                
                pi.x="pi.norm",logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))



# Looks like there's a bug in fityx - it should have added $covariates. Add manually:

fit.n.ip0$covariates = FALSE



# Plot fits

openGraph(h=4,w=9)

par(mfrow=c(1,2))

plot(fit.n.ip0,smooth.fy=TRUE)



# Plot GoF and print GoF statistics

gof.ip0 = gof.LT2D(fit.n.ip0,plot=TRUE)

gof.ip0

# Make GoF .png
# Sample GoF data
gof.ip0 <- list(
  X = c(`K-S` = 0.5429783, `CvM` = 0.6859833),
  Y = c(`K-S` = 0.7491595, `CvM` = 0.7413269)
)

# Convert to data frame for display
gof_df <- data.frame(
  Statistic = names(gof.ip0$X),
  X = as.numeric(gof.ip0$X),
  Y = as.numeric(gof.ip0$Y)
)

# Load gridExtra for table rendering
library(gridExtra)
library(grid)

# Save as PNG
png("gof.png", width = 600, height = 300)
grid.table(gof_df)
dev.off()
