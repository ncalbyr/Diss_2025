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
ystart = 0.05 # forward distance beyond which it is IMPOSSIBLE to observe
b=c(5, 8)
logphi=c(0.02, -4.4)
w=0.03 # this is what is used in the paper
fit.n.ip0=fityx(y[x<=w],x[x<=w],b=b,hr="ip0",ystart=ystart,
                pi.x="pi.norm",logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))

# Looks like there's a bug in fitxy - it should have added $covariates. Add manually:
fit.n.ip0$covariates = FALSE

# Plot fits
openGraph(h=4,w=9)
par(mfrow=c(1,2))
plot(fit.n.ip0,smooth.fy=TRUE)

# Plot GoF and print Go0F ststistics
gof.ip0 = gof.LT2D(fit.n.ip0,plot=TRUE)
gof.ip0

fit.n.ip0$counts


# -------------------------------------------------
# Try fitting with LT2D.fit
# -------------------------------------------------

# First need to make a Distance-like data frame:
# (Put arbitrary numbers in all columns except x, y)
n = dim(primate.dat)[1]
primate.df = data.frame(stratum=rep(1,n), transect=rep(1,n), object=1:n, 
                        size=rep(1,n), area=rep(100,n), L=rep(100,n),
                        x=primate.dat$x, y=primate.dat$y)

# Try fitting:
Fit.n.ip0=LT2D.fit(primate.df,b=b,hr="ip0",ystart=ystart,
                   pi.x="pi.norm",logphi=logphi,w=w,hessian=TRUE,
                   control=list(trace=5,maxit=1000))

## Compare fit to that from fitxy:
# Parameter estimate
fit.n.ip0$par
Fit.n.ip0$fit$par
# Coefficient of variation
fit.n.ip0$CVpar
Fit.n.ip0$fit$CVpar

# From which we conclude that the fits are identical

# Plot fits
openGraph(h=4,w=9)
par(mfrow=c(1,2))
plot(Fit.n.ip0,smooth.fy=TRUE)

# Plot GoF and print Go0F ststistics
openGraph(h=4,w=9)
par(mfrow=c(1,2))
gof.ip0 = gof.LT2D(fit.n.ip0,plot=TRUE)
gof.ip0

# Look at abundance estimates
Fit.n.ip0$ests
names(Fit.n.ip0)
names(Fit.n.ip0$fit)
Fit.n.ip0$fit$par
Fit.n.ip0$fit$pi.x
summary(Fit.n.ip0$fit)
