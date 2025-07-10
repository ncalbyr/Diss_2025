library(devtools)

install_github("david-borchers/LT2Dcal",force=TRUE)

library('LT2D')

# set perp trunc, forward trunc
# total transect length and survey
# area:
w = 0.15 ; ystart = 0.55
L = 10 ; A = 2*w*L

# set value of 'true' parameters for simulated data:

# set up and plot detection function
nx = 50
ny=100
xs = seq(0,w,length=nx)
ys = seq(0,ystart,length=ny)
h.fun.name = "h1"
b=c(-3, 0.8)
h.fun = match.fun(h.fun.name) # make h.fun the function specified via a character variable
par(mfrow=c(1,3))
# plot perp dist detection function
p.vals = p.approx(ys,xs,h.fun,b) # detection function values to plot
plot(xs,p.vals,type='l',ylim=range(0,p.vals),xlab='Perp. distance, x',ylab=expression(p(x)))
# and now the 2D detection function:
pmat = p.approx(ys,xs,h.fun,b,what="pxy")
persp(x=xs,y=ys,z=pmat,theta=120,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="p(det by y)")
# and now the 2D pdf:
pdfmat = p.approx(ys,xs,h.fun,b,what="fxy")
persp(x=xs,y=ys,z=pdfmat,theta=45,phi=25, xlab="Perp dist (x)", ylab="Forward dist (y)",
      zlab="f(y|x)")

# Set up and plot density function 
pi.fun.name = "pi.chnorm" # specify density function name
logphi <- c(-0.01,log(0.03)) # density function parameters
pi.fun = match.fun(pi.fun.name) # make Dfun the function specified via a character variable
D.pdf = pi.fun(xs,logphi,w) # density pdf values to plot
plot(xs,D.pdf,type='l',ylim=range(0,D.pdf),xlab='Perp. distance, x',ylab=expression(pi(x)))

# produce simulated data:
set.seed(3)
simDat = simXY(50, pi.fun.name,
               logphi, 'h1', 
               b, w, 
               ystart)

# Look at data:
plot(simDat)

# create the data.frame:
all.1s <- rep(1,length(simDat$locs$x))
obj <- 1:length(simDat$locs$x)
sim.df <- data.frame(x = simDat$locs$x,
                     y = simDat$locs$y,
                     stratum = all.1s,
                     transect = all.1s,
                     L = L,
                     area = A,
                     object = obj,
                     size = all.1s)

# fit an LT2D model
fit <- LT2D.fit(DataFrameInput = sim.df,
                hr = 'h1',
                # start values for b:
                b = b,
                ystart = ystart,
                pi.x = 'pi.norm',
                # start values for logphi:
                logphi = logphi,
                w = w,
                hessian = TRUE,
                control=list(trace=5))

# plot fitted functions
par(mfrow=c(1,2))
plot(fit)

# print table of point estimates
fit$ests

# look at goodness of fit
par(mfrow=c(1,2))
gof.LT2D(fit, plot=TRUE)

# Look at AIC
fit$fit$AIC

# bootstrap for variance and interval estimates of total abundance
boot <- LT2D.bootstrap(fit,r=999,alpha = 0.05)
boot$ci
# plot bootstrap histogram:
hist(boot$Ns,main='',xlab='Estimate of Abundance')
# add point estimate and CI
mle.N = fit$ests[nrow(fit$ests),ncol(fit$ests)]
points(mle.N,0,pch=19,col="red")
arrows(boot$ci[1],0,boot$ci[2],0,angle=90,code=3,col="red")
