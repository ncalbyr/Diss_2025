library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')

##### LOAD PRIMATE DATA #####
data("primate.dat")
summary(primate.dat) # Gives x-y coordinates of primate observations
plot(primate.dat)

# Set bounds of survey area
w = 0.03 ; ystart = 0.35 # Perpendicular and forward truncation
L = 10 ; A = 2*w*L
# Build the data.frame() to fit to the 2D model
all.1s <- rep(1,length(primate.dat$x))
obj <- 1:length(primate.dat$x)
sim.df <- data.frame(x = primate.dat$x,
                     y = primate.dat$y,
                     stratum = all.1s,
                     transect = all.1s,
                     L = L,
                     area = A,
                     object = obj,
                     size = all.1s)
# Fit this data to LT2D.fit()
LT2D.fit(DataFrameInput = primate.dat,
         hr = )