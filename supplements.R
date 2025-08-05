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
hist(inverse_hazard)
