##### 6) LOOP SIMULATION CODE TO PRODUCE DISTRIBUTION OF ESTIMATES

####### MAKE SURE PARAMETERS FOR BETA AND LPHI FROM "primatefit.R"
######### ARE USING THE SAMING SETTINGS (e.g. perpendicular truncation "w")
# Install packages and establish GitHub connection
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
# LT2D Method (Loop)
library(truncnorm)

# Initialize storage
x_set <- list()
detected_xy <- list()
lt2d_density <- numeric(100)  # However many estimates you want (at least 100)
                            # Even 50 took a long time (6 hours)

### Use parameters from "primatefit.R"
beta1 <- Fit.n.ip0$fit$par[1]
beta2 <- Fit.n.ip0$fit$par[2]

lphi1 <- Fit.n.ip0$fit$par[3]
lphi2 <- Fit.n.ip0$fit$par[4]
###

for(i in 1:100){
  # Generate x-values from the known perpendicular distribution
  x_set[[i]] <- rtruncnorm(n = 210, 
                           a = 0, 
                           b = 0.05, 
                           mean = lphi1, sd = exp(lphi2))
  
  # Generate forward distances (y) ONLY for detected individuals
  detected_xy[[i]] <- detect2DLT(x = x_set[[i]], hr = "ip0", b = c(beta1,beta2),
                                 ystart = 0.05, ny = 1000,
                                 getIDs = TRUE)
  
  # Build data frame for the fitting function
  n <- nrow(detected_xy[[i]])
  new.df <- data.frame(stratum = rep(1, n), 
                       transect = rep(1, n), 
                       object = 1:n, 
                       size = rep(1, n), 
                       area = rep(100, n), 
                       L = rep(100, n),
                       x = detected_xy[[i]]$x, 
                       y = detected_xy[[i]]$y)
  
  # Fit model
  Fit.new <- LT2D.fit(new.df, 
                      b = c(beta1, beta2), 
                      hr = "ip0", 
                      ystart = 0.05,
                      pi.x = "pi.norm", 
                      logphi = c(lphi1, lphi2), 
                      w = 0.04, # test 0.03 and 0.04
                      hessian = TRUE,
                      control = list(trace = 5, maxit = 1000))
  
  # Store Density estimate
  lt2d_density[i] <- Fit.new$ests$D[1]
}

# Show estimates
Fit.new$ests
png(filename = "lt2d_loop_estimates.png", width = 700, height = 100)
grid.table(Fit.new$ests)
dev.off()
## Build plots
par(mfrow=c(1,2))

# Density
mean_density_lt2d <- mean(lt2d_density)
hist(lt2d_density, main = "LT2D Density (w=0.04)",xlab="Primate Density")
abline(v = 210/6, col = "red", lwd = 2)
abline(v = mean_density_lt2d, col = "blue", lwd = 2, lty = 5)
legend(x = "topright", legend = c("True Density","Esimated Density"),col = c("red","blue"),
       lty = c(5,1))

# Abundance
lt2d_abundance <- lt2d_density*6
mean_abund_lt2d <- mean(lt2d_abundance)
hist(lt2d_abundance, main = "LT2D Abundance (w = 0.04)",xlab="Primate Abundance")
abline(v = 210, col = "red", lwd = 2)
abline(v = mean_abund_lt2d, col = "blue", lwd = 2, lty = 5)
legend(x = "topright", legend = c("True Abundance","Esimated Abundance"),col = c("red","blue"),
       lty = c(5,1))
