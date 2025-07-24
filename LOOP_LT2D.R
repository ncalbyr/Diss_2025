##### 6) LOOP SIMULATION CODE TO PRODUCE DISTRIBUTION OF ESTIMATES

# LT2D Method (Loop)
library(truncnorm)

# Initialize storage
x_set <- list()
detected_xy <- list()
abundance <- numeric(100)  # However many estimates you want (at least 100)
                            # Even 50 took a long time (6 hours)
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
                      w = 0.03, 
                      hessian = TRUE,
                      control = list(trace = 5, maxit = 1000))
  
  # Store Density estimate
  density[i] <- Fit.new$ests$D[1]
}
lt2d_density <- abundance # This took approximately 6 hours to produce
                      ## All maxit=1000, numeric(50),1:50
hist(lt2d_density)
abline(v=35,col="red")
par(mfrow=c(1,1))
lt2d_abundance <- lt2d_density*6
hist(lt2d_abundance)
