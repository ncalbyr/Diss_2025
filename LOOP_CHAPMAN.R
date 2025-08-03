################################################################################
# Make sure to run "mismatch.R" before getting here, if you are testing mismatch
################################################################################
# Load packages
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(truncnorm)
library(circular)


# Build function
simulate_chapman <- function(Fit.n.ip0, n_animals = 210, area = 100, 
                            beta = NULL, lphi = NULL, mismatch = FALSE,
                            move = 0) {
  # 1. Extract parameters
  if (is.null(beta)) {
    beta1 <- Fit.n.ip0$fit$par[1]
    beta2 <- Fit.n.ip0$fit$par[2]
  } else {
    beta1 <- beta[1]
    beta2 <- beta[2]
  }
  if (is.null(lphi)) {
    lphi1 <- Fit.n.ip0$fit$par[3]
    lphi2 <- Fit.n.ip0$fit$par[4]
  } else {
    lphi1 <- lphi[1]
    lphi2 <- lphi[2]
  }
  
  # 2. Simulate x-values (perpendicular distances)
  x_set <- rtruncnorm(n = n_animals, a = 0, b = 0.05, mean = lphi1, sd = exp(lphi2))
  
  # 3. Simulate forward distances for detections
  detected_xy <- detect2DLT(x = x_set, hr = "ip0", b = c(beta1, beta2),
                            ystart = 0.05, ny = 1000, getIDs = TRUE)
  if (nrow(detected_xy) < 5) return(NA)
  
  # 4. Move to 2nd observer location
  moving <- move.data(df = detected_xy, move = move, keep_angle = FALSE)

  # 5. Simulate detection outcomes for observer 2
  df <- data.frame(
    id = rep(1:nrow(moving), each = 2),
    obs = rep(1:2, times = nrow(moving)),
    x = c(moving$x, moving$x2),
    y = c(moving$y, moving$y2),
    detect = NA
  )
  
  # Use approximation of detection probability for observer 2
  ys <- seq(0, 0.05, length.out = 100)
  obs2.probs <- p.approx(ys, df$x[df$obs == 2], ip0, b = c(beta1, beta2), what = "px")
  
  n_obs <- nrow(df) / 2
  df$detect[df$obs == 2] <- rbinom(n_obs, 1, obs2.probs)
  df$detect[df$obs == 1] <- 1  
  df$detect[abs(df$x) > 0.05] <- 0
  
  # 6. Apply Chapman estimator
  result <- tryCatch({
    chapman.mr(df = df, mismatch = mismatch)
  }, error = function(e) return(NA))
  
  return(result)
}

###############################################################################

set.seed(42)
n_simulations <- 100
#########################################################
# NO mismatch/NO movement
chapman_nm_none <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = FALSE, move = 2),
  simplify = FALSE
  )
#########################################################
# NO mismatch/AVOIDANT movement
chapman_nm_avoid <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = FALSE, move = 0),
  simplify = FALSE)
#########################################################
# NO mismatch/RANDOM movement
chapman_nm_random <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = FALSE, move = 1),
  simplify = FALSE)
#########################################################
# YES mismatch/AVOIDANT movement
chapman_ym_avoid <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = TRUE, move = 0),
  simplify = FALSE)
#########################################################
# YES mismatch/RANDOM movement
chapman_ym_random <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = TRUE, move = 1),
  simplify = FALSE)
#########################################################
# YES mismatch/NO movement
chapman_ym_none <- replicate(
  n_simulations,
  simulate_chapman(Fit.n.ip0, mismatch = TRUE, move = 2),
  simplify = FALSE)
#########################################################
# Clean and convert to data.frames
chapman_df_nm_none <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_nm_none)))
chapman_df_nm_none$density <- chapman_df_nm_none$V1/6
chapman_df_nm_avoid <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_nm_avoid)))
chapman_df_nm_avoid$density <- chapman_df_nm_avoid$V1/6
chapman_df_nm_random <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_nm_random)))
chapman_df_nm_random$density <- chapman_df_nm_random$V1/6

chapman_df_ym_none <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_ym_none)))
chapman_df_ym_none$density <- chapman_df_ym_none$V1/6
chapman_df_ym_avoid <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_ym_avoid)))
chapman_df_ym_avoid$density <- chapman_df_ym_avoid$V1/6
chapman_df_ym_random <- as.data.frame(do.call(rbind, Filter(Negate(is.na), chapman_ym_random)))
chapman_df_ym_random$density <- chapman_df_ym_random$V1/6
#########################################################
# Name values
colnames(chapman_df_nm_none) <- c("Nhat", "LCL", "UCL","density")
colnames(chapman_df_nm_avoid) <- c("Nhat", "LCL", "UCL","density")
colnames(chapman_df_nm_random) <- c("Nhat", "LCL", "UCL","density")

colnames(chapman_df_ym_none) <- c("Nhat", "LCL", "UCL","density")
colnames(chapman_df_ym_avoid) <- c("Nhat", "LCL", "UCL","density")
colnames(chapman_df_ym_random) <- c("Nhat", "LCL", "UCL","density")

#########################################################
## Compute ABUNDANCE means
#no mismatch
mean_abund_nm_none <- mean(chapman_df_nm_none$Nhat)
mean_abund_nm_avoid <- mean(chapman_df_nm_avoid$Nhat)
mean_abund_nm_random <- mean(chapman_df_nm_random$Nhat)

#yes mismatch
mean_abund_ym_none <- mean(chapman_df_ym_none$Nhat)
mean_abund_ym_avoid <- mean(chapman_df_ym_avoid$Nhat)
mean_abund_ym_random <- mean(chapman_df_ym_random$Nhat)

## Compute DENSITY means
#no mismatch
mean_density_nm_none <- mean(chapman_df_nm_none$density)
mean_density_nm_avoid <- mean(chapman_df_nm_avoid$density)
mean_density_nm_random <- mean(chapman_df_nm_random$density)
#yes mismatch
mean_density_ym_none <- mean(chapman_df_ym_none$density)
mean_density_ym_avoid <- mean(chapman_df_ym_avoid$density)
mean_density_ym_random <- mean(chapman_df_ym_random$density)

# Open plotting window
openGraph(h = 6, w = 10)
par(mfrow = c(1, 3))
###############################################################################
###########################   Abundance   #####################################
###############################################################################
## --- Abundance Estimates (No Mismatch, No Movement) ---
hist(chapman_df_nm_none$Nhat, breaks = 20,
     main = "Chapman Estimator (No mismatch, No movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)  # True value
abline(v = mean_abund_nm_none, col = "blue", lwd = 2, lty = 2)  # Mean estimate
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)

## --- Abundance Estimates (No Mismatch, Avoidant Movement) ---
hist(chapman_df_nm_avoid$Nhat, breaks = 20,
     main = "Chapman Estimator (No mismatch, Avoidant movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)  # True value
abline(v = mean_abund_nm_avoid, col = "blue", lwd = 2, lty = 2)  # Mean estimate
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)

## --- Abundance Estimates (No Mismatch, Random Movement) ---
hist(chapman_df_nm_random$Nhat, breaks = 20,
     main = "Chapman Estimator (No mismatch, Random movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)  # True value
abline(v = mean_abund_nm_random, col = "blue", lwd = 2, lty = 2)  # Mean estimate
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)

###############################################################################
## --- Abundance Estimates (Yes mismatch, No movement) ---
hist(chapman_df_ym_none$Nhat, breaks = 20,
     main = "Chapman Estimator (Yes mismatch, No movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)
abline(v = mean_abund_ym_none, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Abundance Estimates (Yes mismatch, Avoidant movement) ---
hist(chapman_df_ym_avoid$Nhat, breaks = 20,
     main = "Chapman Estimator (Yes mismatch, Avoid movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)
abline(v = mean_abund_ym_avoid, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Abundance Estimates (Yes mismatch, Random movement) ---
hist(chapman_df_ym_random$Nhat, breaks = 20,
     main = "Chapman Estimator (Yes mismatch, Random movement)",
     xlab = "Abundance Estimate", col = "lightgray", border = "white")
abline(v = 210, col = "red", lwd = 2)
abline(v = mean_abund_ym_random, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Abundance", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
###############################################################################
##########################   DENSITY ##########################################
###############################################################################
## --- Density Estimates (No mismatch, No movement) ---
hist(chapman_df_nm_none$density, breaks = 20,
     main = "Density Estimate (No mismatch, No movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_nm_none, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Density Estimates (No mismatch, Avoidant movement) ---
hist(chapman_df_nm_avoid$density, breaks = 20,
     main = "Density Estimate (No mismatch, Avoidant movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_nm_avoid, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Density Estimates (No mismatch, Random movement) ---
hist(chapman_df_nm_random$density, breaks = 20,
     main = "Density Estimate (No mismatch, Random movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_nm_random, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
###############################################################################
## --- Density Estimates (Yes mismatch, No movement) ---
hist(chapman_df_ym_none$density, breaks = 20,
     main = "Density Estimate (Yes mismatch, No movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_ym_none, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Density Estimates (Yes mismatch, Avoidant movement) ---
hist(chapman_df_ym_avoid$density, breaks = 20,
     main = "Density Estimate (Yes mismatch, Avoidant movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_ym_avoid, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
## --- Density Estimates (Yes mismatch, Random movement) ---
hist(chapman_df_ym_random$density, breaks = 20,
     main = "Density Estimate (Yes mismatch, Random movement)",
     xlab = "Density Estimate", col = "lightgray", border = "white")
abline(v = 35, col = "red", lwd = 2)
abline(v = mean_density_ym_random, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("True Density", "Mean Estimate"),
       col = c("red", "blue"), lty = c(1, 2), lwd = 2)
