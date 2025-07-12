# Install packages and establish GitHub connection
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(circular)

# Build function
sim.data.real <- function(primate, move) {
  n <- nrow(primate)
  
  df <- data.frame(id = rep(1:n, each = 2),
                   obs = rep(1:2, n),
                   x = NA, y = NA, 
                   forw.dist = NA, detect = NA, 
                   angle = NA)
  
  # Assign observed positions from primate data
  df$x[df$obs == 1] <- primate$x
  df$y[df$obs == 1] <- primate$y
  
  # First observer detection
  detect.2d <- detect2DLT(primate$x, hr = ip0, b = c(4.9, 0.036), ystart = 1700, ny = 1000)
  df$detect[df$obs == 1] <- as.numeric(primate$x %in% detect.2d$x)
  df$forw.dist[df$obs == 1 & df$detect == 1] <- detect.2d$y
  
  # Movement model
  if (move == 0) { # avoidance movement
    nleft <- sum(df$x[df$obs == 1] <= 0)
    
    angleleft <- rwrappedcauchy(nleft, mu = circular(pi), rho = 0.9)
    df$angle[df$x <= 0 & df$obs == 1] <- as.numeric(angleleft)
    
    angleright <- rwrappedcauchy(n - nleft, mu = circular(0), rho = 0.9)
    df$angle[df$x > 0 & df$obs == 1] <- as.numeric(angleright)
    
    dist <- (max(abs(df$x[df$obs == 1])) - abs(df$x[df$obs == 1])) * runif(n, 0, 0.1) * 2.5
  } else if (move == 1) { # random movement
    df$angle <- rwrappedcauchy(n, mu = circular(0), rho = 0)
    dist <- rlnorm(n, log(12) + 1.5, 1.5)
  }
  
  # Apply movement to second observation
  df$x[df$obs == 2] <- df$x[df$obs == 1] + dist * cos(df$angle[df$obs == 1])
  df$y[df$obs == 2] <- df$y[df$obs == 1] + dist * sin(df$angle[df$obs == 1]) / 1000
  
  # Second observer detection
  ys <- seq(0, 1700, length.out = 100)
  obs2.probs <- p.approx(ys, df$x[df$obs == 2], ip0, b = c(4.9, 0.036), what = "px")
  df$detect[df$obs == 2] <- rbinom(n, 1, obs2.probs)
  
  df$detect[abs(df$x) > 1600] <- 0
  
  # Filter individuals not detected at all
  df$keep <- rep(!(df$detect[df$obs == 1] == 0 & df$detect[df$obs == 2] == 0), each = 2)
  df <- subset(df, keep == TRUE)
  df[, c("angle", "keep")] <- NULL
  
  return(df)
}

# Test function
sim.data.real(primate, move = 0)

# Visualize
library(ggplot2)

# Assuming df is the output of sim.data.real()
df <- sim.data.real(primate, move = 1)

# Reshape for path plotting
library(tidyr)
library(dplyr)

df_path <- df %>%
  pivot_wider(names_from = obs, values_from = c(x, y, detect), names_prefix = "obs") %>%
  rename(x1 = x_obs1, y1 = y_obs1, detect1 = detect_obs1,
         x2 = x_obs2, y2 = y_obs2, detect2 = detect_obs2)
df_path$detect2 <- factor(df_path$detect2, levels = c(0, 1), labels = c("Not Detected", "Detected"))

ggplot(df_path) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = detect2),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.6) +
  geom_point(aes(x = x1, y = y1), color = "blue", size = 2, alpha = 0.6) +
  geom_point(aes(x = x2, y = y2), color = "red", size = 2, alpha = 0.6) +
  scale_color_manual(values = c("Not Detected" = "gray", "Detected" = "green")) +
  theme_minimal() +
  labs(title = "Animal Movement and Detection",
       x = "X Coordinate", y = "Y Coordinate",
       color = "Detection at Obs 2")
