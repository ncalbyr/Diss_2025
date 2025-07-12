# Install packages and establish GitHub connection
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(circular)

# Load primate data
data(primate.dat)
primate <- primate.dat
head(primate)

# Movement function/second occasion
move.data <- function(df, move = 0) {
  # df must have columns x and y
  n <- nrow(df)
  
  if (move == 0) {
    # Avoidance movement (away from x = 0)
    nleft <- sum(df$x <= 0)
    
    df$angle <- NA
    angleleft <- rwrappedcauchy(nleft, mu = circular(pi), rho = 0.9)
    df$angle[df$x <= 0] <- as.numeric(angleleft)
    
    angleright <- rwrappedcauchy(n - nleft, mu = circular(0), rho = 0.9)
    df$angle[df$x > 0] <- as.numeric(angleright)
    
    # Distance increases with distance from edge
    dist <- (max(abs(df$x)) - abs(df$x)) * runif(n, 0, 0.1) * 2.5
    
  } else if (move == 1) {
    # Random movement
    df$angle <- as.numeric(rwrappedcauchy(n, mu = circular(0), rho = 0))
    dist <- rlnorm(n, meanlog = -3, sdlog = 0.5)
    
  } else {
    stop("Invalid move value. Use 0 (avoidance) or 1 (random).")
  }
  
  # Compute new x-y locations
  df$x2 <- df$x + dist * cos(df$angle)
  df$y2 <- df$y + dist * sin(df$angle) / 1000  # small forward shift
  
  df$angle <- NULL
  return(df)
}

# Test function
set.seed(124)
primate_moved<- move.data(df = primate,
                                move = 1)
head(primate_moved)

# Visualize movement with truncation lines
plot(primate$x, primate$y, pch = 16, col = "blue", 
     xlim = range(primate_moved$x, primate_moved$x2),
     ylim = range(primate_moved$y, primate_moved$y2), 
     xlab = "Perpendicular (x)", ylab = "Forward (y)")

# Add movement arrows
arrows(primate$x, primate$y, primate_moved$x2, primate_moved$y2, 
       length = 0.05, col = "darkgray")

# Add second location points
points(primate_moved$x2, primate_moved$y2, pch = 17, col = "red")

# Add truncation lines
abline(v = c(-0.03, 0.03), col = "black", lty = 2)  # vertical truncation
abline(h = 0.15, col = "black", lty = 2)           # horizontal truncation

# Add legend
legend("topright", legend = c("Occasion 1", "Occasion 2"), 
       pch = c(16, 17), col = c("blue", "red"))
