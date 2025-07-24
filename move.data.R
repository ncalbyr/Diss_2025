# Install packages and establish GitHub connection
library(devtools)
install_github("david-borchers/LT2Dcal",force=TRUE)
library('LT2D')
library(circular)
library(gridExtra)
library(ggplot2)
library(grid)

# Load primate data
data(primate.dat)
primate <- primate.dat
head(primate)

# Movement function/second occasion
move.data <- function(df, move = 0, keep_angle = FALSE) {
  # df must have columns x and y (for occasion 1 only)
  n <- nrow(df)
  
  # Initialize angle and distance
  df$angle <- NA
  dist <- rep(NA, n)
  
  if (move == 0) {
    # Avoidance movement
    nleft <- sum(df$x <= 0)
    
    # Angles away from transect line
    angleleft <- rwrappedcauchy(nleft, mu = circular(pi), rho = 0.9)
    df$angle[df$x <= 0] <- as.numeric(angleleft)
    
    angleright <- rwrappedcauchy(n - nleft, mu = circular(0), rho = 0.9)
    df$angle[df$x > 0] <- as.numeric(angleright)
    
    # Distance depends on distance from center
    dist <- (max(abs(df$x)) - abs(df$x)) * runif(n, 0, 0.1) * 2.5
    
  } else if (move == 1) {
    # Random movement
    angle1 <- rwrappedcauchy(n, mu = circular(0), rho = 0)
    df$angle <- as.numeric(angle1)
    
    dist <- rlnorm(n, log(0.01), 0.001) # (n, meanlog, sdlog)
    
  } else if (move == 2) {
    # no movement
    df$angle <- 0
    dist <- 0
  } else {
    stop("Invalid move value. Use 0 (avoidance) or 1 (random).")
  }
  
  # Compute new x and y positions for second occasion
  df$x2 <- df$x + dist * cos(df$angle)
  df$y2 <- df$y + dist * sin(df$angle)  # scale forward movement
  
  if (!keep_angle) df$angle <- NULL
  return(df)
}


# Test function
set.seed(126)
primate_moved_2<- move.data(df = primate,
                            move = 2)

head(primate_moved_2)


# Visualize movement with truncation lines
plot(primate$x, primate$y, pch = 16, col = "blue", 
     xlim = range(primate_moved_0$x, primate_moved_0$x2),
     ylim = range(primate_moved_0$y, primate_moved_0$y2), 
     xlab = "Perpendicular (x)", ylab = "Forward (y)")

# Add movement arrows
arrows(primate$x, primate$y, primate_moved_0$x2, primate_moved_0$y2, 
       length = 0.05, col = "darkgray")

# Add second location points
points(primate_moved_0$x2, primate_moved_0$y2, pch = 17, col = "red")

# Add truncation lines
abline(v = c(-0.03, 0.03), col = "black", lty = 2)  # vertical truncation
abline(h = 0.15, col = "black", lty = 2)           # horizontal truncation

# Add legend
legend("topright", legend = c("Occasion 1", "Occasion 2"), 
       pch = c(16, 17), col = c("blue", "red"))


# Other movement option
set.seed(124)
primate_moved_1<- move.data(df = primate,
                                move = 1)
head(primate_moved_1)

# Visualize movement with truncation lines
plot(primate$x, primate$y, pch = 16, col = "blue", 
     xlim = range(primate_moved_1$x, primate_moved_1$x2),
     ylim = range(primate_moved_1$y, primate_moved_1$y2), 
     xlab = "Perpendicular (x)", ylab = "Forward (y)")

# Add movement arrows
arrows(primate$x, primate$y, primate_moved_1$x2, primate_moved_1$y2, 
       length = 0.05, col = "darkgray")

# Add second location points
points(primate_moved_1$x2, primate_moved_1$y2, pch = 17, col = "red")

# Add truncation lines
abline(v = c(-0.03, 0.03), col = "black", lty = 2)  # vertical truncation
abline(h = 0.15, col = "black", lty = 2)           # horizontal truncation

# Add legend
legend("topright", legend = c("Occasion 1", "Occasion 2"), 
       pch = c(16, 17), col = c("blue", "red"))



### PNG Generation
# Save head of table as PNG
png("primate_moved.png", width = 600, height = 250)
grid.table(head(primate_moved_0, 6), rows = NULL)
dev.off()
