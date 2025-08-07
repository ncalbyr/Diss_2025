library(circular)
library(gridExtra)
library(ggplot2)
library(grid)
# Build detections into moved data

move.detect <- function(move, sigma) {
  # 1. Use existing initial positions
  primate <- data.frame(x = primate.dat$x, y = primate.dat$y)
  
  # 2. Move animals using your existing movement code
  moved <- move.data(df = primate, move = move)
  
  # 3. Create long-format data
  df <- data.frame(
    id = rep(1:nrow(moved), each = 2),
    obs = rep(1:2, times = nrow(moved)),
    x = c(moved$x, moved$x2),
    y = c(moved$y, moved$y2),
    detect = NA
  )
  
  # 4. Simulate detection
  df$detect <- rbinom(nrow(moved) * 2, 1, prob = exp(-(df$x^2) / (2 * sigma^2)))
  
  # 5. Optional: remove undetected animals
  detected_ids <- unique(df$id[df$detect == 1])
  df <- subset(df, id %in% detected_ids)
  
  return(df)
}

set.seed(224)
detection_0 <- move.detect(move = 0, sigma = 0.02)
detection_1 <- move.detect(move = 1, sigma = 0.02)
detection_2 <- move.detect(move = 2, sigma = 0.02)
# View tables
head(detection_0)
head(detection_1)

# Save head of table as PNG
png("detection_0.png", width = 600, height = 250)
grid.table(head(detection_0, 6), rows = NULL)
dev.off()

# Estimate Abundance with Chapman Estimator
chapman.mr <- function(df, mismatch){
  if (mismatch==TRUE){df <- mismatch(df)}
  
  S1 <- nrow(df[df$obs==1 & df$detect==1, ])  # first occasion
  S2 <- nrow(df[df$obs==2 & df$detect==1, ])  # second occasion
  B <- df$detect[df$obs==1]==1 & df$detect[df$obs==2]==1  
  B <- length(B[B==TRUE])  # caught by both occasions
  N.hat <- (S1+1)*(S2+1)/(B+1)-1  # abundance estimate
  var.N <- (S1+1)*(S2+1)*(S1-B)*(S2-B)/(((B+1)^2)*(B+2))
  d <- exp(1.96*sqrt(log(1+(var.N/(N.hat^2)))))
  lcl <- N.hat/d; ucl <- N.hat*d
  return(c(N.hat, lcl, ucl))
}

first_chap_est <- chapman.mr(df = detection_0, mismatch = F)
first_chap_est/6 # Density

## The next step is to introduce a DUPLICATE PROCESS (IMPERFECT MATCHING)
#### This function already assumes PERFECT matching between occasions