# Build detections into moved data

detect.data <- function(move, sigma) {
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

detect.data(move = 1, sigma = 0.02)
