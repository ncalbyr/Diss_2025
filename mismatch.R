library(fields)  # for rdist()
library(tidyverse)
sim.mismatch.long <- function(df) {
  # Ensure we only consider individuals seen at least once
  seen_once <- df %>%
    group_by(id) %>%
    summarise(n_detected = sum(detect)) %>%
    filter(n_detected > 0) %>%
    pull(id)
  
  df_filtered <- df %>% filter(id %in% seen_once)
  
  # Spread into detection history per individual
  detect_by_obs <- df_filtered %>%
    select(id, obs, detect) %>%
    tidyr::pivot_wider(names_from = obs, values_from = detect, names_prefix = "obs") %>%
    mutate(
      mismatch = ifelse(obs1 != obs2, 1, 0)
    )
  
  # Return summary
  mismatch_rate <- mean(detect_by_obs$mismatch, na.rm = TRUE)
  
  return(list(
    mismatch_df = detect_by_obs,
    mismatch_rate = mismatch_rate
  ))
}
# Test
sim.mismatch.long(detection)

sim.mismatch <- function(df) {
  df1 <- subset(df, obs == 1 & detect == 1)[, c("x", "y")]
  df2 <- subset(df, obs == 2 & detect == 1)[, c("x", "y")]
  
  dist.pair <- as.data.frame(rdist(df1, df2))  # n1 x n2
  dist.pair$unique <- 1:nrow(dist.pair)
  
  df1$id <- 1:nrow(df1)
  df1$detect <- 1
  df1$obs <- 1
  
  # Set realistic distance thresholds (assuming km units)
  min_thresh <- 0.03   # ~30 meters
  max_thresh <- 0.10   # ~100 meters
  
  for (i in 1:(ncol(dist.pair) - 1)) {
    min.index <- which.min(dist.pair[, i])
    
    if (length(min.index) > 0) {
      d <- dist.pair[min.index, i]
      if (d < min_thresh) {
        detect2 <- 1
      } else if (d > max_thresh) {
        detect2 <- 0
      } else {
        # Placeholder detection probability (or use p.approx if you define it)
        detect2 <- rbinom(1, 1, prob = exp(-(d^2) / (2 * 0.02^2)))
      }
    } else {
      detect2 <- 0
    }
    
    if (detect2 == 1) {
      # Matched: add as occasion 2 for same ID
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], df1$id[min.index], 2, 1)
      dist.pair <- dist.pair[-min.index, ]
    } else {
      # Not matched: create new ID for this occasion 2 detection
      new_id <- max(df1$id) + 1
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], new_id, 1, 0)
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], new_id, 2, 1)
    }
  }
  
  # Fill in unmatched occasion 1 entries at occasion 2 as not detected
  for (i in dist.pair$unique) {
    df1[nrow(df1) + 1, ] <- c(df1$x[df1$id == i], df1$y[df1$id == i], i, 2, 0)
  }
  
  df1 <- df1[order(df1$id), ]
  rownames(df1) <- NULL
  return(df1)
}

# Test function

sim.mismatch(detection)
