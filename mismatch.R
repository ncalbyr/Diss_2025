library(fields)  # for rdist()
mismatch <- function(df){
  #browser()
  
  ## 1. get distance between every pair of detected objects
  df1 <- subset(df, df$obs==1 & df$detect==1)[ , c("x", "y")]  
  df2 <- subset(df, df$obs==2 & df$detect==1)[ , c("x", "y")]  
  dist.pair <- as.data.frame(rdist(df1, df2))
  
  dist.pair$unique <- 1:nrow(dist.pair)
  df1$id <- 1:nrow(df1); df1$detect <- df1$obs <- rep(1, nrow(df1))
  
  ## 2. using min distance to decide mismatching
  for (i in 1:(ncol(dist.pair)-1)) {
    min.index <- which.min(dist.pair[, i])
    if (length(min.index)>0){  
      if (dist.pair[min.index, i] < 0.03) {detect2 <- 1}
      else if (dist.pair[min.index, i] > 0.1) {detect2 <- 0}
      else{detect2 <- rbinom(1, 1, p.approx(ys <- seq(0, 0.35, length.out=100), dist.pair[min.index, i], ip0, b=c(6, 0.000005), what = "px"))}
    }else{detect2 <- 0}  # if no obs1 detection to match
    
    if (detect2==1){  # if matched
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], dist.pair$unique[min.index], 2, 1)
      dist.pair <- dist.pair[-min.index, ]
    }else{  # if no match
      id <- max(unique(df1$id))+1
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], id, 1, 0)
      df1[nrow(df1) + 1, ] <- c(df2$x[i], df2$y[i], id, 2, 1)}}
  
  for (i in dist.pair$unique){ 
    df1[nrow(df1) + 1, ] <- c(df1$x[df1$id==i], df1$y[df1$id==i], i, 2, 0)}
  
  ## 3. return new dataset
  df1 <- df1[order(df1$id), ]
  return(df1)
}

# Test
matched <- mismatch(detection)
head(matched)
# Save head of table as PNG
png("matched.png", width = 600, height = 250)
grid.table(head(matched, 6), rows = NULL)
dev.off()
