fit.mrds <- function(df, mismatch){
  if(mismatch){
    df <- sim.mismatch(df)
    names(df) <- c( "x","y","object","observer","detected")
  }else{ #change names to be in format for ddf
    names(df) <- c("object","observer", "x","y","forw.dist","detected")
  }
  df$distance <- abs(df$x)
  df$distance[df$distance > 1600] <- 1600
  df$Region.Label = rep(1,dim(df)[1])
  df$Sample.Label = rep(1,dim(df)[1])
  try({
    model <- ddf(method = "io", dsmodel =~cds(key ="hr"),
                 mrmodel =~glm(link = "logit", formula = ~distance),
                 data = df, meta.data = list(width = 1600), control = list(refit = T, nrefit = 5, debug = T))
    ests <- dht(model, region.table = data.frame(Region.Label = 1, Area = 600000*2000*2),
                sample.table = data.frame(Region.Label = 1, Sample.Label = 1,Effort = 600000))
    
    N <- ests$individuals$N$Estimate
    lci <- ests$individuals$N$lcl
    uci <- ests$individuals$N$ucl
    return(c(N, lci, uci))
  })
}

fit.mrds(df = simul, mismatch = TRUE)