
n = 1000
set.seed(23)
moran_statistic <- 
  mclapply(1:100, mc.cores = 20, function(i){
    x <- rnorm(n, sd = 2)
    y <- rnorm(n, sd = 2)
    v <- rnorm(n, sd = 2)
    dat <- tibble(x, y, v)
    d <- as.matrix(dist(dat[, c('x', 'y')]))
    cov_mat <- exp(-d^2 / sigma)
    
    # # Check for negative eigenvalues
    # eigenvalues <- eigen(cov_mat)$values
    # if(any(eigenvalues < 0)) message("Negative eigenvalue detected in covariance matrix")
    
    # Spatial confounder s, sampled from multivariate normal distribution
    # for now multiply by 10 
    s <- mvtnorm::rmvnorm(1, sigma = 2 * cov_mat)
    
    # Calculate Moran's I statistic to see how spatially correlated it is
    location_dists_inv <- 1/d
    diag(location_dists_inv) <- 0
    moran_test <- Moran.I(as.vector(s), location_dists_inv)
    return(moran_test$observed)
  })

moran_dist <- unlist(moran_statistic)
hist(moran_dist)
