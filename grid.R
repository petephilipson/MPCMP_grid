mulong <- seq(0.01, 10, 0.01) # No real solution for mu >= 10 
nulong <- seq(0.01, 10, 0.01)
grid_1K <- matrix(0, length(mulong), length(nulong))
unique <- 0:10 
logfac <- lfactorial(unique)
library(parallel)
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores) 
grid_1K <- foreach (i = 1:(length(mulong) - 1), .combine = rbind) %dopar% {
  temp <- rep(0, length(nulong))
  for (j in 1:length(nulong))  {
    roots <- polyroot(((0:10) - mulong[i])/(factorial(0:10)^nulong[j]))
    lam <- Re(roots[Re(roots) >= 0 & zapsmall(Im(roots), 2) ==  0])
    if (length(lam) == 0){lam <- grid_1K[i - 1, j]}
    temp[j] <- lam
  }
  temp
}