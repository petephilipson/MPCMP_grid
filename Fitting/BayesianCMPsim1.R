# Code to fit Bayesian MPCMP model for simulation 1
library(mpcmp)
library(mvtnorm)

source("sim_studies.R")

#### Simulation loop across three methods ####
run_sims <- function(Nsims = 2, n = 150){
  q <- c(0.005, 0.995, 0.025, 0.975, 0.05, 0.95)
  fit1res <- matrix(NA, nrow = Nsims, ncol = 13)
  fit2res <- matrix(NA, nrow = Nsims, ncol = 13)
  fit3res <- matrix(NA, nrow = Nsims, ncol = 13)
  fit1cred <- matrix(NA, nrow = Nsims, ncol = 11*6)
  fit2cred <- matrix(NA, nrow = Nsims, ncol = 11*6)
  fit3cred <- matrix(NA, nrow = Nsims, ncol = 11*6)
  for (i in 1:Nsims){
    print(i)
    simdat <- sim1(n = n)$data
    fit1 <- fit_huang_sim1(simdat, niters = 1000)
    fit2 <- fit_grid_sim1(simdat, niters = 5000)
    #fit3 <- fit_grid_2(simdat)
    fit1res[i,] <- c(i, colMeans(fit1$paras), fit1$cpu_time)
    fit1cred[i, ] <- as.vector(t(summary(as.mcmc(fit1$paras), quantiles = q)[[2]]))
    fit2res[i,] <- c(i, colMeans(fit2$paras), fit2$cpu_time)
    fit2cred[i, ] <- as.vector(t(summary(as.mcmc(fit2$paras), quantiles = q)[[2]]))
    #fit3res[i,] <- c(i, colMeans(fit3$paras), fit3$cpu_time)
    #fit3cred[i, ] <- as.vector(t(summary(as.mcmc(fit3$paras), quantiles = q)[[2]]))
  }
  list(fit1res = fit1res, fit2res = fit2res, fit3res = fit3res,
       fit1cred = fit1cred, fit2cred = fit2cred, fit3cred = fit3cred)
}

#### Code to fit Huang approach ####
fit_huang_sim1 <- function(simdat, niters = 1000){
fit.original <- glm.cmp(y ~ . - 1, simdat)

# Extract design matrix and response
x.matrix <- fit.original$x
y <- simdat$y

# Fitted betas
Beta.hat <- fit.original$coefficients

# Fitted dispersion
nu.hat <- fit.original$nu

# Covariance matrix for proposal distribution
Sigma_B <- 0.6*fit.original$variance_beta

# Prior mean and variance for Beta
p <- dim(x.matrix)[2]
mu0 <- rep(0,p)
Sigma0 <- diag(p)*10^5

# Initial values for Beta
Beta0 <- Beta.hat

# Initial value for nu
nu0 <- nu.hat

means0 <- exp(x.matrix%*%Beta0)
loglike_0 <- sum(dcomp(y, means0, nu0, log.p = T, summax = 11)) #+ dmvnorm(Beta0, mu0, Sigma0, log=T)
# remove comment above for normal prior, otherwise flat prior

# K is number of MCMC samples
K <- niters
Betas <- matrix(, nrow = K, ncol = p)
nus <- matrix(,nrow=K, ncol = 1)
colnames(Betas) <-  paste("x", 1:10, sep = "")
colnames(nus) <- "dispersion"
iter_times <- rep(0, K)
ptm <- proc.time()
for(k in 1:K){
  start <- proc.time()[3]
  if(k%%100==0) {
    cat(paste0("iteration: ", k, "\n"))
  }
  # Beta update
  Beta1 <- rmvnorm(1, Beta0, Sigma_B)
  means1 <- exp(x.matrix%*%t(Beta1))
  loglike_1 <- sum(dcomp(y, means1, nu0, log.p = T, summax = 11)) #+ dmvnorm(Beta1, mu0, Sigma0, log=T)
  # remove comment above for normal prior, otherwise flat prior
  loglikediff <- loglike_1 - loglike_0
  logprob <- log(runif(1))
  if(logprob < loglikediff){
    Beta0 <- Beta1
    loglike_0 <- loglike_1
    means0 <- means1
  }
  Betas[k,] <- Beta0
  
  # Dispersion update
  nu1 <- rexp(1, 1/nu0)
  loglike_1 <- sum(dcomp(y, means0, nu1, log.p = T))
  nu.ratio <- nu1/nu0
  loglikediff <- loglike_1 - loglike_0 + log(nu.ratio) + nu.ratio - 1/nu.ratio
  logprob <- log(runif(1))
  if(logprob < loglikediff){
    nu0 <- nu1
    loglike_0 <- loglike_1
  }
  nus[k,] <- nu0
  iter_times[k] <- proc.time()[3] - start
}
paras <- cbind(Betas, nus, iter_times)
cpu_time <- proc.time() - ptm
list(cpu_time = cpu_time[3], paras = paras)
}

#### Fit grid for lambda only ####
fit_grid_sim1 <- function(simdat, niters = 1000){
  ## Fit full model
  fit.original <- glm.cmp(y ~ . - 1, simdat)

  # Extract design matrix and response
  x.matrix <- fit.original$x
  y <- simdat$y
  
  # Fitted betas
  Beta.hat <- fit.original$coefficients
  
  # Fitted dispersion
  nu.hat <- fit.original$nu
  
  # Covariance matrix for proposal distribution
  Sigma_B <- 0.7*fit.original$variance_beta
  
  # Prior mean and variance for Beta
  p <- dim(x.matrix)[2]
  mu0 <- rep(0,p)
  Sigma0 <- diag(p)*10^5
  
  # Initial value for Beta
  Beta0 <- Beta.hat
  mu0 <- exp(x.matrix %*% Beta0)
  if (max(mu0) >= 10){mu0 <- pmin(9.9, mu0)}
  if (min(mu0) < 0.001){mu0 <- pmax(0.001, mu0)}
  # Initial value for nu
  nu0 <- nu.hat
  
  # K is number of MCMC samples
  K <- niters
  Betas <- matrix(, nrow = K, ncol = p)
  nus <- matrix(,nrow=K, ncol = 1)
  colnames(Betas) <-  paste("x", 1:10, sep = "")
  colnames(nus) <- "dispersion"
  iter_times <- rep(0, K)
  ptm <- proc.time()
  
  # Some additional preparation
  summ <- summarydata(y)
  lfac <- summ$lfac
  lfac_y <- summ$lfac_y
  N <- summ$N
  unique <- summ$unique
  sum_y <- summ$sum_y
  
  for(k in 1:K){
    start <- proc.time()[3]
    if(k%%1000==0) {
    cat(paste0("iteration: ", k, "\n"))
    }
    # Beta update
    Beta1 <- rmvnorm(1, Beta0, Sigma_B)
    mu1 <- exp(x.matrix%*%t(Beta1))
    if (max(mu1) >= 10){mu1 <- pmin(9.9, mu1)}
    if (min(mu1) < 0.001){mu1 <- pmax(0.001, mu1)}
    #Grid part
    log_lambda <- log(lam_grid_poly_10K[mu0*1000, nu0*1000]) 
    log_lambda1 <- log(lam_grid_poly_10K[mu1*1000, nu0*1000])
    # Calculate ln(G) directly (depends on lambda)
    if(length(log_lambda) != N){
      print(mu0)
      print(nu0)
      print(log_lambda)
    }
    lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(rep(nu0, N), lfac)))) 
    lG1 <- log(rowSums(exp(tcrossprod(log_lambda1, unique) - tcrossprod(rep(nu0, N), lfac)))) 
    denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      Beta0 <- Beta1
      mu0 <- mu1
    }
    Betas[k,] = Beta0
    
    # Dispersion update
    nu1 <- rexp(1, 1/nu0)
    if (nu1 < 0.001){nu1 <- 0.001}
    if (nu1 > 9.9){nu1 <- 9.9}
    log_lambda <- log(lam_grid_poly_10K[mu0*1000, nu0*1000]) 
    log_lambda1 <- log(lam_grid_poly_10K[mu0*1000, nu1*1000])
    lG <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(rep(nu0, N), lfac))))
    lG1 <- log(rowSums(exp(tcrossprod(log_lambda1, unique) - tcrossprod(rep(nu1, N), lfac))))
    denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG) + (nu0 - nu1)*lfac_y)
    nu.ratio <- nu1/nu0
    denrat <- denrat + log(nu.ratio) + nu.ratio - 1/nu.ratio
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      nu0 <- nu1
    }
    nus[k,] <- nu0
    iter_times[k] <- proc.time()[3] - start
  }
  
  paras <- cbind(Betas, nus)
  cpu_time <- proc.time() - ptm
  list(cpu_time = cpu_time[3], paras = paras, iter_times = iter_times)
}

#### Function to extract some data summaries needed in updates
summarydata <- function(y)
{
  N <- length(y)
  unique <- 0:max(y) #sort(unique(data[, col]))
  lfac <- lfactorial(unique)
  lfac_y <- lfactorial(y)
  sum_lfac_y <- sum(lfactorial(y))
  sum_y <- sum(y)
  list(N = N, unique = unique, lfac = lfac, lfac_y = lfac_y, sum_lfac_y = sum_lfac_y, 
       sum_y = sum_y)
}


#### Fit grid for lambda and normalising constant ####
fit_grid_2 <- function(simdat, niters = 1000){
  # Grid with lambda and lnG
  # Fit full model
  fit.original <- glm.cmp(y ~ . - 1, simdat)
  
  # Extract design matrix and response
  x.matrix = fit.original$x
  y = simdat$y
  
  # Fitted betas
  Beta.hat = fit.original$coefficients
  
  # Fitted dispersion
  nu.hat = fit.original$nu
  
  # Covariance matrix for proposal distribution
  Sigma_B = fit.original$variance_beta
  
  # Prior mean and variance for Beta
  p = dim(x.matrix)[2]
  mu0 = rep(0,p)
  Sigma0 = diag(p)*10^5
  
  # Initial value for Beta
  Beta0 = Beta.hat
  mu0 <- exp(x.matrix %*% Beta0)
  if (max(mu0) >= 10){mu0 <- pmin(9.9, mu0)}
  if (min(mu0) < 0.001){mu0 <- pmax(0.001, mu0)}
  # Initial value for nu
  nu0 = nu.hat
  
  K = niters
  Betas = matrix(, nrow = K, ncol = p)
  nus = matrix(,nrow=K, ncol = 1)
  colnames(Betas) =  paste("x", 1:10, sep = "")
  colnames(nus) = "dispersion"
  ptm <- proc.time()
  
  # Some additional preparation
  summ <- summarydata(simdat, col = 1)
  lfac <- summ$lfac
  lfac_y <- summ$lfac_y
  N <- summ$N
  unique <- summ$unique
  sum_y <- summ$sum_y
  
  for(k in 1:K){
    
    if(k%%1000==0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    # Beta update
    Beta1 = rmvnorm(1, Beta0, Sigma_B)
    mu1 = exp(x.matrix%*%t(Beta1))
    if (max(mu1) >= 10){mu1 <- pmin(9.9, mu1)}
    if (min(mu1) < 0.001){mu1 <- pmax(0.001, mu1)}
    
    # Grid part
    log_lambda <- log(lam_grid_poly_10K[mu0*1000, nu0*1000]) # 0.001 is stepsize (1000 = 1/0.001)
    log_lambda1 <- log(lam_grid_poly_10K[mu1*1000, nu0*1000])
    # Find ln(G) on grid - depends on lambda, nu
    lG <- lnG_grid_poly[mu0*1000, nu0*1000]
    lG1 <- lnG_grid_poly[mu1*1000, nu0*1000]
    # End of grid part
    denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      Beta0 <- Beta1
      mu0 <- mu1
    }
    Betas[k,] = Beta0
    
    # Dispersion update
    nu1 = rexp(1, 1/nu0)
    if (nu1 < 0.001){nu1 <- 0.001}
    if (nu1 > 9.9){nu1 <- 9.9}
    log_lambda <- log(lam_grid_poly_10K[mu0*1000, nu0*1000]) 
    log_lambda1 <- log(lam_grid_poly_10K[mu0*1000, nu1*1000])
    # Find ln(G) on grid - depends on lambda, nu
    lG <- lnG_grid_poly[mu0*1000, nu0*1000]
    lG1 <- lnG_grid_poly[mu0*1000, nu1*1000]
    denrat <- sum(y*(log_lambda1 - log_lambda) - (lG1 - lG) + (nu0 - nu1)*lfac_y)
    nu.ratio = nu1/nu0
    denrat <- denrat + log(nu.ratio) + nu.ratio - 1/nu.ratio
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      nu0 <- nu1
    }
    nus[k,] = nu0
      }
  paras <- cbind(Betas, nus)
  cpu_time <- proc.time() - ptm
  list(cpu_time = cpu_time[3], paras = paras)
}
  
#### Debugging computation of lambda in Huang approach (via mpcmp package) ####
comp_lambdas <- function(x, mu, nu, lambdalb = 1e-10, lambdaub = 1000,
                         maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1, summax = 10) {
  df <- CBIND(mu = mu, nu = nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  mu <- df[, 1]
  nu <- df[, 2]
  lambda <- df[, 3]
  lambdalb <- df[, 4]
  lambdaub <- df[, 5]
  lambda.ok <-
    comp_lambdas_fixed_ub(x, mu, nu,
                          lambdalb = lambdalb, lambdaub = lambdaub,
                          maxlambdaiter = maxlambdaiter, tol = tol,
                          lambdaint = lambda, summax = summax
    )
  lambda <- lambda.ok$lambda
  lambdaub <- lambda.ok$lambdaub
  lambdaub.err.ind <- which(is.nan(lambda))
  sub_iter1 <- 1
  while (length(lambdaub.err.ind) > 0 && sub_iter1 <= 100) {
    lambdaub[lambdaub.err.ind] <- 0.5 * lambdaub[lambdaub.err.ind]
    lambda.ok <-
      comp_lambdas_fixed_ub(x, mu[lambdaub.err.ind], nu[lambdaub.err.ind],
                            lambdalb = lambdalb[lambdaub.err.ind],
                            lambdaub = lambdaub[lambdaub.err.ind],
                            maxlambdaiter = maxlambdaiter, tol = tol,
                            lambdaint = lambda[lambdaub.err.ind], summax = summax
      )
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1 + 1
    lambdaub.err.ind <- which(is.nan(lambda))
  }
  lambdaub.err.ind <- which(lambda / lambdaub >= 1 - tol)
  sub_iter1 <- 1
  while (length(lambdaub.err.ind) > 0 && sub_iter1 <= 100) {
    #print(lambdaub.err.ind)
    #print(c(length(lambdaub.err.ind), sub_iter1))
    lambdaub[lambdaub.err.ind] <- 2^(sub_iter1) * lambdaub[lambdaub.err.ind]
    lambda.ok <-
      comp_lambdas_fixed_ub(x, mu[lambdaub.err.ind], nu[lambdaub.err.ind],
                            lambdalb = lambdalb[lambdaub.err.ind],
                            lambdaub = lambdaub[lambdaub.err.ind],
                            maxlambdaiter = maxlambdaiter, tol = tol,
                            lambdaint = lambda[lambdaub.err.ind],
                            summax = summax
      )
    lambda[lambdaub.err.ind] <- lambda.ok$lambda
    lambdaub[lambdaub.err.ind] <- lambda.ok$lambdaub
    sub_iter1 <- sub_iter1 + 1
    lambdaub.err.ind <- which(lambda / lambdaub >= 1 - tol)
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- max(lambdaub)
  return(out)
}

comp_lambdas_fixed_ub <- function(x, mu, nu, lambdalb = 1e-10, lambdaub = 1000,
                                  maxlambdaiter = 1e3, tol = 1e-6, lambdaint = 1,
                                  summax = 10) {
  # df <- CBIND(mu=mu, nu=nu, lambda = lambdaint, lb = lambdalb, ub = lambdaub)
  # mu <- df[,1]
  # nu <- df[,2]
  # lambda <- df[,3]
  # lb <- df[,4]
  # ub <- df[,5]
  lambda <- lambdaint
  lb <- lambdalb
  ub <- lambdaub
  iter <- 1
  log.Z <- logZ(x, log(lambda), nu, summax = summax)
  mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
  not.converge.ind <- which(abs(mean1 - mu) > tol)
  while (length(not.converge.ind) > 0 && iter < 200) {
    still.above.target.ind <- which((mean1[not.converge.ind]
                                     > mu[not.converge.ind]))
    still.below.target.ind <- which(mean1[not.converge.ind] < mu[not.converge.ind])
    lb[not.converge.ind[still.below.target.ind]] <-
      lambda[not.converge.ind[still.below.target.ind]]
    ub[not.converge.ind[still.above.target.ind]] <-
      lambda[not.converge.ind[still.above.target.ind]]
    lambda <- (lb + ub) / 2
    log.Z <- logZ(x, log(lambda), nu, summax = summax)
    mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    while (sum(mean1 == 0) > 0) {
      ub[not.converge.ind[mean1 == 0]] <- ub[not.converge.ind[mean1 == 0]] / 2
      lambdaub <- lambdaub / 2
      lambda <- (lb + ub) / 2
      log.Z <- logZ(x, log(lambda), nu, summax = summax)
      mean1 <- comp_means(lambda, nu, log.Z = log.Z, summax = summax)
    }
    not.converge.ind <- which((1 - (((abs(mean1 - mu) <= tol) + (lambda == lb) + (lambda == ub)
                                     + (ub == lb)) >= 1)) == 1)
    iter <- iter + 1
  }
  while (length(not.converge.ind) > 0 && iter < maxlambdaiter) {
    # basically the content of comp_variances without recalculating Z and mean1
    term <- matrix(0, nrow = length(mu), ncol = summax)
    for (y in 1:summax) {
      term[, y] <- exp(2 * log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
    }
    var1 <- apply(term, 1, sum) - mean1^2
    ## newton raphson update
    newtonsteps <- -lambda[not.converge.ind] * mean1[not.converge.ind] /
      (var1[not.converge.ind])^2 * (log(mean1[not.converge.ind]) - log(mu[not.converge.ind]))
    
    
    lambda.new <- lambda[not.converge.ind] + newtonsteps
    ## if newton raphson steps out of bound, use bisection method
    out.of.bound.ind <- which((lambda.new < lb[not.converge.ind])
                              + (lambda.new > ub[not.converge.ind]) == 1)
    if (length(out.of.bound.ind > 0)) {
      lambda.new[out.of.bound.ind] <-
        (lb[not.converge.ind[out.of.bound.ind]] + ub[not.converge.ind[out.of.bound.ind]]) / 2
      # any out of bound updates are replaced with mid-point of ub and lb
    }
    lambda[not.converge.ind] <- lambda.new
    log.Z <- logZ(x, log(lambda), nu, summax)
    term <- matrix(0, nrow = length(mu), ncol = summax)
    for (y in 1:summax) {
      term[, y] <- exp(log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
    }
    mean1 <- apply(term, 1, sum)
    if (length(out.of.bound.ind) > 0) {
      still.above.target.ind <- which(mean1[not.converge.ind[out.of.bound.ind]]
                                      > mu[not.converge.ind[out.of.bound.ind]])
      still.below.target.ind <- which(mean1[out.of.bound.ind] < mu[out.of.bound.ind])
      if (length(still.below.target.ind) > 0) {
        lb[not.converge.ind[out.of.bound.ind[still.below.target.ind]]] <-
          lambda[not.converge.ind[out.of.bound.ind[still.below.target.ind]]]
      }
      if (length(still.above.target.ind) > 0) {
        ub[not.converge.ind[out.of.bound.ind[still.above.target.ind]]] <-
          lambda[not.converge.ind[out.of.bound.ind[still.above.target.ind]]]
        # any out of bound updates are replaced with mid-point of ub and lb
      }
    }
    not.converge.ind <- which((1 - (((abs(mean1 - mu) <= tol) + (lambda == lb) + (lambda == ub)
                                     + (ub == lb)) >= 1)) == 1)
    iter <- iter + 1
  }
  out <- list()
  out$lambda <- lambda
  out$lambdaub <- lambdaub
  return(out)
}

logZ <- function(y, log_lambda, nu, summax = 10) {
  # approximates normalizing constant for COMP distributions
  # lambda, nu are recycled to match the length of each other.
  df <- CBIND(log_lambda = log_lambda, nu = nu)
  log_lambda <- df[, 1]
  nu <- df[, 2]
  summ <- summarydata(y)
  lfac <- summ$lfac
  unique <- summ$unique
  logZ <- log(rowSums(exp(tcrossprod(log_lambda, unique) - tcrossprod(nu, lfac))))
  return(logZ)
  #return(logZ_c(log_lambda, nu, summax = summax))
}

dcomp <- function(x, mu, nu = 1, lambda, log.p = FALSE, lambdalb = 1e-10,
                  lambdaub = 1000, maxlambdaiter = 1e3, tol = 1e-6, summax = 10, dat) {
  # compute the pmf/density for COMP distirbution with mean mu and dispersion nu
  # x, mu, nu are recycled to match the length of each other.
  # lambdaub will be scaled down/up  if there is
  # over-/under-dispersion so that the correct lambda can be found
  if (missing(mu) && missing(lambda)) {
    stop('argument "mu" is missing, with no default')
  }
  if (!missing(mu) && !missing(lambda)) {
    stop("specify 'mu' or 'lambda' but not both")
  }
  not.miss.mu <- !missing(mu)
  if (missing(mu)) {
    mu <- Inf
  }
  if (missing(lambda)) {
    lambda <- Inf
  }
  df <- CBIND(x = x, mu = mu, nu = nu, lambda = lambda)
  x <- df[, 1]
  mu <- df[, 2]
  nu <- df[, 3]
  lambda <- df[, 4]
  warn <- FALSE
  if (not.miss.mu) {
    if (missing(summax)) {
      summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
    }
    mu.ok.ind <- which(mu > 0)
    mu.err.ind <- which(mu <= 0)
    if (length(mu.err.ind) > 0) {
      lambda[mu.err.ind] <- mu[mu.err.ind]
    }
    if (length(mu.ok.ind) > 0) {
      lambda.ok <- comp_lambdas(x, mu[mu.ok.ind], nu[mu.ok.ind],
                                lambdalb = lambdalb,
                                lambdaub = lambdaub,
                                # lambdaub = min(lambdaub,2*max(lambdaold))
                                maxlambdaiter = maxlambdaiter, tol = tol,
                                summax = summax
      )
      lambda[mu.ok.ind] <- lambda.ok$lambda
      lambdaub <- lambda.ok$lambdaub
    }
  } else {
    # A <- (8*nu^2+12*nu+3)/(96*nu^2*lambda^(1/nu))
    # B <- (1+6*nu)/(144*nu^3*lambda^(2/nu))
    # D <- 1+(nu-1)*(A+B)
    # mu <- rep(max(lambda^(1/nu)-(nu-1)/(2*nu)+1/D*((nu-1)*(-A/nu+2*B/nu))),
    # length(lambda))
    # mu_error <- which(is.nan(mu)>0 | mu< 0)
    mu <- comp_means(lambda, nu, summax = 10)
    if (missing(summax)) {
      #summax <- ceiling(max(c(mu + 20 * sqrt(mu / nu), 100)))
      cat("As you do not specify mu nor summax, summax will be calculated based on\n")
      cat("mu which is calcualted by truncated sum at 500.\n")
      cat("If you believe the mean of the distribution is somewhat close to 500,\n")
      cat("you may want to do some experiment with the comp_means() and\n")
      cat("specify summax instead to improve the accuracy.\n")
    }
  }
  # at a vector of yvalues
  pmf <- rep(0, length(x))
  for (i in 1:length(x)) {
    if ((mu[i] == 0 || lambda[i] == 0) && x[i] == 0) {
      pmf[i] <- 0 # log(1), 1 as the distribution is degenerated at 0
    } else if (mu[i] < 0 | lambda[i] < 0 | nu[i] <= 0) {
      pmf[i] <- NaN
      warn <- TRUE
    } else {
      if (!is.wholenumber(x[i])) {
        warning(paste("non-integer x =", x[i]))
        pmf[i] <- -Inf # log(0)
      } else {
        if (x[i] < 0) {
          pmf[i] <- -Inf
        } else { # log(0)
          # pmf <- log(density)
          pmf[i] <- x[i] * log(lambda[i]) - (nu[i] * lfactorial(x[i])) -
            logZ(x, log(lambda[i]), nu[i], summax)
        }
      }
    }
  }
  if (!log.p) {
    pmf <- exp(pmf)
  }
  if (warn) {
    warning("NaN(s) produced")
  }
  return(pmf)
}

CBIND <- function(..., deparse.level = 1) {
  dots <- list(...)
  len <- sapply(dots, length)
  dots <- lapply(seq_along(dots),
                 function(i, x, len) rep(x[[i]], length.out = len),
                 x = dots, len = max(len)
  )
  do.call(cbind, c(dots, deparse.level = deparse.level))
}

comp_means <- function(lambda, nu, log.Z, summax = 10) {
  # approximates mean by truncation of COMP distributions
  # lambda, nu, mu.bd are recycled to match the length of each other.
  if (missing(log.Z)) {
    df <- CBIND(lambda = lambda, nu = nu)
    log.Z <- logZ(log(df[, 1]), df[, 2], summax)
  }
  df <- CBIND(lambda = lambda, nu = nu, log.Z = log.Z)
  lambda <- df[, 1]
  nu <- df[, 2]
  log.Z <- df[, 3]
  term <- matrix(0, nrow = length(lambda), ncol = summax)
  for (y in 1:summax) {
    term[, y] <- exp(log(y - 1) + (y - 1) * log(lambda) - nu * lgamma(y) - log.Z)
  }
  mean1 <- apply(term, 1, sum)
  return(mean1)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
  

}

#### Poisson fit ####
fit_pois_sim1 <- function(simdat, niters = 1000){
  ## Fit full Poisson model
  fit.original <- glm(y ~ . - 1, family = "poisson", x = TRUE, simdat)
  
  # Extract design matrix and response
  x.matrix = fit.original$x
  y = simdat$y
  
  # Fitted betas
  Beta.hat = fit.original$coefficients
  
  # Covariance matrix for proposal distribution
  Sigma_B = 0.7*vcov(fit.original)
  
  # Prior mean and variance for Beta
  p = dim(x.matrix)[2]
  mu0 = rep(0,p)
  Sigma0 = diag(p)*10^5
  
  # Initial value for Beta
  Beta0 = Beta.hat
  mu0 <- exp(x.matrix %*% Beta0)
  if (max(mu0) >= 10){mu0 <- pmin(9.9, mu0)}
  if (min(mu0) < 0.001){mu0 <- pmax(0.001, mu0)}
  
  # K is number of MCMC samples
  K = niters
  Betas = matrix(, nrow = K, ncol = p)
  colnames(Betas) =  paste("x", 1:10, sep = "")
  iter_times <- rep(0, K)
  ptm <- proc.time()
  
  # Some additional preparation
  summ <- summarydata(simdat, col = 1)
  lfac <- summ$lfac
  lfac_y <- summ$lfac_y
  N <- summ$N
  unique <- summ$unique
  sum_y <- summ$sum_y
  
  for(k in 1:K){
    start <- proc.time()[3]
    if(k%%1000==0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    # Beta update
    Beta1 = rmvnorm(1, Beta0, Sigma_B)
    mu1 = exp(x.matrix%*%t(Beta1))
    denrat <- sum(y*(log(mu1) - log(mu0)) - (mu1 - mu0))
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      Beta0 <- Beta1
      mu0 <- mu1
    }
    Betas[k,] = Beta0
    
    iter_times[k] <- proc.time()[3] - start
  }
  
  paras <- cbind(Betas)
  cpu_time <- proc.time() - ptm
  list(cpu_time = cpu_time[3], paras = paras, iter_times = iter_times)
}

#### Negative binomial fit ####
fit_nb_sim1 <- function(simdat, niters = 1000){
  ## Fit full NB model
  library(MASS)
  fit.original <- glm.nb(y ~ . - 1, simdat, x = TRUE)
  
  # Extract design matrix and response
  x.matrix = fit.original$x
  y = simdat$y
  
  # Fitted betas
  Beta.hat = fit.original$coefficients
  # Fitted dispersion
  eta.hat = fit.original$theta
  
  # Covariance matrix for proposal distribution
  Sigma_B = 0.1*vcov(fit.original)
  
  # Prior mean and variance for Beta
  p = dim(x.matrix)[2]
  mu0 = rep(0,p)
  Sigma0 = diag(p)*10^5
  
  # Initial values for Beta and lambda
  Beta0 <- Beta.hat
  lam0 <- exp(x.matrix %*% Beta0)
  # Initial value for eta
  eta0 <- eta.hat
  
  # K is number of MCMC samples
  K = niters
  Betas = matrix(, nrow = K, ncol = p)
  etas = matrix(, nrow = K, ncol = 1)
  mus <- matrix(, nrow = K, ncol = nrow(x.matrix))
  colnames(Betas) =  paste("x", 1:10, sep = "")
  colnames(etas) = "dispersion"
  iter_times <- rep(0, K)
  ptm <- proc.time()
  
  # Some additional preparation
  summ <- summarydata(simdat, col = 1)
  lfac <- summ$lfac
  lfac_y <- summ$lfac_y
  N <- summ$N
  unique <- summ$unique
  sum_y <- summ$sum_y
  
  for(k in 1:K){
    start <- proc.time()[3]
    if(k%%1000==0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    # mu update (closed form)
    mu <- rgamma(N, eta0 + y, 1 + eta0/lam0)
    mus[k, ] <- mu
    
    # Beta update
    Beta1 <- rmvnorm(1, Beta0, Sigma_B)
    lam1 <- exp(x.matrix%*%t(Beta1))
    denrat <- sum(eta0*mu*(1/lam0 - 1/lam1) - (lam1 - lam0)*eta0)
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      Beta0 <- Beta1
      lam0 <- lam1
    }
    Betas[k,] = Beta0
      
    # Dispersion update
    eta1 = rexp(1, 1/eta0)
    denrat <- (eta1 - eta0)*(log(mu/lam0) - mu/lam0) + eta1*log(eta1) - eta0*log(eta0) -
      lgamma(eta1) + lgamma(eta0)
    eta.ratio = eta1/eta0
    denrat <- denrat + log(eta.ratio) + eta.ratio - 1/eta.ratio
    laccept <- min(0, denrat)
    accept <- (log(runif(1)) < laccept)
    if (accept){
      eta0 <- eta1
    }
    etas[k,] = eta0
    iter_times[k] <- proc.time()[3] - start
  }
  
  paras <- cbind(Betas, etas, mus)
  cpu_time <- proc.time() - ptm
  list(cpu_time = cpu_time[3], paras = paras, iter_times = iter_times)
}
