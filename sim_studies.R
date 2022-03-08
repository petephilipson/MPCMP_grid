#### Simulation studies for MPCMP grid paper ####
library(mpcmp)
#### Simulation 1: n = 150 with 9 covariates, scalar dispersion - similar to bids data
# 5 binary covariates (x1 - x4, x9), 4 continuous (x5 - x8)
# Count response in (0-10), mode, median of 1 mean of 1.74, variance of 2.05
sim1 <- function(n){
# Generate all binary covariates with 50/50 chance of 0/1
xbin <- rbinom(n*5, 1, 0.5)
x1 <- xbin[1:n]
x2 <- xbin[(n+1):(2*n)]
x3 <- xbin[(2*n+1):(3*n)]
x4 <- xbin[(3*n+1):(4*n)]
x9 <- xbin[(4*n+1):(5*n)]
x5 <- rnorm(n, 1.35, 0.05)
x6 <- rlnorm(n, -1.7, 1)
x7 <- rlnorm(n, -1.1, 1.5)
x8 <- x7^2
X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9)

beta <- c(1, 0.25, -0.20, 0.05, 0.50, -0.70, -0.40, 0.20, -0.01, -0.05)
mu <- exp(X%*%beta)
y <- rcomp(n, mu, nu = 1.5)
data <- data.frame(y = y, X = X)
colnames(data)[-1] =  paste("x", 1:10, sep = "")
list(data = data)
}

#### Simulation 2: n = 2000
sim2 <- function(n = 2000, nref = 20){
  # Generate two binary covariates (home and no fans) with 50/50 chance of 0/1
  # Referee-specific dispersions
  xbin <- rbinom(n*2, 1, 0.5)
  # Partition into two variables and create their interaction
  x1 <- xbin[1:n]
  x2 <- xbin[(n+1):(2*n)]
  x3 <- x1 * x2
  X <- cbind(1, x1, x2, x3)
  # Simulate referees
  ref <- rev(seq(0, 0.95, length.out = nref))
  # Corner constraint by setting first ref as zero but fit as an STZ constraint
  ref_long <- rep(ref, n/nref)
  beta <- c(0.525, -0.10, -0.20, 0.10)
  mu <- exp(X%*%beta + ref_long)
  # Dispersion by referee: half under- and half overdispersed
  nu_ref <- c(rep(1.25, nref/2), rep(0.8, nref/2))
  nu_ref_long <- rep(nu_ref, n/nref)
  y <- rcomp(n, mu, nu = nu_ref_long)
  ref_id <- as.factor(rep(1:nref, n/nref))
  data <- data.frame(y = y, X = X, ref_id = ref_id, ref_abil = ref_long, ref_disp = nu_ref_long)
  colnames(data)[2:5] =  paste("x", 1:4, sep = "")
  list(data = data)
}

#### Simulation 3: n = 15000 - similar to cricket data
# Taking players with >= 150 Test wickets gives 111 players with average of 125 observations
# Include home, innings effects
# Player-specific dispersions
sim3 <- function(nplayers = 100, ninns = 150){
  n <- nplayers * ninns
  # Generate binary covariate (home) with 50/50 chance of 0/1
  x1 <- rbinom(n, 1, 0.5)
  # Generate factor for innings
  inns_lab <- sample(1:4, n, replace = TRUE)
  # Collect in design matrix
  X <- model.matrix(~ x1 + as.factor(inns_lab))
  # Simulate players
  player <- seq(0.95, 0, length.out = nplayers) 
  # Set as an STZ constraint
  player <- player - mean(player)
  player_long <- rep(player, ninns)
  beta <- c(0.8, -0.10, 0, 0.05, 0.10) #c(log(2.5), -0.10, 0.05, 0.10, 0.15)
  mu <- exp(X%*%beta + player_long)
  # Dispersion by player  - now cross the dispersion line in a sequence
  nu_player <- seq(0.8, 1.25, length.out = nplayers)
  nu_player_long <- rep(nu_player, ninns)
  y <- rcomp(n, mu, nu = nu_player_long)
  player_id <- as.factor(rep(1:nplayers, ninns))
  data <- data.frame(y = y, X = X, player_id = player_id, player_abil = player_long, 
                     player_disp = nu_player_long)
  colnames(data)[2:6] =  paste("x", 0:4, sep = "")
  list(data = data)
}


