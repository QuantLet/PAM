# -------------------------------------------------------------------------------
# Excess bond risk premia modelling with PAM
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx",
              "lars", "scales", "tilting")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("OnestepSCADpen_withBIC_norm.r")

# Simulation setup
n.obs       = 500               # No of observations
n.par       = 10                # No  of parameters
n.sim       = 1000              # No of simulations
seed1       = 20171110          # Seed simulation X
seed2       = 20170602          # Seed simulation epsilon
r           = 0.5               # Correlation parameter of X
sd.eps      = 1                 # Standard deviation of epsilon
a           = 3.7               # Recommended value of parameter a for SCAD
n.boot      = 1000 
M           = 250
K           = n.obs/M

# Initiate cluster for parallel computing
n.cores = detectCores()           # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of scenarios in every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h.simcores = n.sim %% n.cores
if (h.simcores != 0){
  n.simcores[1:h.simcores] = n.simcores[1:h.simcores] + 1
}

# True beta coefficients (homogeneous for all t = 1, ..., 1000)
tmp1  = rep(1, 3)
tmp2  = rep(0, n.par - length(tmp1))
b     = c(tmp1, tmp2)

# Simulation of the design matrix
mu    = rep(0, n.par)
r     = 0.5
Sigma = matrix(0, nrow = n.par, ncol = n.par)

for (im in 1:n.par) {
  for (jm in 1:n.par) {
    if (im == jm){
      Sigma[im, jm] = 1
    }else {
      Sigma[im, jm] = r^abs(im - jm)
    }
  }
}

# Signal-to-noise ratio
# (t(b) %*% Sigma %*% b) / (sd.eps^2)

X = list()
for (ix in 1:n.sim){
  set.seed(seed1)
  X[[ix]] = mvrnorm(n = n.obs, mu, Sigma)
}  

# Simulation of the error term for t = 1, ..., n.obs
eps  = list()
set.seed(seed2)
for (ie in 1:n.sim){
  eps[[ie]] = rnorm(n.obs, mean = 0, sd = sd.eps)
} 

# Computation of Y for t = 1, ..., n.obs
Y    = list()
for (iy in 1:n.sim){
  Y.tmp = numeric(0)
  for (jy in 1:n.obs){
    Y.tmp = c(Y.tmp, b %*% X[[iy]][jy, ] + eps[[iy]][jy])
  }
  Y[[iy]] = Y.tmp
}

# Simulating mutlipliers from a bounded distribution
runif.mod = function(n){
  unif1    = runif(n, min = 0, max = 12)
  unif.a   = unif1[which(unif1 <= 1)]
  index.a  = which(unif1 <= 1)
  unif.b   = unif1[which((1 < unif1) & (unif1 <= 4))]
  index.b  = which((1 < unif1) & (unif1 <= 4))
  unif.c   = (unif1[which(unif1 > 4)] - 4)/8
  index.c  = which(unif1 > 4)
  unif2    = cbind(c(index.a, index.b, index.c), c(unif.a, unif.b, unif.c))
  unif.fin = unif2[order(unif2[, 1], decreasing = FALSE), 2]
  unif.fin
}

# Funtion to compute bootstrapped penalized log-likelihood function
loglik.pen.MB = function(beta.0, beta, lambda, X, Y, multipl){
  n.obs   = length(Y)
  
  X = X * sqrt(multipl)
  Y = Y * sqrt(multipl)
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm
  Y       = Y.norm
  
  loglik1 = -(- (1/(n.obs * 2 * sd.eps^2)  *  t(Y - X.cent %*% beta) %*% (Y - X.cent %*% beta))
              - (t(pen.prime(beta.0 * X.norm0, lambda, a)) %*% abs(beta * X.norm0)))
  loglik1
}

# Function to compute bootstrapped unpenalized log-likelihood function
loglik.MB = function(beta, lambda, X, Y, multipl){
  n.obs   = length(Y)
  
  X = X * sqrt(multipl)
  Y = Y * sqrt(multipl)
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm
  Y       = Y.norm
  
  loglik1 = -(- (1/(n.obs * 2 * sd.eps^2) *  t(Y - X.cent %*% beta) %*% (Y - X.cent %*% beta)))
  loglik1
}

# Function to compute penalized log-likelihood function
loglik.pen = function(beta.0, beta, lambda, X, Y){
  n.obs   = length(Y)
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm
  Y       = Y.norm
  
  loglik1 = -(- (1/(n.obs * 2 * sd.eps^2) *  t(Y - X.cent %*% beta) %*% (Y - X.cent %*% beta))
              - (t(pen.prime(beta.0 * X.norm0, lambda, a)) %*% abs(beta * X.norm0)))
  loglik1
}

# Function to compute unpenalized log-likelihood function
loglik = function(beta, lambda, X, Y){
  n.obs   = length(Y)
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm
  Y       = Y.norm
  
  loglik1 = -(- (1/(n.obs * 2 * sd.eps^2) * t(Y - X.cent %*% beta) %*% (Y - X.cent %*% beta)))
  loglik1
}

# Function returning real likelihood ratio and its bootstrapped counterpart 
# (whole distribution)
sim.LR.MBLR = function(h){
  s1 = ifelse(h == 1, 1, sum(n.simcores[1:(h - 1)],1))
  s2 = sum(n.simcores[1:h])
  real.LR         = numeric(0)
  real.LR.pen     = numeric(0)
  boot.LR.fix     = list()
  boot.LR.fix.pen = list()
  for (s in s1:s2){
    k           = 1
    k1          = 0
    k2          = 1
    beta        = list()
    lambda      = list()
    b.pts.tmp   = numeric(0)
    sim.no      = s - s1 + 1
    
    # Fit the model for the left interval I_t^(k)
    X.tmp.L      = X[[s]][((k1 * M) + 1):(k2 * M), ]
    Y.tmp.L      = Y[[s]][((k1 * M) + 1):(k2 * M)]
    n.obs.tmp.L  = length(Y.tmp.L)
    object.tmp.L = Adapt.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
    beta.tmp.L   = object.tmp.L$beta
    lambda.tmp.L = object.tmp.L$lambda
    beta0.tmp.L  = object.tmp.L$beta.0
    
    # Fit the model over the right interval I_t^(k + 1) - I_t^(k)
    X.tmp.R      = X[[s]][((k2 * M) + 1):((k2 + 1) * M), ]
    Y.tmp.R      = Y[[s]][((k2 * M) + 1):((k2 + 1) * M)]
    n.obs.tmp.R  = length(Y.tmp.R)
    object.tmp.R = Adapt.SCAD(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R)
    beta.tmp.R   = object.tmp.R$beta
    lambda.tmp.R = object.tmp.R$lambda
    beta0.tmp.R  = object.tmp.R$beta.0
    
    # Fit the model over the whole interval I_t^(k + 1)
    X.tmp.T      = rbind(X.tmp.L, X.tmp.R)
    Y.tmp.T      = c(Y.tmp.L, Y.tmp.R)
    n.obs.tmp.T  = n.obs.tmp.L + n.obs.tmp.R
    object.tmp.T = Adapt.SCAD(X.tmp.T, Y.tmp.T, a, n.obs.tmp.T)
    beta.tmp.T   = object.tmp.T$beta
    lambda.tmp.T = object.tmp.T$lambda
    beta0.tmp.T  = object.tmp.T$beta.0
    
    # Evaluate the test statistic (real penalized likelihood ratio)
    lik.L.pen = -loglik.pen(beta0.tmp.L, beta.tmp.L, lambda.tmp.L, X.tmp.L, Y.tmp.L)
    lik.R.pen = -loglik.pen(beta0.tmp.R, beta.tmp.R, lambda.tmp.R, X.tmp.R, Y.tmp.R)
    lik.T.pen = -loglik.pen(beta0.tmp.T, beta.tmp.T, lambda.tmp.T, X.tmp.T, Y.tmp.T)
    
    # Evaluate the test statistic (real unpenalized likelihood ratio)
    lik.L     = -loglik(beta.tmp.L, lambda.tmp.L, X.tmp.L, Y.tmp.L)
    lik.R     = -loglik(beta.tmp.R, lambda.tmp.R, X.tmp.R, Y.tmp.R)
    lik.T     = -loglik(beta.tmp.T, lambda.tmp.T, X.tmp.T, Y.tmp.T)
    
    real.LR.pen[sim.no] = ((n.obs.tmp.L/n.obs.tmp.T) * lik.L.pen 
                          + (n.obs.tmp.R/n.obs.tmp.T) * lik.R.pen - lik.T.pen)
    
    real.LR[sim.no] = ((n.obs.tmp.L/n.obs.tmp.T) * lik.L 
                              + (n.obs.tmp.R/n.obs.tmp.T) * lik.R - lik.T)
 
    # Find quantiles of LR with use of MB
    lik.ratio.fix     = numeric(0)
    # lik.ratio.var     = numeric(0)
    lik.ratio.fix.pen = numeric(0)
    # lik.ratio.var.pen = numeric(0)
    
    # Simulation of multipliers u_i (i = 1, ..., n.obs.tmp.T)
    set.seed      = 20170424 * s * 10
    multipliers   = matrix(runif.mod(n.boot * n.obs.tmp.T),
                           ncol = n.obs.tmp.T, nrow = n.boot)
    
    # multipliers = matrix(rexp(n.boot * n.obs.tmp.T, rate = 1), # exponential
    #                      ncol = n.obs.tmp.T, nrow = n.boot)
    
    # multipliers = matrix(rpois(n.boot * n.obs.tmp.T, lambda = 1), # exponential
    #                      ncol = n.obs.tmp.T, nrow = n.boot)
    
    multipliers.L = multipliers[, 1:n.obs.tmp.L]
    multipliers.R = multipliers[, (n.obs.tmp.L + 1):n.obs.tmp.T]

    for (l in 1:(n.boot)){
      # Multiplier bootstrap for left-hand interval I_t^(k)
      multipl.L         = multipliers.L[l, ]
      
      boot.object.fix.L = Adapt.SCAD.MB(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L, 
                                    as.numeric(lambda.tmp.L), multipl.L)
      
      boot.beta.fix.L   = boot.object.fix.L$beta
      boot.beta0.fix.L  = boot.object.fix.L$beta.0
      
      
      lik.boot.fix.L     = -loglik.MB(boot.beta.fix.L, as.numeric(lambda.tmp.L), 
                                      X.tmp.L, Y.tmp.L, multipl.L)
      lik.boot.fix.L.pen = -loglik.pen.MB(boot.beta0.fix.L, boot.beta.fix.L, as.numeric(lambda.tmp.L), 
                                       X.tmp.L, Y.tmp.L, multipl.L)
      
      # Multiplier bootstrap for right-hand interval I_t^(k + 1) - I_t^(k)
      multipl.R         = multipliers.R[l, ]
      
      boot.object.fix.R = Adapt.SCAD.MB(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R, 
                                    as.numeric(lambda.tmp.R), multipl.R)
      
      boot.beta.fix.R   = boot.object.fix.R$beta
      boot.beta0.fix.R  = boot.object.fix.R$beta.0
      
      lik.boot.fix.R     = -loglik.MB(boot.beta.fix.R, as.numeric(lambda.tmp.R), 
                                      X.tmp.R, Y.tmp.R, multipl.R)
      lik.boot.fix.R.pen = -loglik.pen.MB(boot.beta0.fix.R, boot.beta.fix.R, as.numeric(lambda.tmp.R), 
                                          X.tmp.R, Y.tmp.R, multipl.R)
      
      # Multiplier bootstrap for the whole interval I_t^(k + 1) (with shift)
      multipl            = multipliers[l, ]
      
      object.boot.fix.T  = Adapt.SCAD.MB.shift(X.tmp.L, Y.tmp.L, X.tmp.R, Y.tmp.R, a, 
                                              lambda.tmp.L, lambda.tmp.R, multipl.L, 
                                              multipl.R, beta.tmp.L, beta.tmp.R)
      boot.beta.fix.T    = object.boot.fix.T$beta
      boot.lambda.fix.T  = object.boot.fix.T$lambda
      boot.beta0.fix.T   = object.boot.fix.T$beta.0

      one.R              = rep(1, n.obs.tmp.R)
      X.tmp.R.mean       = drop(one.R %*% X.tmp.R)/n.obs.tmp.R
      X.tmp.R.cent       = scale(X.tmp.R, X.tmp.R.mean, FALSE)
      
      Y.boot.T           = c(Y.tmp.L, (Y.tmp.R - X.tmp.R.cent %*% (beta.tmp.R - beta.tmp.L)))
      
      lik.boot.fix.T     = -loglik.MB(boot.beta.fix.T, boot.lambda.fix.T, 
                                      X.tmp.T, Y.boot.T, multipl)
      lik.boot.fix.T.pen = -loglik.pen.MB(boot.beta0.fix.T, boot.beta.fix.T, boot.lambda.fix.T, 
                                          X.tmp.T, Y.boot.T, multipl)
      
      lik.ratio.fix      = c(lik.ratio.fix, ((n.obs.tmp.L/n.obs.tmp.T) * lik.boot.fix.L + 
                                           (n.obs.tmp.R/n.obs.tmp.T) *lik.boot.fix.R - 
                                           lik.boot.fix.T))
      
      lik.ratio.fix.pen  = c(lik.ratio.fix.pen, ((n.obs.tmp.L/n.obs.tmp.T) * lik.boot.fix.L.pen + 
                                            (n.obs.tmp.R/n.obs.tmp.T) *lik.boot.fix.R.pen - 
                                             lik.boot.fix.T.pen))
    }
    
    boot.LR.fix[[sim.no]]     = lik.ratio.fix
    boot.LR.fix.pen[[sim.no]] = lik.ratio.fix.pen
  }
  
  values = list(real.LR, real.LR.pen, boot.LR.fix, boot.LR.fix.pen)
  names(values) = c("real.LR", "real.LR.pen", "boot.LR.fix", "boot.LR.fix.pen")
  return(values)
}

# All scenarios paralelly
Sys.time()
real.boot.LR = foreach(isim = 1:n.cores, .packages = c('glmnet', 'MASS')) %dopar% sim.LR.MBLR(isim)
Sys.time()

# Close cluster
stopCluster(cl)

# Collect results
real.LR.final         = numeric(0)
real.LR.pen.final     = numeric(0)
boot.LR.fix.final     = list()
boot.LR.fix.pen.final = list()
boot.LR.var.final     = list()
boot.LR.var.pen.final = list()
for (icor in 1:n.cores){
  real.LR.final         = c(real.LR.final, real.boot.LR[[icor]]$real.LR)
  real.LR.pen.final     = c(real.LR.pen.final, real.boot.LR[[icor]]$real.LR.pen)
  boot.LR.fix.final     = c(boot.LR.fix.final, real.boot.LR[[icor]]$boot.LR.fix)
  boot.LR.fix.pen.final = c(boot.LR.fix.pen.final, real.boot.LR[[icor]]$boot.LR.fix.pen)
  boot.LR.var.final     = c(boot.LR.var.final, real.boot.LR[[icor]]$boot.LR.var)
  boot.LR.var.pen.final = c(boot.LR.var.pen.final, real.boot.LR[[icor]]$boot.LR.var.pen)
}

# Evaluate number of events real.LR <= alpha-quantile of boot.LR (90% quantile)
corr.90.fix     = 0
corr.90.fix.pen = 0
for (s in 1:n.sim){
  if (real.LR.final[s] <= quantile(boot.LR.fix.final[[s]], probs = 0.9)){
    corr.90.fix = corr.90.fix + 1
  }
  if (real.LR.pen.final[s] <= quantile(boot.LR.fix.pen.final[[s]], probs = 0.9)){
    corr.90.fix.pen = corr.90.fix.pen + 1
  }
}
corr.90.fix
corr.90.fix.pen

# Evaluate number of events real.LR <= alpha-quantile of boot.LR (95% quantile)
corr.95.fix     = 0
corr.95.fix.pen = 0
for (s in 1:n.sim){
  if (real.LR.final[s] <= quantile(boot.LR.fix.final[[s]], probs = 0.95)){
    corr.95.fix = corr.95.fix + 1
  }
  if (real.LR.pen.final[s] <= quantile(boot.LR.fix.pen.final[[s]], probs = 0.95)){
    corr.95.fix.pen = corr.95.fix.pen + 1
  }
}
corr.95.fix
corr.95.fix.pen

# Evaluate number of events real.LR <= alpha-quantile of boot.LR (97.5% quantile)
corr.975.fix     = 0
corr.975.fix.pen = 0
for (s in 1:n.sim){
  if (real.LR.final[s] <= quantile(boot.LR.fix.final[[s]], probs = 0.975)){
    corr.975.fix = corr.975.fix + 1
  }
  if (real.LR.pen.final[s] <= quantile(boot.LR.fix.pen.final[[s]], probs = 0.975)){
    corr.975.fix.pen = corr.975.fix.pen + 1
  }
}
corr.975.fix
corr.975.fix.pen
