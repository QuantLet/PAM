# -------------------------------------------------------------------------------
# Excess bond risk premia modelling with PAM
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx", "lars", 
              "scales", "tilting", "VGAM")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("OnestepSCADpen_withBIC_norm.r")

# Load data
tmpdata     = read.csv("BRP_data.csv",sep=",") # Bond Risk Premium data
data        = sapply(subset(tmpdata, select = c(2:(dim(tmpdata)[2]))), as.numeric)
dates       = as.Date(as.character(as.POSIXct(as.character(tmpdata[, 1]), format = "%Y%m%d")))
covariates0 = as.matrix(data[, 5:dim(data)[2]])
n.monthly   = dim(covariates0)[1]
n.weekly    = n.monthly * 4

seq.tmp = 1:n.monthly
covariates1 = matrix(0, nrow = n.weekly, ncol = dim(covariates0)[2])
for (icov in 1:dim(covariates0)[2]){
  covariates1[, icov] = spline(seq.tmp, covariates0[, icov], n = n.weekly)$y
}
covariates2 = covariates1[, -c(12,13,19,20,25)]
covariates  = covariates2[, -c(6,7,8,9,10)]
covariates  = covariates[, c(1:5)]
BRP2        = spline(seq.tmp, data[, 1], n = n.weekly)$y
BRP3        = spline(seq.tmp, data[, 2], n = n.weekly)$y
BRP4        = spline(seq.tmp, data[, 3], n = n.weekly)$y
BRP5        = spline(seq.tmp, data[, 4], n = n.weekly)$y

# Define pre-specified parameters
n.par       = dim(covariates)[2]           # no of parameters
M           = 48 * 4                       # increment of observations 
K           = dim(covariates)[1]%/%M       # number of subsamples
a           = 3.7                          # recommended value of parameter a for SCAD
sd.eps      = 1 
n.boot      = 1000

# Simulating mutlipliers from a bounded distribution
runif.mod = function(n){
  unif1 = runif(n, min = 0, max = 12)
  unif.a = unif1[which(unif1 <= 1)]
  index.a  = which(unif1 <= 1)
  unif.b = unif1[which((1 < unif1) & (unif1 <= 4))]
  index.b  = which((1 < unif1) & (unif1 <= 4))
  unif.c = (unif1[which(unif1 > 4)] - 4)/8
  index.c  = which(unif1 > 4)
  unif2   = cbind(c(index.a, index.b, index.c), c(unif.a, unif.b, unif.c))
  unif.fin = unif2[order(unif2[, 1], decreasing = FALSE), 2]
  unif.fin
}

# Function to compute bootstrapped penalized log-likelihood function
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


Y           = BRP2 
X           = covariates

b.pts       = numeric(0)
k           = 1
k0          = 0
k1          = 0
k2          = 1
beta        = list()
beta.0      = list()
lambda      = list()
bic         = list()

# 1.step: Fit the model for assumed homogeneous interval I_t^(1)
X.tmp.L      = X[((k1 * M) + 1):(k2 * M), ]
Y.tmp.L      = Y[((k1 * M) + 1):(k2 * M)]
n.obs.tmp.L  = length(Y.tmp.L)
object.tmp.L = Adapt.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
beta[[k]]    = object.tmp.L$beta
lambda[[k]]  = object.tmp.L$lambda
beta.0[[k]]  = object.tmp.L$beta.0
bic[[k]]     = object.tmp.L$bic
boot.L = numeric(0)
boot.R = numeric(0)
boot.T = numeric(0)

Sys.time()
while (k < 10){
  # 1.step: Fit the model for assumed homogeneous interval I_t^(k)
  X.tmp.L       = X[((k1 * M) + 1):(k2 * M), ]
  Y.tmp.L       = Y[((k1 * M) + 1):(k2 * M)]
  n.obs.tmp.L   = length(Y.tmp.L)
  object.tmp.L  = Adapt.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
  beta.tmp.L    = object.tmp.L$beta
  lambda.tmp.L  = object.tmp.L$lambda
  beta0.tmp.L   = object.tmp.L$beta.0
  bic.tmp.L     = object.tmp.L$bic

  # 2. a) step: Fit the model over the next interval I_t^(k + 1) - I_t^(k)
  X.tmp.R      = X[((k2 * M) + 1):((k2 + 1) * M), ]
  Y.tmp.R      = Y[((k2 * M) + 1):((k2 + 1) * M)]
  n.obs.tmp.R  = length(Y.tmp.R)
  object.tmp.R = Adapt.SCAD(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R)
  beta.tmp.R   = object.tmp.R$beta
  lambda.tmp.R = object.tmp.R$lambda
  beta0.tmp.R  = object.tmp.R$beta.0
  bic.tmp.R     = object.tmp.R$bic
  
  # 2. b) step: Fit the model over the whole interval I_t^(k + 1)
  X.tmp.T      = rbind(X.tmp.L, X.tmp.R)
  Y.tmp.T      = c(Y.tmp.L, Y.tmp.R)
  n.obs.tmp.T  = n.obs.tmp.L + n.obs.tmp.R
  object.tmp.T = Adapt.SCAD(X.tmp.T, Y.tmp.T, a, n.obs.tmp.T)
  beta.tmp.T   = object.tmp.T$beta
  lambda.tmp.T = object.tmp.T$lambda
  beta0.tmp.T  = object.tmp.T$beta.0
  bic.tmp.T     = object.tmp.T$bic
 
  # 3.step: Evaluate the test statistic
  # Real penalized likelihood ratio
  lik.L.pen = -loglik.pen(beta0.tmp.L, beta.tmp.L, lambda.tmp.L, X.tmp.L, Y.tmp.L)
  lik.R.pen = -loglik.pen(beta0.tmp.R, beta.tmp.R, lambda.tmp.R, X.tmp.R, Y.tmp.R)
  lik.T.pen = -loglik.pen(beta0.tmp.T, beta.tmp.T, lambda.tmp.T, X.tmp.T, Y.tmp.T)
  
  # Real unpenalized likelihood ratio
  lik.L     = -loglik(beta.tmp.L, lambda.tmp.L, X.tmp.L, Y.tmp.L)
  lik.R     = -loglik(beta.tmp.R, lambda.tmp.R, X.tmp.R, Y.tmp.R)
  lik.T     = -loglik(beta.tmp.T, lambda.tmp.T, X.tmp.T, Y.tmp.T)
  
  t.stat.pen = ((n.obs.tmp.L/n.obs.tmp.T) * lik.L.pen 
                + (n.obs.tmp.R/n.obs.tmp.T) * lik.R.pen - lik.T.pen)
  
  t.stat     = ((n.obs.tmp.L/n.obs.tmp.T) * lik.L 
                + (n.obs.tmp.R/n.obs.tmp.T) * lik.R - lik.T)
  
  # 4.step: Find quantiles of LR with use of MB
  
  lik.ratio.fix.pen   = numeric(0)
  
  # Simulation of multipliers u_i (i = 1, ..., n.obs.tmp.T)
  set.seed      = 20170424 * k * 10
  # multipliers   = matrix(runif.mod(n.boot * n.obs.tmp.T),
  #                        ncol = n.obs.tmp.T, nrow = n.boot)
  
  # multipliers = matrix(rexp(n.boot * n.obs.tmp.T, rate = 1), # exponential
  #                      ncol = n.obs.tmp.T, nrow = n.boot)
  
  multipliers = matrix(rpois(n.boot * n.obs.tmp.T, lambda = 1), # poisson
                       ncol = n.obs.tmp.T, nrow = n.boot)
  
  multipliers.L = multipliers[, 1:n.obs.tmp.L]
  multipliers.R = multipliers[, (n.obs.tmp.L + 1):n.obs.tmp.T]
  
  for (l in 1:(n.boot)){
    # Multiplier bootstrap for left-hand interval I_t^(k)
    multipl.L          = multipliers.L[l, ]
    
    boot.object.fix.L  = Adapt.SCAD.MB(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L, 
                                       as.numeric(lambda.tmp.L), multipl.L)
    
    boot.beta.fix.L    = boot.object.fix.L$beta
    boot.beta0.fix.L   = boot.object.fix.L$beta.0
    
    lik.boot.fix.L.pen = -loglik.pen.MB(boot.beta0.fix.L, boot.beta.fix.L, as.numeric(lambda.tmp.L), 
                                        X.tmp.L, Y.tmp.L, multipl.L)
    
    # Multiplier bootstrap for right-hand interval I_t^(k + 1) - I_t^(k)
    multipl.R         = multipliers.R[l, ]
    
    boot.object.fix.R = Adapt.SCAD.MB(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R, 
                                      as.numeric(lambda.tmp.R), multipl.R)
    
    boot.beta.fix.R   = boot.object.fix.R$beta
    boot.beta0.fix.R  = boot.object.fix.R$beta.0
    
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
    
    lik.boot.fix.T.pen = -loglik.pen.MB(boot.beta0.fix.T, boot.beta.fix.T, boot.lambda.fix.T, 
                                        X.tmp.T, Y.boot.T, multipl)
    
    lik.ratio.fix.pen  = c(lik.ratio.fix.pen, ((n.obs.tmp.L/n.obs.tmp.T) * lik.boot.fix.L.pen + 
                                                 (n.obs.tmp.R/n.obs.tmp.T) *lik.boot.fix.R.pen - 
                                                 lik.boot.fix.T.pen))
    print(l)
  }
  q50.med     = quantile(lik.ratio.fix.pen, probs = 0.5)
  q90         = quantile(lik.ratio.fix.pen, probs = 0.9)
  q95         = quantile(lik.ratio.fix.pen, probs = 0.95)
  q97.5       = quantile(lik.ratio.fix.pen, probs = 0.975)
  
  k1 = k1 + 1
  k2 = k2 + 1
  k  = k + 1
 
  # 5.step: Test for homogeneity and evaluate current beta estimator
  if (t.stat.pen <= q90){
    X.hom      = X[((k0 * M) + 1):(k2 * M), ]
    Y.hom      = Y[((k0 * M) + 1):(k2 * M)]
    n.obs.hom  = length(Y.hom)
    object.hom = Adapt.SCAD(X.hom, Y.hom, a, n.obs.hom)
    for (k3 in (k0 + 1):k){
      beta[[k3]]   = object.hom$beta
      lambda[[k3]] = object.hom$lambda
    }
  } else {
    b.pts       = c(b.pts, (((k - 1) * M) + 1))
    k0          = k2 - 1
    beta[[k]]   = beta.tmp.R
    lambda[[k]] = lambda.tmp.R
  } 
  print(k)
}
Sys.time()

Y.fitted = numeric(0)
# Plot Y.fitted vs real data
for (k in 1:K){
  n.obs   = M 
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X[(((k - 1) * M) + 1):(k * M), ] )/n.obs
  X.cent  = scale(X[(((k - 1) * M) + 1):(k * M), ], X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y[(((k - 1) * M) + 1):(k * M)])/n.obs
  Y.norm  = scale(Y[(((k - 1) * M) + 1):(k * M)], Y.mean, FALSE)
  
  Y.fitted = c(Y.fitted, X.cent %*% beta[[k]] + Y.mean)
}
length(Y.fitted)

plot(Y[1:length(Y.fitted)], type = "l", ylim = c(-0.05, 0.05))
lines(Y.fitted, col = "red")

