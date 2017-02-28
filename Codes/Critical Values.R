# -------------------------------------------------------------------------------
# Critical values calibration (according to Chen and Niu (2014))
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Simulation setup
n.obs      = 1200            # no of observations
n.par      = 12              # no of parameters
n.sim      = 1000            # no of simulations
seed1      = 20150206        # seed simulation X
seed2      = 20150602        # seed simulation epsilon
M          = 50              # increment of observations between successive subsamples
K          = 24              # number of subsamples
sig        = 1               # st. dev. of the error term

# True beta coefficients (homogeneous for all t = 1, ..., 1000)
tmp1  = c(1, 1.5, 1, 1, 2, -3, -1.5, 1, 2, 5, 3, 1)
tmp2  = rep(0, n.par - length(tmp1))
b     = c(tmp1, tmp2)

# Simulation of the design matrix
mu    = rep(0, n.par)
r     = 0.5
Sigma = matrix(0, nrow = n.par, ncol = n.par)

for (i in 1:n.par) {
  for (j in 1:n.par) {
    if (i == j){
      Sigma[i, j] = 1
    }else {
      Sigma[i, j] = r^abs(i - j)
    }
  }
}

X = list()
set.seed(seed1)
for (i in 1:n.sim){
  X[[i]] = mvrnorm(n = n.obs, mu, Sigma)
}  

# Simulation of the error term for t = 1, ..., n.obs
eps  = list()
set.seed(seed2)
for (i in 1:n.sim){
  eps[[i]] = rnorm(n.obs, mean = 0, sd = sig)
} 

# Computation of Y for t = 1, ..., 1000
Y    = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in 1:n.obs){
    Y.tmp = c(Y.tmp, b %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y[[i]] = Y.tmp
}

s = 1

# Find OLS (also MLE) for every subsample k = 1, ..., K

betas = list()
for (s in 1:n.sim){
  betas.tmp = numeric(0)
  for (k in 1:K){
    X.tmp = X[[s]][1:(k*M),]
    Y.tmp = Y[[s]][1:(k*M)]
    beta.fit  = solve(t(X.tmp) %*% X.tmp) %*% t(X.tmp) %*% Y.tmp
    betas.tmp = cbind(betas.tmp, beta.fit)
  }
  betas[[s]] = betas.tmp 
}

test.stat = function(s, k){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  loglik1   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
               %*% (Y.tmp - X.tmp %*% beta.tmp1)))

  beta.tmp2 = as.matrix(betas[[s]][, (k - 1)])
  loglik2   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
               %*% (Y.tmp - X.tmp %*% beta.tmp2)))
     
  test.stat = sqrt(abs(loglik1 - loglik2))
  test.stat
}

t.stat = numeric(0)
for (s in 1:n.sim){
  tstat.tmp = 0
  for (k in 2:K){
    tstat.tmp = c(tstat.tmp, test.stat(s, k))
  }
  t.stat = cbind(t.stat, tstat.tmp)
}

est.error = function(s, k){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  beta.true = solve(t(X[[s]]) %*% X[[s]]) %*% t(X[[s]]) %*% Y[[s]]
  loglik1   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                  %*% (Y.tmp - X.tmp %*% beta.tmp1)))

  loglik2   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.true)
                  %*% (Y.tmp - X.tmp %*% beta.true)))
  
  est.error = sqrt(abs(loglik1 - loglik2))
  est.error
}

est.err = numeric(0)
for (s in 1:n.sim){
  err.tmp = numeric(0)
  for (k in 1:K){
    err.tmp = c(err.tmp, est.error(s, k))
  }
  est.err = cbind(est.err, err.tmp)
}

err.exp = apply(est.err, 1, mean)

dist.fct = function(s, l, k){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  loglik1   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                  %*% (Y.tmp - X.tmp %*% beta.tmp1)))
  
  beta.tmp2 = as.matrix(betas[[s]][, l])
  loglik2   = (-(k * M)/2 * (2 * pi * sig) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
                  %*% (Y.tmp - X.tmp %*% beta.tmp2)))
  
  test.stat = sqrt(abs(loglik1 - loglik2))
  test.stat
}

# Set starting values of expected stochastic distances and critical values
exp.dist = apply(t.stat, 1, mean)
zeta     = c(rep(0, (K - 1)))
incr     = 10^(-3)
dist     = t.stat

k = 1
# Loop for subsamples k = 2, ..., K
for (k in 1:(K - 1)){
  while (exp.dist[k] > err.exp[k]) {
    zeta.tmp = zeta[k] + incr
    dist.tmp = numeric(0)
    for (s in 1:n.sim){
      dist.tmp[s] = ifelse(t.stat[k, s] <= zeta.tmp, 0, dist[k, s])
    }
    exp.dist[k] = mean(dist.tmp)
    zeta[k] = zeta.tmp
  }
  for (s in 1:n.sim){
    for (l in (k + 1):(K - 1)) {
      dist[l, s] = ifelse(dist.tmp[s] == 0, dist[l, s], dist.fct(s, k, l))
    }
  }
}








