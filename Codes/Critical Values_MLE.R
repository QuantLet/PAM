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
  eps[[i]] = rnorm(n.obs, mean = 0, sd = sd.eps)
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
plot(Y[[i]])

risk.bound = gamma(1/2)

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

test.stat = function(s, k, l){
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

dist.fct = function(s, k, l){
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
  
  stoch.dist = sqrt(abs(loglik1 - loglik2))
  stoch.dist
}

cv.calib = function(zeta, incr){
  while (is.element(FALSE, dist.test)){
    zeta = zeta + incr
    for (s in 1:n.sim){
      if (dist[[(m - 1)]][s, (m - 1)] == 0 && test.stat(s, m, (m - 1)) <= zeta) {
        dist[[m]][s, (m : K)] = 0
      } else if (dist[[(m - 1)]][s, (m - 1)] == 0 && test.stat(s, m, (m - 1)) > zeta) {
        dist[[m]][s, (m : K)] = dist.non0[s, (m : K)]
      } else {
        dist[[m]][s, ] = dist[[(m - 1)]][s, ]
      }
    }  
    exp.dist  = apply(as.matrix(dist[[m]][, (m : K)]), 2, mean)
    dist.test = (exp.dist <= ((m - 1)/(K - 1) * risk.bound))
    # if (all(dist.test)) break
    # zeta = zeta + incr
  }
  zeta
}


cv.calib.ee = function(zeta, incr){
  while (is.element(FALSE, dist.test)){
    zeta = zeta + incr
    for (s in 1:n.sim){
      if (dist[[(m - 1)]][s, (m - 1)] == 0 && test.stat(s, m, (m - 1)) <= zeta) {
        dist[[m]][s, (m : K)] = 0
      } else if (dist[[(m - 1)]][s, (m - 1)] == 0 && test.stat(s, m, (m - 1)) > zeta) {
        dist[[m]][s, (m : K)] = dist.non0[s, (m : K)]
      } else {
        dist[[m]][s, ] = dist[[(m - 1)]][s, ]
      }
    }  
    exp.dist  = apply(as.matrix(dist[[m]][, (m : K)]), 2, mean)
    dist.test = (exp.dist <= err.exp[(m : K)])
    # if (all(dist.test)) break
    # zeta = zeta + incr
  }
  zeta
}




# Hladam zeta_m, cize l = 1 (ked zamietnem v 2. kroku homo, tak moj odhad je vzdy theta.hat_1)

dist      = list()
dist[[1]] = matrix(0, ncol = 24, nrow = 1000)
zeta      = c(rep(Inf, K))
incr      = 10^(-2)

Sys.time()
for (m in (2 : K)){
  zeta[m]  = 0
  dist.test = FALSE
  
  # Matrix of stochastic distances D_t^(k)
  dist[[m]]   = matrix(0, ncol = 24, nrow = 1000)
  dist.non0   = matrix(0, ncol = 24, nrow = 1000)
  for (s in 1:n.sim){
    for (k in (m : K)){
      dist.non0[s, k] = dist.fct(s, k, (m - 1))
    }
  }
  zeta[m] = cv.calib(zeta[m], incr)
}
Sys.time()

zeta.riskbound.mle = zeta
zeta.esterr.mle    = zeta
