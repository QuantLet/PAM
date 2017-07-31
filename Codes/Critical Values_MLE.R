# -------------------------------------------------------------------------------
# Critical values calibration (according to Chen and Niu (2014)) using MLE
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "lars")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Simulation setup
n.obs      = 1000            # no of observations
n.par      = 10              # no of parameters
n.sim      = 1000            # no of simulations
seed1      = 20150206        # seed simulation X
seed2      = 20150602        # seed simulation epsilon
M          = 50              # increment of observations between successive subsamples
K          = 20              # number of subsamples
sd.eps     = 1               # st. dev. of the error term
risk.bound = gamma(1/2)      # define riskbound

# True beta coefficients (homogeneous for all t = 1, ..., n.obs)
tmp1  = c(1, 1, 1, 1, 1)
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

# Find OLS (also MLE) for every subsample k = 1, ..., K
betas = list()
for (s in 1:n.sim){
  betas.tmp = numeric(0)
  for (k in 1:K){
    X.tmp = X[[s]][1:(k * M),]
    Y.tmp = Y[[s]][1:(k * M)]
    beta.fit  = solve(t(X.tmp) %*% X.tmp) %*% t(X.tmp) %*% Y.tmp
    betas.tmp = cbind(betas.tmp, beta.fit)
  }
  betas[[s]] = betas.tmp 
}

# Function for computing test statistic
test.stat = function(s, k, l){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  loglik1   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
               %*% (Y.tmp - X.tmp %*% beta.tmp1)))

  beta.tmp2 = as.matrix(betas[[s]][, l])
  loglik2   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
               %*% (Y.tmp - X.tmp %*% beta.tmp2)))
     
  test.stat = sqrt(2*abs(loglik1 - loglik2))
  # test.stat = (abs(loglik1 - loglik2))^(1/3) # Different norm
  test.stat
}

# Function for computing estimation error
est.error = function(s, k){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  beta.true = solve(t(X[[s]]) %*% X[[s]]) %*% t(X[[s]]) %*% Y[[s]]
  loglik1   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                  %*% (Y.tmp - X.tmp %*% beta.tmp1)))

  loglik2   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.true)
                  %*% (Y.tmp - X.tmp %*% beta.true)))
  
  est.error = sqrt(2*abs(loglik1 - loglik2))
  # est.error = (abs(loglik1 - loglik2))^(1/3) # Different norm
  est.error
}

# Compute estimation error
est.err = numeric(0)
for (s in 1:n.sim){
  err.tmp = numeric(0)
  for (k in 1:K){
    err.tmp = c(err.tmp, est.error(s, k))
  }
  est.err = cbind(est.err, err.tmp)
}

# Expected estimation error
err.exp = apply(est.err, 1, mean)

# Function for computing stochastic distance
dist.fct = function(s, k, l){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  loglik1   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                  %*% (Y.tmp - X.tmp %*% beta.tmp1)))
  
  beta.tmp2 = as.matrix(betas[[s]][, l])
  loglik2   = (- (1 / (2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
                  %*% (Y.tmp - X.tmp %*% beta.tmp2)))
  
  stoch.dist = sqrt(2*abs(loglik1 - loglik2))
  # stoch.dist = (abs(loglik1 - loglik2))^(1/3) # Different norm
  stoch.dist
}

# Function for calibration of critical values based on the risk bound defined
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
  }
  zeta
}

# Function for calibration of critical values based on the expected estimation error
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
  }
  zeta
}

# Find critical values based on Chen&Niu(2014)
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
  zeta[m] = cv.calib.ee(zeta[m], incr)
}
Sys.time()


# C.v.'s 
zeta.esterr.mle = c(Inf, 6.84, 5.18, 3.61, 3.27, 2.86, 2.79, 2.17, 2.17,
                    2.04, 2.09, 1.78, 1.74, 1.71, 1.68, 1.62, 1.50, 1.34, 1.63, 1.81)

# C.v.'s based on different norm
zeta.esterr.mle.13 = c(Inf, 3.90, 2.69, 2.40, 2.07, 1.91, 1.75, 1.60, 1.57, 1.53, 1.44, 1.41, 1.36,
                       1.36, 1.33, 1.19, 1.23, 1.23, 1.13, 1.11, 1.15, 1.12, 1.08, 1.14)
dim(as.matrix(c(Inf, 3.90, 2.69, 2.40, 2.07, 1.91, 1.75, 1.60, 1.57, 1.53, 1.44, 1.41, 1.36,
  1.36, 1.33, 1.19, 1.23, 1.23, 1.13, 1.11, 1.15, 1.12, 1.08, 1.14)))[2]
