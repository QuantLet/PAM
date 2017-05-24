# -------------------------------------------------------------------------------
# Critical values calibration (according to Chen and Niu (2014)) using Adapt.SCAD
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("Adapt.SCAD_withCV.r")

# Simulation setup
n.obs      = 1200             # no of observations
n.par      = 20               # no of parameters
n.sim      = 1000             # no of simulations
seed1      = 20150206         # seed simulation X
seed2      = 20150602         # seed simulation epsilon
M          = 50               # increment of observations between successive subsamples
K          = 24               # number of subsamples
sd.eps     = 1                # st. dev. of the error term
risk.bound = gamma(1/2)       # Define risk bound 
max.steps   = 50              # Max no of iterations in Adapt.SCAD
lambda.grid = seq(1 * n.obs^(-1/3), 30 * n.obs^(-1/3), 1 * n.obs^(-1/3))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD

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

# Computation of Y for t = 1, ..., n.obs
Y    = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in 1:n.obs){
    Y.tmp = c(Y.tmp, b %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y[[i]] = Y.tmp
}
plot(Y[[i]])

# Initiate cluster for parallel computing
n.cores = detectCores()         # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of simulations for every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h          = n.sim %% n.cores
if (h != 0){
  n.simcores[1:h] = n.simcores[1:h] + 1
}

# Funtion for estimation of betas with help of adaptive SCAD method
find.betas = function(l){
  betas = list()
  wghts = list()
  pens  = list()
  s1 = ifelse(l == 1, 1, sum(n.simcores[1:(l - 1)],1))
  s2 = sum(n.simcores[1:l])
  
  for (s in s1:s2){
    ad.beta    = numeric(0)
    ad.penalty = numeric(0)
    ad.weights = numeric(0)
    for (k in 1:K){
      X.tmp       = X[[s]][1:(k * M),]
      Y.tmp       = Y[[s]][1:(k * M)]
      adapt.estim = Adapt.SCAD(X.tmp, Y.tmp, a, (k * M), max.steps)
      beta.fit    = adapt.estim$beta
      pen.fit     = adapt.estim$penalty
      w.fit       = adapt.estim$weights
      
      ad.beta     = cbind(ad.beta, beta.fit)
      ad.penalty  = cbind(ad.penalty, pen.fit)
      ad.weights  = cbind(ad.weights, w.fit)
    }
    stmp          = s - s1 + 1
    betas[[stmp]] = ad.beta
    pens[[stmp]]  = ad.penalty
    wghts[[stmp]] = ad.weights
  }
  values          = list(betas, pens, wghts)
  names(values)   = c("betas", "pens", "wghts") 
  return(values)
}

# Collect results over all used cores
res.betas = function(n.cores, input){
  betas = list()
  pens  = list()
  wghts = list()
  for (i in 1:n.cores){
    betas = c(betas, input[[i]]$betas)
    pens  = c(pens, input[[i]]$pens)
    wghts = c(wghts, input[[i]]$wghts)
  }
  values         = list(betas, pens, wghts)
  names(values)  = c("betas", "pens", "wghts") 
  return(values)
}

# Compute the estimated values
Sys.time()
betas.tmp   = foreach(i = 1:n.cores, .packages = "glmnet") %dopar% find.betas(i)   # Parallel computing
Sys.time()
final.betas = res.betas(n.cores, betas.tmp) 

betas = final.betas$betas
pens  = final.betas$pens
wghts = final.betas$wghts
betas.scad = betas
pens.scad  = pens
wghts.scad = wghts

# Test statistic
test.stat = function(s, k, l){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  pen.tmp1  = as.matrix(pens[[s]][, k])
  wght.tmp1 = as.matrix(wghts[[s]][, k])
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp1)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp1))
               - ((k * M) * t(pen.tmp1)
                  %*% abs(wght.tmp1)))
  
  beta.tmp2 = as.matrix(betas[[s]][, l])
  pen.tmp2  = as.matrix(pens[[s]][, l])
  wght.tmp2 = as.matrix(wghts[[s]][, l])
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp2)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp2))
               - ((k * M) * t(pen.tmp2)
                  %*% abs(wght.tmp2)))
  
  test.stat = sqrt(abs(loglik1 - loglik2))
  test.stat
}

# Estimation error
est.error = function(s, k){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  pen.tmp1  = as.matrix(pens[[s]][, k])
  wght.tmp1 = as.matrix(wghts[[s]][, k])
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp1)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp1))
               - ((k * M) * t(pen.tmp1)
                  %*% abs(wght.tmp1)))
  
  beta.tmp2 = as.matrix(betas[[s]][, K])
  pen.tmp2  = as.matrix(pens[[s]][, K])
  wght.tmp2 = as.matrix(wghts[[s]][, K])
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp2)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp2))
               - ((k * M) * t(pen.tmp2)
                  %*% abs(wght.tmp2)))
  
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

# Expected estimation error
err.exp = apply(est.err, 1, mean)

# Function for computing stochastic distance
dist.fct = function(s, k, l){
  X.tmp     = X[[s]][1:(k * M),]
  Y.tmp     = Y[[s]][1:(k * M)]
  beta.tmp1 = as.matrix(betas[[s]][, k])
  pen.tmp1  = as.matrix(pens[[s]][, k])
  wght.tmp1 = as.matrix(wghts[[s]][, k])
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp1)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp1))
               - ((k * M) * t(pen.tmp1)
                  %*% abs(wght.tmp1)))
  
  beta.tmp2 = as.matrix(betas[[s]][, l])
  pen.tmp2  = as.matrix(pens[[s]][, l])
  wght.tmp2 = as.matrix(wghts[[s]][, l])
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp2)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp2))
               - ((k * M) * t(pen.tmp2)
                  %*% abs(wght.tmp2)))
  
  stoch.dist = sqrt(abs(loglik1 - loglik2))
  stoch.dist
}

# Calibration of C.v.'s based on the defined risk bound
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

# Calibration of C.v.'s based on the expected error
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
        simchp = c(simchp, s)
      }
    }  
    exp.dist  = apply(as.matrix(dist[[m]][, (m : K)]), 2, mean)
    dist.test = (exp.dist <= err.exp[(m : K)])
  }
  zeta
}

# Find critical values based on Chen&Niu(2014)
dist      = list()
dist[[1]] = matrix(0, ncol = K, nrow = n.sim)
zeta      = c(rep(Inf, K))
incr      = 10^(-2)

Sys.time()
for (m in (2 : K)){
  zeta[m]  = 0
  dist.test = FALSE
  
  # Matrix of stochastic distances D_t^(k)
  dist[[m]]   = matrix(0, ncol = K, nrow = n.sim)
  dist.non0   = matrix(0, ncol = K, nrow = n.sim)
  for (s in 1:n.sim){
    for (k in (m : K)){
      dist.non0[s, k] = dist.fct(s, k, (m - 1))
    }
  }
  zeta[m] = cv.calib.ee(zeta[m], incr)
}
Sys.time()

# Close cluster
stopCluster(cl)

