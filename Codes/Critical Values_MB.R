# -------------------------------------------------------------------------------
# Critical values calibration using Multiplier Bootstrap and SCAD
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
n.obs       = 1200            # no of observations
n.par       = 20              # no of parameters
n.sim       = 1               # no of simulations
seed1       = 20150206        # seed simulation X
seed2       = 20150602        # seed simulation epsilon
M           = 50              # increment of observations between successive subsamples
K           = 24              # number of subsamples
sd.eps      = 1               # st. dev. of the error term
max.steps   = 50                # Max no of iterations in Adapt.SCAD
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
n.cores = detectCores()        # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of simulations for every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h          = n.sim %% n.cores
if (h != 0){
  n.simcores[1:h] = n.simcores[1:h] + 1
}

# Function to find betas
find.betas = function(l){
  betas  = list()
  wghts  = list()
  pens   = list()
  lambda = list()
  s1 = ifelse(l == 1, 1, sum(n.simcores[1:(l - 1)],1))
  s2 = sum(n.simcores[1:l])
  
  for (s in s1:s2){
    ad.beta    = numeric(0)
    ad.penalty = numeric(0)
    ad.weights = numeric(0)
    ad.lambda  = numeri(0)
    for (k in 1:K){
      X.tmp          = X[[s]][1:(k * M),]
      Y.tmp          = Y[[s]][1:(k * M)]
      adapt.estim    = Adapt.SCAD(X.tmp, Y.tmp, a, (k * M), max.steps)
      beta.fit       = adapt.estim$beta
      pen.fit        = adapt.estim$penalty
      w.fit          = adapt.estim$weights
      l.fit          = adapt.estim$lambda
      
      ad.beta        = cbind(ad.beta, beta.fit)
      ad.penalty     = cbind(ad.penalty, pen.fit)
      ad.weights     = cbind(ad.weights, w.fit)
      ad.lambda      = cbind(ad.lambda, l.fit)
    }
    stmp             = s - s1 + 1
    betas[[stmp]]    = ad.beta
    pens[[stmp]]     = ad.penalty
    wghts[[stmp]]    = ad.weights
    lambdass[[stmp]] = ad.lambda
  }
  values          = list(betas, pens, wghts, lambdas)
  names(values)   = c("betas", "pens", "wghts", "lambdas") 
  return(values)
}

# Collect results from all used cores
res.betas = function(n.cores, input){
  betas   = list()
  pens    = list()
  wghts   = list()
  lambdas = list()
  for (i in 1:n.cores){
    betas   = c(betas, input[[i]]$betas)
    pens    = c(pens, input[[i]]$pens)
    wghts   = c(wghts, input[[i]]$wghts)
    lambdas = c(lambdas, input[[i]]$lambdas)
  }
  values         = list(betas, pens, wghts, lambdas)
  names(values)  = c("betas", "pens", "wghts", "lambdas") 
  return(values)
}

# Compute estimated betas
Sys.time()
betas.tmp   = foreach(i = 1:n.cores, .packages = "glmnet") %dopar% find.betas(i)   # Parallel computing
Sys.time()
final.betas = res.betas(n.cores, betas.tmp) 

betas   = final.betas$betas
pens    = final.betas$pens
wghts   = final.betas$wghts
lambdas = final.betas$lambdas

# Take 1 simulated scenario and create lists of betas, X's and Y's
s = 1
beta.new  = list()
X.tmp.new = list()
Y.tmp.new = list()
for (k in 1:K){
 nonzero        = which(betas[[s]][, k] != 0)
 beta.new[[k]]  = as.matrix(betas[[s]][nonzero, k])
 X.tmp.new[[k]] = X[[s]][1:(k * M), nonzero]
 Y.tmp.new[[k]] = Y[[s]][1:(k * M)]
}

# For each k, simulate multipliers from N(1,1)
n.boot      = 10000
multipliers = list()
for (k in 1:K){
  set.seed    = 20170424
  n.obs.multi = k * M
  multipliers[[k]] = matrix(rnorm(n.boot * n.obs.multi, sd = 1, mean = 1), 
                            ncol = n.obs.multi, nrow = n.boot)
}

mutli.boot = function(beta, X, Y, k, l){
  loglik.boot.fitted = (-(k * M)/2 * log(2 * pi * sd.eps^2) * sum(multipliers[[k]][l, ]) - (1/(2 * sd.eps^2) 
                        * sum((Y - X %*% beta)^2 * multipliers[[k]][l, ])))
  
}

# Delete the zero parameters and corresponding regressors
loglik = function(beta){
  X = X.tmp.new[[k]]
  Y = Y.tmp.new[[k]]
  loglik1 = -(-(k * M)/2 * log(2 * pi * sd.eps^2) * sum(multipliers[[k]][l, ]) - (1/(2 * sd.eps^2) 
                                                                            * sum((Y - X %*% beta)^2 * multipliers[[k]][l, ])))
  loglik1
}

q90       = numeric(0)
q95       = numeric(0)
q97.5     = numeric(0)
lik.ratio = numeric(0)
Sys.time()
for (k in 1:K){
  for (l in 1:n.boot){
    boot.beta     = optim(beta.new[[k]], loglik)$par
    lik.boot.beta = -loglik(boot.beta)
    lik.beta      = -loglik(beta.new[[k]])
    lik.ratio     = c(lik.ratio, sqrt(2 * (lik.boot.beta - lik.beta)))
  }
  q90[k]   = quantile(lik.ratio, probs = 0.9)
  q95[k]   = quantile(lik.ratio, probs = 0.95)
  q97.5[k] = quantile(lik.ratio, probs = 0.975)
  print(k)
}
Sys.time()

zeta90   = q90    # C.v.'s defined by 90 % quantile
zeta95   = q95    # C.v.'s defined by 95 % quantile
zeta97.5 = q97.5  # C.v.'s defined by 97.5 % quantile

