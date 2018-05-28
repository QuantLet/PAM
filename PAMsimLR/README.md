[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **PAMsimLR** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: PAMsimLR

Published in: submitted to N/A 

Description: ‘Performs the Penalized Adaptive Method (PAM), a combination of propagation-separation approach and a SCAD penalty, to evaluate empirical coverage probability of bootstrapped confidence regions of penalized likelihood ratio.’

Keywords: ‘linear model, regression, SCAD penalty, bic, simulations, bootstrap, quantile’

See also: ‘PAMsimCP, PAMCocPia, PAMinsam, PAMoutsam’

Author: Lenka Zboňáková

Submitted:  23 May 2018 by Lenka Zboňáková

Input: 
- n.sim   : Number of simulated scenarios
- r       : Correlation parameter of the design matrix X
- sd.eps  : Standard deviation of the error term
- n.boot  : Number of bootstrap loops
- n.obs   : Number of observations
- q.par   : Number of nonzero parameters
- n.par   : Number of parameters
- mb.type : Distribution of multipliers (Bound, Exp or Pois)
```

### R Code
```r

# Clear all variables
rm(list = ls(all = TRUE))
graphics.off()

# Set directory
setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx",
              "lars", "scales", "tilting")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("PAMfunc.r")

# Initiate cluster for parallel computing
n.cores = detectCores()           # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

obs.sim = function(n.obs, q.par, n.par, n.sim, r, sd.eps){
  # True beta coefficients (homogeneous for all t = 1, ..., 1000)
  tmp1  = rep(1, q.par)
  tmp2  = rep(0, n.par - length(tmp1))
  b     = c(tmp1, tmp2)
  
  # Simulation of the design matrix
  mu    = rep(0, n.par)
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
  
  values = list(b, X, Y)
  names(values) = list("b", "X", "Y")
  return(values)
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

# Function returning real likelihood ratio and its bootstrapped counterpart 
# (whole distribution)
sim.LR.MBLR = function(h, mb.type = c("Bound", "Exp", "Pois")){
  K       = n.obs/2
  M       = n.obs/K
  s1      = ifelse(h == 1, 1, sum(n.simcores[1:(h - 1)],1))
  s2      = sum(n.simcores[1:h])
  real.LR = numeric(0)
  boot.LR = list()
  
  for (s in s1:s2){
    k           = 1
    k1          = 0
    k2          = 1
    beta        = list()
    lambda      = list()
    b.pts.tmp   = numeric(0)
    sim.no      = s - s1 + 1
   
    # Fit the model over the whole interval I_t^(k + 1)
    X.tmp.T      = X[[s]][((k1 * K) + 1):((k2 + 1) * K), ]
    Y.tmp.T      = Y[[s]][((k1 * K) + 1):((k2 + 1) * K)]
    n.obs.tmp.T  = length(Y.tmp.T)
    object.tmp.T = Onestep.SCAD(X.tmp.T, Y.tmp.T, a, n.obs.tmp.T)
    beta.tmp.T   = object.tmp.T$beta
    a.tmp.T      = object.tmp.T$a
    lambda.tmp.T = object.tmp.T$lambda
    beta0.tmp.T  = object.tmp.T$beta.0

    # Evaluate the real penalized likelihood over the whole interval I_t^(k + 1)
    lik.T        = loglik.pen(a.tmp.T, beta0.tmp.T, beta.tmp.T, lambda.tmp.T, X.tmp.T,
                              Y.tmp.T)
    
    # Find quantiles of LR with use of MB
    # Simulation of multipliers u_i (i = 1, ..., n.obs.tmp.T)
    set.seed     = 20170424 * s * 10
    
    if (mb.type == "Bound"){
      multipliers   = matrix(runif.mod(n.boot * n.obs.tmp.T),      
                             ncol = n.obs.tmp.T, nrow = n.boot) # Bounded distribution
    }
    if (mb.type == "Exp"){
      multipliers   = matrix(rexp(n.boot * n.obs.tmp.T, rate = 1), 
                             ncol = n.obs.tmp.T, nrow = n.boot) # Exp(1) distribution
    }
    if (mb.type == "Pois"){
      multipliers   = matrix(rpois(n.boot * n.obs.tmp.T, lambda = 1),
                             ncol = n.obs.tmp.T, nrow = n.boot) # Pois(1) distribution
    }
    
    MB.ratios = numeric(0)
    lik.sel   = numeric(0)
    
    # Compute real and bootstrapped likelihood for all positions of s 
    for (ind in 1:length(sel)){
      lik.ratio    = numeric(0)
      
      # Fit the model over the interval I_t^(k,s)
      X.sel.L      = X[[s]][((k1 * K) + 1):(k2 * K + sel[ind]), ]
      Y.sel.L      = Y[[s]][((k1 * K) + 1):(k2 * K + sel[ind])]
      n.obs.sel.L  = length(Y.sel.L)
      object.sel.L = Onestep.SCAD(X.sel.L, Y.sel.L, a, n.obs.sel.L)
      beta.sel.L   = object.sel.L$beta
      a.sel.L      = object.sel.L$a
      lambda.sel.L = object.sel.L$lambda
      beta0.sel.L  = object.sel.L$beta.0
      
      # Fit the model over the interval I_t^(k + 1,s) - I_t^(k,s)
      X.sel.R      = X[[s]][((k2 * K) + 1 + sel[ind]):((k2 + 1) * K), ]
      Y.sel.R      = Y[[s]][((k2 * K) + 1 + sel[ind]):((k2 + 1) * K)]
      n.obs.sel.R  = length(Y.sel.R)
      object.sel.R = Onestep.SCAD(X.sel.R, Y.sel.R, a, n.obs.sel.R)
      beta.sel.R   = object.sel.R$beta
      a.sel.R      = object.sel.R$a
      lambda.sel.R = object.sel.R$lambda
      beta0.sel.R  = object.sel.R$beta.0
      
      X.sel.T      = X.tmp.T
      Y.sel.T      = Y.tmp.T
      n.obs.sel.T  = n.obs.tmp.T
      
      # Evaluate the real penalized likelihood ratio over the whole interval I_t^(k + 1,s)
      lik.sel.L    = loglik.pen(a.sel.L, beta0.sel.L, beta.sel.L, lambda.sel.L, X.sel.L, 
                                   Y.sel.L)
      lik.sel.R    = loglik.pen(a.sel.R, beta0.sel.R, beta.sel.R, lambda.sel.R, X.sel.R, 
                                   Y.sel.R)
      lik.sel.T    = lik.T
      
      lik.sel      = c(lik.sel, ((n.obs.sel.L/n.obs.sel.T) * lik.sel.L 
                       + (n.obs.sel.R/n.obs.sel.T) * lik.sel.R - lik.sel.T))
      
      for (l in 1:(n.boot)){
        # Multiplier bootstrap for the left-hand interval I_t^(k, s)
        multipl.L     = multipliers[l, 1:n.obs.sel.L]
        
        boot.object.L = Onestep.SCAD.MB(X.sel.L, Y.sel.L, a, n.obs.sel.L, 
                                        as.numeric(lambda.sel.L), multipl.L)
        
        boot.beta.L   = boot.object.L$beta
        boot.a.L      = boot.object.L$a
        boot.a0.L     = boot.object.L$a.0
        boot.beta0.L  = boot.object.L$beta.0
        
        lik.boot.L    = loglik.pen.MB(boot.a.L, boot.beta0.L, boot.beta.L, 
                                      as.numeric(lambda.sel.L), X.sel.L, Y.sel.L, 
                                      multipl.L)
        
        # Multiplier bootstrap for the right-hand interval I_t^(k + 1,s) - I_t^(k,s)
        multipl.R     = multipliers[l, (n.obs.sel.L + 1):n.obs.sel.T]
        
        boot.object.R = Onestep.SCAD.MB(X.sel.R, Y.sel.R, a, n.obs.sel.R, 
                                          as.numeric(lambda.sel.R), multipl.R)
        
        boot.beta.R   = boot.object.R$beta
        boot.a.R      = boot.object.R$a
        boot.a0.R     = boot.object.R$a.0
        boot.beta0.R  = boot.object.R$beta.0
        
        lik.boot.R    = loglik.pen.MB(boot.a.R, boot.beta0.R, boot.beta.R, 
                                      as.numeric(lambda.sel.R), X.sel.R, Y.sel.R, 
                                      multipl.R)
        
        # Multiplier bootstrap for the whole interval I_t^(k + 1,s) (with shift)
        multipl        = multipliers[l, ]
        
        boot.object.T  = Onestep.SCAD.MB.shift(X.sel.L, Y.sel.L, X.sel.R, Y.sel.R, a, 
                                               lambda.sel.L, lambda.sel.R, multipl.L, 
                                               multipl.R, beta.sel.L, beta.sel.R, a.sel.L, 
                                               a.sel.R)
        
        boot.beta.T    = boot.object.T$beta
        boot.a.T       = boot.object.T$a
        boot.a0.T      = boot.object.T$a.0
        boot.lambda.T  = boot.object.T$lambda
        boot.beta0.T   = boot.object.T$beta.0
        
        Y.boot.T       = c(Y.sel.L, (Y.sel.R - rep((a.sel.R - a.sel.L), n.obs.sel.R) 
                                     - X.sel.R %*% (beta.sel.R - beta.sel.L)))
        
        lik.boot.T     = loglik.pen.MB(boot.a.T, boot.beta0.T, boot.beta.T, boot.lambda.T, 
                                       X.sel.T, Y.boot.T, multipl)
        
        lik.ratio      = c(lik.ratio, ((n.obs.sel.L/n.obs.sel.T) * lik.boot.L
                                       + (n.obs.sel.R/n.obs.sel.T) * lik.boot.R 
                                        - lik.boot.T))
      }
      
      MB.ratios        = cbind(MB.ratios, lik.ratio)
    }
    
    # Find maximum over the real and bootstrapped likelihood ratios
    real.LR[sim.no]    = max(lik.sel)
    boot.LR[[sim.no]]  = apply(MB.ratios, 1, max)
  }
  
  values = list(real.LR, boot.LR)
  names(values) = c("real.LR", "boot.LR")
  return(values)
}

# Evaluate number of events real.LR <= alpha-quantile of boot.LR (alpha - % quantile)
correct.cp = function(alphaquant, real.LR.input, boot.LR.input){
  corr.cp = 0
  for (s in 1:n.sim){
    if (real.LR.input[s] <= quantile(boot.LR.input[[s]], probs = alphaquant)){
      corr.cp = corr.cp + 1
    }
  }
  corr.cp
}

coll.res = function(input){
  
  real.LR     = numeric(0)
  boot.LR     = list()
  for (icor in 1:n.cores){
    real.LR = c(real.LR, input[[icor]]$real.LR)
    boot.LR = c(boot.LR, input[[icor]]$boot.LR)
  }
  
  corr.50  = correct.cp(0.5, real.LR, boot.LR)/n.sim * 100
  corr.75  = correct.cp(0.75, real.LR, boot.LR)/n.sim * 100
  corr.80  = correct.cp(0.8, real.LR, boot.LR)/n.sim * 100
  corr.85  = correct.cp(0.85, real.LR, boot.LR)/n.sim * 100
  corr.90  = correct.cp(0.9, real.LR, boot.LR)/n.sim * 100
  corr.95  = correct.cp(0.95, real.LR, boot.LR)/n.sim * 100
  corr.975 = correct.cp(0.975, real.LR, boot.LR)/n.sim * 100
  corr.99  = correct.cp(0.99, real.LR, boot.LR)/n.sim * 100
  
  values   = list(real.LR, boot.LR, corr.50, corr.75, 
                  corr.80, corr.85, corr.90, corr.95,
                  corr.975, corr.99)
  names(values) = c("real.LR", "boot.LR", "corr.50", "corr.75", 
                    "corr.80", "corr.85", "corr.90", "corr.95",
                    "corr.975", "corr.99")
  return(values) 
}

# Simulation setup
n.sim      = 1000           # Number of simulations
seed1      = 20171110       # Seed simulation X
seed2      = 20170602       # Seed simulation epsilon
r          = 0.5            # Correlation parameter of X
sd.eps     = 1              # Standard deviation of epsilon
a          = 3.7            # Recommended value of parameter a for SCAD
n.boot     = 1000           # Number of bootstrap loops
sel        = 0              # Position of s   
n.obs      = 100            # Number of observations
q.par      = 3              # Number of nonzero parameters
n.par      = 10             # Number of parameters
mb.type    = "Bound"        # Bound, Exp or Pois

# Define number of scenarios in every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h.simcores = n.sim %% n.cores
if (h.simcores != 0){
  n.simcores[1:h.simcores] = n.simcores[1:h.simcores] + 1
}

object.sim = obs.sim(n.obs, q.par, n.par, n.sim, r, sd.eps)
b          = object.sim$b
X          = object.sim$X
Y          = object.sim$Y

# All scenarios
Sys.time()
real.boot.LR = foreach(isim = 1:n.cores, 
                       .packages = c('glmnet', 'MASS')) %dopar% sim.LR.MBLR(isim, mb.type)
Sys.time()

# Close cluster
stopCluster(cl)

```

automatically created on 2018-05-28