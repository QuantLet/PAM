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
n.obs       = 500              # No of observations
n.par       = 10                # No  of parameters
n.sim       = 1              # No of simulations
seed1       = 20171110          # Seed simulation X
seed2       = 20170602          # Seed simulation epsilon
r           = 0.5               # Correlation parameter of X
sd.eps      = 1                 # Standard deviation of epsilon
cp1.seq     = c(50, 100, 200, 400)  # Define change points
a           = 3.7             # Recommended value of parameter a for SCAD
n.boot      = 1000
M           = 50
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

# True beta coefficients with change in t = cp1
tmp1.1  = rep(1, 5)
tmp1.2  = rep(0, n.par - length(tmp1.1))
b1      = c(tmp1.1, tmp1.2)

tmp2.1  = seq(1, 0.2, -0.2)
tmp2.2  = rep(0, n.par - length(tmp2.1))
b2      = c(tmp2.1, tmp2.2)

# Simulation of the design matrix
mu    = rep(0, n.par)
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
for (i in 1:n.sim){
  set.seed(seed1)
  X[[i]] = mvrnorm(n = n.obs, mu, Sigma)
}  

# Simulation of the error term for t = 1, ..., n.obs
eps  = list()
set.seed(seed2)
for (i in 1:n.sim){
  eps[[i]] = rnorm(n.obs, mean = 0, sd = sd.eps)
} 

# Function for computing Y for t = 1, ..., n.obs
Y.sim = function(cp1){
  # Computation of Y for t = 1, ..., cp1
  Y1    = list()
  for (i in 1:n.sim){
    Y.tmp = numeric(0)
    for (j in 1:cp1){
      Y.tmp = c(Y.tmp, b1 %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y1[[i]] = Y.tmp
  }
  
  # Computation of Y for t = cp1, ..., cp2
  Y2    = list()
  for (i in 1:n.sim){
    Y.tmp = numeric(0)
    for (j in (cp1 + 1):n.obs){
      Y.tmp = c(Y.tmp, b2 %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y2[[i]] = Y.tmp
  }
  
  Y   = list()
  for (i in 1:n.sim){
    Y[[i]] = c(Y1[[i]], Y2[[i]])
  }
  Y
}

# Signal-to-noise ratio
(t(b1) %*% Sigma %*% b1) / (sd.eps^2)
(t(b2) %*% Sigma %*% b2) / (sd.eps^2)

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

sim.change.points = function(h){
  s1 = ifelse(h == 1, 1, sum(n.simcores[1:(h - 1)],1))
  s2 = sum(n.simcores[1:h])
  b.pts = list()
  for (s in s1:s2){
    k           = 1
    k1          = 0
    k2          = 1
    beta        = list()
    beta0.tmp   = list()
    lambda      = list()
    b.pts.tmp   = 0
 
    # 1.step: Fit the model for assumed homogeneous interval I_t^(1)
    X.tmp.L        = X[[s]][((k1 * M) + 1):(k2 * M), ]
    Y.tmp.L        = Y[[s]][((k1 * M) + 1):(k2 * M)]
    n.obs.tmp.L    = length(Y.tmp.L)
    object.tmp.L   = Adapt.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
    beta[[k]]      = object.tmp.L$beta
    lambda[[k]]    = object.tmp.L$lambda
    beta0.tmp[[k]] = object.tmp.L$beta.0
    
    Sys.time()
    while (k < K){
      # 1.step: Fit the model for assumed homogeneous interval I_t^(k)
      X.tmp.L       = X[[s]][((k1 * M) + 1):(k2 * M), ]
      Y.tmp.L       = Y[[s]][((k1 * M) + 1):(k2 * M)]
      n.obs.tmp.L   = length(Y.tmp.L)
      object.tmp.L  = Adapt.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
      beta.tmp.L    = object.tmp.L$beta
      lambda.tmp.L  = object.tmp.L$lambda
      beta0.tmp.L   = object.tmp.L$beta.0
     
      # 2. a) step: Fit the model over the next interval I_t^(k + 1) - I_t^(k)
      X.tmp.R      = X[[s]][((k2 * M) + 1):((k2 + 1) * M), ]
      Y.tmp.R      = Y[[s]][((k2 * M) + 1):((k2 + 1) * M)]
      n.obs.tmp.R  = length(Y.tmp.R)
      object.tmp.R = Adapt.SCAD(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R)
      beta.tmp.R   = object.tmp.R$beta
      lambda.tmp.R = object.tmp.R$lambda
      beta0.tmp.R  = object.tmp.R$beta.0
      
      # 2. b) step: Fit the model over the whole interval I_t^(k + 1)
      X.tmp.T      = rbind(X.tmp.L, X.tmp.R)
      Y.tmp.T      = c(Y.tmp.L, Y.tmp.R)
      n.obs.tmp.T  = n.obs.tmp.L + n.obs.tmp.R
      object.tmp.T = Adapt.SCAD(X.tmp.T, Y.tmp.T, a, n.obs.tmp.T)
      beta.tmp.T   = object.tmp.T$beta
      lambda.tmp.T = object.tmp.T$lambda
      beta0.tmp.T  = object.tmp.T$beta.0
      
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
      set.seed      = 20170424 * k * s * 10
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
      }
      q50.med     = quantile(lik.ratio.fix.pen, probs = 0.5)
      q90         = quantile(lik.ratio.fix.pen, probs = 0.9)
      q95         = quantile(lik.ratio.fix.pen, probs = 0.95)
      q97.5       = quantile(lik.ratio.fix.pen, probs = 0.975)
         
      k2 = k2 + 1
      k  = k + 1
      
      # 5.step: Test for homogeneity and evaluate current beta estimator
      if (t.stat.pen <= q95){
        for (k3 in (k1 + 1):k){
          beta[[k3]]   = beta.tmp.T
          lambda[[k3]] = lambda.tmp.T
        }
      } else {
        b.pts.tmp   = c(b.pts.tmp, (((k - 1) * M) + 1))
        k1          = k2 - 1
        beta[[k]]   = beta.tmp.R
        lambda[[k]] = lambda.tmp.R
      } 
      print(k)
    }
    Sys.time()
    sim.no          = s - s1 + 1
    b.pts[[sim.no]] = b.pts.tmp 
  }
  return(b.pts)
}

# Simulation output summary
output.summary = function(final.chpts, cp1){
  temp  = seq((M + 1), (n.obs - M + 1), M)
  
  # chpts defines number of change points identified for every interval I_1, ..., I_M
  chpts = rep(0, length(temp))
  for (j in 1:length(temp)){
    for (s in 1:n.sim){
      if (temp[j] %in% final.chpts[[s]]){
        chpts[j] = chpts[j] + 1
      }
    }
  }
  
  # icorr defines number of correctly identified change points 
  # (i.e. whether the change point was identified and whether it was 
  #  the first one to be found, which would correspond the simulations setup)
  icorr = 0
  for (i in 1:n.sim){
    if (final.chpts[[i]][2] == (cp1 + 1) && !is.na(final.chpts[[i]][2]))  {
      icorr = icorr + 1
    }
  }
  
  values = list(chpts, icorr)
  names(values) = c("chpts", "icorr")
  return(values)
}


# Find the change points, chpts and icorr
out.sum    = list()
final.list = list()
Sys.time()
# n.sim = 2
# num = 1
# K = 3
for (num in 1:length(cp1.seq)){
  cp1 = cp1.seq[num]
  Y = Y.sim(cp1)    
  final.chpts       = foreach(ipar = 1:n.cores, .packages = c('glmnet', 'optimx'), .combine = 'c') %dopar% sim.change.points(ipar)   
  final.list[[num]] = final.chpts
  out.sum[[num]]    = output.summary(final.chpts, cp1)
  print(num)
}
Sys.time()

# Close cluster
stopCluster(cl)

