# -------------------------------------------------------------------------------
# Simulation of structural changes in a model and an adaptive estimation
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
n.obs       = 1200              # No of observations
n.par       = 20                # No of parameters
n.sim       = 20                # No of simulations
seed1       = 20150206          # Seed simulation X
seed2       = 20150602          # Seed simulation epsilon
r           = 0.5               # Correlation parameter of X
sd.eps      = 1                 # Standard deviation of epsilon
cp1.seq     = seq(50, 100, 50)  # Define change points
max.steps   = 50                # Max no of iterations in Adapt.SCAD
lambda.grid = seq(1 * n.obs^(-1/3), 30 * n.obs^(-1/3), 1 * n.obs^(-1/3))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD

# Initiate cluster for parallel computing
n.cores = detectCores()         # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of scenarios in every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h          = n.sim %% n.cores
if (h != 0){
  n.simcores[1:h] = n.simcores[1:h] + 1
}

# True beta coefficients with change in t = cp1
tmp1.1  = c(1, 1.5, 1, 1, 2, -3, -1.5, 1, 0, 5, 0, 1)
tmp1.2  = rep(0, n.par - length(tmp1.1))
b1      = c(tmp1.1, tmp1.2)

tmp2.1  = c(1, 1.5, 0, 0, 2, -3, -1.5, 0, -4, 5)
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

# Adaptive method for multistep SCAD estimators
adapt.fct = function(t, s, k){
  X.tmp         = X[[s]][t:(t + k * M - 1), ]
  Y.tmp         = Y[[s]][t:(t + k * M - 1)]
  adapt.tmp     = Adapt.SCAD(X.tmp, Y.tmp, a, k * M, max.steps) # Adaptive SCAD fit
  values        = list(adapt.tmp$beta, adapt.tmp$penalty, adapt.tmp$weights)
  names(values) = c("beta", "penalty", "weights")
  return(values)
}

# Function for test statistic computation
test.stat = function(t, s, k){
  X.tmp     = as.matrix(X[[s]][t:(t + k * M - 1), ])
  Y.tmp     = as.matrix(Y[[s]][t:(t + k * M - 1)])
  tmp1      = adapt.fct(t, s, k)
  tmp2      = adapt.fct(t, s, k - 1)
  beta.tmp1 = tmp1$beta
  beta.tmp2 = tmp2$beta
  
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp1)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp1))
               # - ((k * M) * sum(tmp1$penalty))
               )
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                                      * t(Y.tmp - X.tmp %*% beta.tmp2)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp2))
               # - ((k * M) * sum(tmp2$penalty))
               )
  test.stat = sqrt(2 * abs(loglik1 - loglik2))
  test.stat
}

# Function for finding change points
find.chpts = function(h, zeta){
  s1 = ifelse(h == 1, 1, sum(n.simcores[1:(h - 1)],1))
  s2 = sum(n.simcores[1:h])
  b.pts = list()
  for (s in s1:s2){
    t = 1
    k = 2
    m = 2
    l = 1
    theta = list()
    theta[[1]] = adapt.fct(t, s, 1)$beta
    b.pts.tmp = 0
    while (m <= K && t < (K * M)) {
      if (test.stat(t, s, k) <= zeta[k]) {
        l = m
        theta[[m]] = adapt.fct(t, s, k)$beta
        k = k + 1
        m = m + 1
      } else if ((l * M + 1) < n.obs) {
        t = l * M + 1
        l = m
        b.pts.tmp = c(b.pts.tmp, t)
        k = 2
        theta[[m]] = adapt.fct(t, s, 1)$beta
        m = m + 1
      } else {
        m = m + 1
      }
    }
    stmp = s - s1 + 1
    b.pts[[stmp]] = b.pts.tmp
  }
  return(b.pts)
}

# Simulation output summary
output.summary = function(final.chpts){
  temp  = seq(51, 1151, 50)
  
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
  # (i.e. whether the change point was identifiedand whether it was 
  #  the first one to be found, which would correspond the simulations setup)
  icorr = 0
  for (i in 1:n.sim){
    if (!is.na(final.chpts[[i]][2]) && final.chpts[[i]][2] == (cp1 + 1))  {
      icorr = icorr + 1
    }
  }
  
  values = list(chpts, icorr)
  names(values) = c("chpts", "icorr")
  return(values)
}

zeta.values = zeta95    # Choose the set of critical values to be used

# Find the change points, chpts and icorr
out.sum = list()
for (num in 1:length(cp1.seq)){
  cp1 = cp1.seq[num]
  Y = Y.sim(cp1)
  final.chpts   = foreach(i = 1:n.cores, .packages = 'glmnet', .combine = 'c') %dopar% find.chpts(i, zeta.values)   
  out.sum[[num]] = output.summary(final.chpts)
  print(num)
}


# Close cluster
stopCluster(cl)


# In order to save time, here are some of the identified critical values

# C.v.'s identified by multiplier bootstrap method
zeta95 = c(1770.02275, 49.36315, 14.64020, 14.17828, 13.87526, 13.70207, 13.67552, 13.71178, 13.88155, 13.80665, 
           13.75727, 13.74442, 13.68032, 13.62109, 13.56758, 13.51338, 13.46749, 13.41742, 13.37000, 13.33120,
           13.29143, 13.25637, 13.22196, 13.18458)

# C.v.'s identified by Chen&Niu(2014) method, where the model was fitted with the MLE method
zeta.esterr.mle = c(Inf, 7.70, 4.41, 3.71, 2.97, 2.64, 2.32, 2.02, 1.96, 1.88, 1.73, 1.66, 
                    1.58, 1.58, 1.53, 1.29, 1.36, 1.35, 1.20, 1.17, 1.23, 1.18, 1.12, 1.22)

# C.v.'s identified by Chen&Niu(2014) method, where the model was fitted with the Adapt.SCAD method
zeta.esterr.scad.mle = c(Inf, 6.65, 4.09, 3.32, 2.68, 2.33, 1.93, 1.81, 1.69, 1.88, 1.54, 
                         1.56, 1.56, 1.34, 1.36, 1.20, 1.31, 1.30, 1.17, 1.07, 1.07, 1.05, 1.12, 0.92)

