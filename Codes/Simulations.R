# -------------------------------------------------------------------------------
# Simulation of structural changes in a model and an adaptive estimation
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("/Users/Lenka/Documents/IRTG 1792/Penalized Adaptive Method/PAM")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("Adapt.SCAD_withCV.r")

# Simulation setup
n.obs       = 1000              # No of observations
n.par       = 10                # No  of parameters
n.sim       = 500               # No of simulations
seed1       = 20150206          # Seed simulation X
seed2       = 20150602          # Seed simulation epsilon
r           = 0.5               # Correlation parameter of X
sd.eps      = 1                 # Standard deviation of epsilon
cp1.seq     = seq(50, 950, 50)  # Define change points
max.steps   = 50                # Max no of iterations in Adapt.SCAD
lambda.grid = seq(1 * n.obs^(-1/3), 30 * n.obs^(-1/3), 1 * n.obs^(-1/3))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD
M = 50
K = 20

# Initiate cluster for parallel computing
n.cores = 4 #detectCores()          # Number of cores to be used
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
tmp1.1  = c(1, 1, 1, 1, 1)
tmp1.2  = rep(0, n.par - length(tmp1.1))
b1      = c(tmp1.1, tmp1.2)

tmp2.1  = c(2, 2, 2, 2, 2)
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
  values        = list(adapt.tmp$beta, adapt.tmp$penaltyprime, adapt.tmp$weights)
  names(values) = c("beta", "penalty", "weights")
  return(values)
}

# Function for test statistic computation
test.stat = function(t, s, k){
  X.tmp      = as.matrix(X[[s]][t:(t + k * M - 1), ])
  Y.tmp      = as.matrix(Y[[s]][t:(t + k * M - 1)])
  tmp1       = adapt.fct(t, s, k)
  tmp2       = adapt.fct(t, s, k - 1)
  beta.tmp1  = tmp1$beta
  beta.tmp2  = tmp2$beta
  pen.prime1 = tmp1$penalty
  pen.prime2 = tmp2$penalty
  
  loglik1   = (- (1/(2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp1)
               - ((k * M) * t(pen.prime1) %*% abs(beta.tmp1)))
  )
  loglik2   = (- (1/(2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
                                                      %*% (Y.tmp - X.tmp %*% beta.tmp2)
               - ((k * M) * t(pen.prime2) %*% abs(beta.tmp2)))
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

zeta.values = q95.NMB    # Choose the set of critical values to be used

# Find the change points, chpts and icorr
out.sum = list()
Sys.time()
for (num in 1:length(cp1.seq)){
  cp1 = cp1.seq[num]
  Y = Y.sim(cp1)
  final.chpts   = foreach(i = 1:n.cores, .packages = 'glmnet', .combine = 'c') %dopar% find.chpts(i, zeta.values)   
  out.sum[[num]] = output.summary(final.chpts)
  print(num)
}
Sys.time()

# Close cluster
stopCluster(cl)

# C.v.'s identified by Chen&Niu(2014) method, where the model was fitted with the MLE method
zeta.esterr.mle = c(Inf, 7.70, 4.41, 3.71, 2.97, 2.64, 2.32, 2.02, 1.96, 1.88, 1.73, 1.66, 
                    1.58, 1.58, 1.53, 1.29, 1.36, 1.35, 1.20, 1.17, 1.23, 1.18, 1.12, 1.22)

# C.v.'s identified by Chen&Niu(2014) method, where the model was fitted with the Adapt.SCAD method
zeta.esterr.scad = c(Inf, 8.69, 4.80, 2.91, 2.68, 2.12, 2.52, 1.82, 1.88, 1.62, 1.51, 1.61, 1.48,
                     1.51, 1.50, 1.28, 1.26, 1.16, 1.14, 1.40)

length(zeta95scad)
plot(zeta95scad[3:24], type = "l", ylim = c(0, max(zeta95scad[3:24])))
lines(zeta95mle, col = "blue")
lines(zeta.esterr.scad.mle, col = "red")
lines(sqrt(2) * zeta.esterr.mle, col = "green")
