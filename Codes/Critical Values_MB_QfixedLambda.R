# -------------------------------------------------------------------------------
# Critical values calibration using Multiplier Bootstrap and SCAD
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("/Users/Lenka/Documents/IRTG 1792/Penalized Adaptive Method/PAM")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("Adapt.SCAD_withCV.r")

# Simulation setup
n.obs       = 1000            # no of observations
n.par       = 10              # no of parameters
n.sim       = 1            # no of simulations
seed1       = 20150206        # seed simulation X
seed2       = 20150602        # seed simulation epsilon
M           = 50              # increment of observations between successive subsamples
K           = 20              # number of subsamples
sd.eps      = 1               # st. dev. of the error term
max.steps   = 50                # Max no of iterations in Adapt.SCAD
lambda.grid = seq(1 * n.obs^(-1/3), 30 * n.obs^(-1/3), 1 * n.obs^(-1/3))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD


# True beta coefficients (homogeneous for all t = 1, ..., 1000)
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

# Computation of Y for t = 1, ..., n.obs
Y    = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in 1:n.obs){
    Y.tmp = c(Y.tmp, b %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y[[i]] = Y.tmp
}

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
  betas   = list()
  wghts   = list()
  pens    = list()
  lambdas = list()
  s1 = ifelse(l == 1, 1, sum(n.simcores[1:(l - 1)],1))
  s2 = sum(n.simcores[1:l])
  
  for (s in s1:s2){
    ad.beta    = numeric(0)
    ad.penalty = numeric(0)
    ad.weights = numeric(0)
    ad.lambda  = numeric(0)
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
    lambdas[[stmp]]  = ad.lambda
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
betas.tmp   = foreach(i = 1, .packages = "glmnet") %dopar% find.betas(i)   # Parallel computing
Sys.time()
final.betas = res.betas(1, betas.tmp) 


betas   = final.betas$betas
pens    = final.betas$pens
wghts   = final.betas$wghts
lambdas = final.betas$lambdas

# Take 1 simulated scenario and create lists of betas, X's and Y's
s = 1

# For each k, simulate multipliers from N(1,1)
n.boot      = 10000
multipliers.norm = list()
multipliers.exp  = list()
multipliers.bern = list()
for (k in 1:K){
  set.seed    = 20170424 + k * 10
  n.obs.multi = k * M
  multipliers.norm[[k]] = matrix(rnorm(n.boot * n.obs.multi, mean = 1, sd = 1), 
                                 ncol = n.obs.multi, nrow = n.boot)
  multipliers.exp[[k]] = matrix(rexp(n.boot * n.obs.multi, rate = 1), 
                                ncol = n.obs.multi, nrow = n.boot)
  multipliers.bern[[k]] = matrix(2 * rbern(n.boot * n.obs.multi, 0.5), 
                                 ncol = n.obs.multi, nrow = n.boot)
}

# Compute penalized log-likelihood function
loglik.pen = function(beta){
  X = X[[s]][1:(k * M), ]
  Y = Y[[s]][1:(k * M)]
  lambda  = lambdas[[s]][k]
  loglik1 = -(- (1/(2 * sd.eps^2) 
              * t((Y - X %*% beta)^2) %*% multipliers[[k]][l, ])
              - ((k * M) * t(pen.prime(beta, lambda, a)) %*% abs(beta)))
  loglik1
}

# Compute log-likelihood function
loglik = function(beta){
  X = X[[s]][1:(k * M), ]
  Y = Y[[s]][1:(k * M)]
  loglik1 = -( - (1/(2 * sd.eps^2) 
               * t((Y - X %*% beta)^2) %*% multipliers[[k]][l, ]))
  loglik1
}

# Compute log-likelihood function
loglik.pen.estim = function(beta){
  X.tmp   = X[[s]][1:(k * M), ]
  Y.tmp   = Y[[s]][1:(k * M)]
  lambda  = lambdas[[s]][k]
  loglik1 = -(- (1/(2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta) %*% (Y.tmp - X.tmp %*% beta))
             - ((k * M) * t(pen.prime(beta, lambda, a)) %*% abs(beta)))
  loglik1
}


# Compute likelihood ratio
q90           = numeric(0)
q95           = numeric(0)
q97.5         = numeric(0)
q50.med       = numeric(0)
lik.ratio.tmp = numeric(0)
lik.ratio     = list()
multipliers   = multipliers.exp
Sys.time()
for (k in 1:K){
  for (l in 1:10000){
    X.tmp = X[[s]][1:(k * M), ]
    Y.tmp = Y[[s]][1:(k * M)]
    multipl       = multipliers[[k]][l, ]
    boot.beta     = Adapt.SCAD.MB(X.tmp, Y.tmp, a, (k * M), max.steps, lambdas[[s]][k], multipl)$beta
    lik.boot.beta = -loglik.pen(boot.beta)
    lik.beta      = -loglik.pen(betas[[s]][, k])
    # lik.boot.beta = -loglik.pen.estim(betas[[s]][, k])
    # lik.beta      = -loglik.pen.estim(b)
    lik.ratio.tmp = c(lik.ratio.tmp, sqrt(2 * abs(lik.boot.beta - lik.beta)))
  }
  q50.med[k]     = quantile(lik.ratio.tmp, probs = 0.5)
  q90[k]         = quantile(lik.ratio.tmp, probs = 0.9)
  q95[k]         = quantile(lik.ratio.tmp, probs = 0.95)
  q97.5[k]       = quantile(lik.ratio.tmp, probs = 0.975)
  lik.ratio[[k]] = lik.ratio.tmp
  print(k)
}
Sys.time()

# LR without using MB
q50.med.WMB.SCAD = c(2.468439, 2.351712, 2.276810, 2.233094, 2.207786, 2.183878,
                     2.168143, 2.155139, 2.146212, 2.137692, 2.131062, 2.122851,
                     2.115663, 2.109319, 2.106777, 2.104328, 2.102323, 2.099543,
                     2.098331, 2.096179, 2.094335, 2.093272, 2.092677, 2.092000)

q90.WMB.SCAD = c(3.614706, 3.479302, 3.372542, 3.302880, 3.270328, 3.240878, 3.221200,
                 3.203104, 3.194546, 3.183214, 3.172650, 3.167305, 3.155796, 3.149053,
                 3.141527, 3.137294, 3.133713, 3.128158, 3.124915, 3.120908, 3.114370,
                 3.110735, 3.107748, 3.105128)

q95.WMB.SCAD = c(3.923716, 3.786261, 3.686837, 3.626303, 3.589884, 3.568400, 3.542780,
                 3.519332, 3.508908, 3.500250, 3.497979, 3.488158, 3.478392, 3.473312,
                 3.467957, 3.460819, 3.454821, 3.443767, 3.438317, 3.432866, 3.423380,
                 3.418768, 3.414895, 3.412546)

q97.5.WMB.SCAD = c(4.181730, 4.095195, 3.986217, 3.920089, 3.876683, 3.833833, 3.805137,
                   3.782415, 3.766147, 3.760804, 3.748642, 3.736573, 3.732357, 3.723796,
                   3.717515, 3.713100, 3.706501, 3.701588, 3.698627, 3.694716, 3.691255,
                   3.690148, 3.687200, 3.684327)

# Exponential bootstrap
q50.med.EMB.SCAD = c(2.164130, 2.296380, 2.327368, 2.275602, 2.335753, 2.390573, 2.437987,
                     2.446342, 2.454642, 2.461037, 2.466647, 2.464767, 2.461160, 2.452747,
                     2.422218, 2.406367, 2.391837, 2.379031, 2.366472, 2.354068)

q90.EMB.SCAD = c(3.348625, 3.544231, 3.503432, 3.424936, 3.498878, 3.561059, 3.617529,
                 3.613662, 3.609958, 3.608140, 3.610149, 3.597470, 3.584900, 3.568074,
                 3.540151, 3.521611, 3.502642, 3.487727, 3.472270, 3.456735)

q95.EMB.SCAD = c(3.712119, 3.954081, 3.905580, 3.817001, 3.904789, 3.962769, 4.018149,
                 4.009528, 3.999589, 3.994551, 3.992841, 3.982880, 3.968954, 3.948051,
                 3.916918, 3.895561, 3.876417, 3.858285, 3.841037, 3.823421)

q97.5.EMB.SCAD = c(4.049922, 4.315286, 4.278419, 4.181767, 4.269420, 4.320798, 4.392081,
                   4.371225, 4.359111, 4.355521, 4.349058, 4.334140, 4.316168, 4.295016,
                   4.269420, 4.245842, 4.222856, 4.201594, 4.181855, 4.163401)

# 2*Bernoulli bootstrap p = 10, q = 5
q50.med.BMB.SCAD = c(2.816530, 2.805769, 2.711674, 2.587187, 2.617995, 2.638324, 2.665609,
                     2.655011, 2.644045, 2.636327, 2.632091, 2.619231, 2.607464, 2.591056,
                     2.553871, 2.531563, 2.510215, 2.489497, 2.471838, 2.456186)

q90.BMB.SCAD = c(4.005971, 4.043978, 3.937047, 3.815622, 3.825688, 3.826862, 3.842862,
                 3.816240, 3.790899, 3.774387, 3.764889, 3.748545, 3.731108, 3.709807,
                 3.678856, 3.655294, 3.634070, 3.613339, 3.594382, 3.576308)

q95.BMB.SCAD = c(4.352748, 4.406802, 4.314485, 4.199027, 4.199027, 4.198332, 4.219030,
                 4.193563, 4.163640, 4.141536, 4.129233, 4.109343, 4.088396, 4.065234,
                 4.035715, 4.011674, 3.990454, 3.968733, 3.949009, 3.929815)

q97.5.BMB.SCAD = c(4.651579, 4.711605, 4.639481, 4.533354, 4.530692, 4.534132, 4.565654,
                   4.537856, 4.505367, 4.480148, 4.465392, 4.443598, 4.422051, 4.399504,
                   4.369513, 4.343949, 4.316690, 4.293845, 4.273463, 4.256487)

