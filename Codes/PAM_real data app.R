# -------------------------------------------------------------------------------
# Excess bond risk premia modelling with PAM
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx", "lars", 
              "scales", "tilting")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("Adapt.SCAD.ridge_withCV.r")

# Load data
tmpdata    = read.csv("BRP_data.csv",sep=",") # Bond Risk Premium data
data       = sapply(subset(tmpdata, select = c(2:(dim(tmpdata)[2]))), as.numeric)
dates      = as.Date(as.character(as.POSIXct(as.character(tmpdata[, 1]), format = "%Y%m%d")))
covariates = as.matrix(data[, 5:dim(data)[2]])
BRP2       = data[, 1]
BRP3       = data[, 2]
BRP4       = data[, 3]
BRP5       = data[, 4]

# Define pre-specified parameters
n.par       = dim(covariates)[2]           # no of parameters
n.sim       = 1                            # no of simulations
M           = 6                            # increment of observations (1y)
K           = dim(covariates)[1]/M  # number of subsamples
a           = 3.7                          # recommended value of parameter a for SCAD
sd.eps      = 1 # ZMENIT

# Simulating mutlipliers from a bounded distribution
runif.mod = function(n){
  unif1 = runif(n, min = 0, max = 12)
  unif.a = unif1[which(unif1 <= 1)]
  index.a  = which(unif1 <= 1)
  unif.b = unif1[which((1 < unif1) & (unif1 <= 4))]
  index.b  = which((1 < unif1) & (unif1 <= 4))
  unif.c = (unif1[which(unif1 > 4)] - 4)/8
  index.c  = which(unif1 > 4)
  unif2   = cbind(c(index.a, index.b, index.c), c(unif.a, unif.b, unif.c))
  unif.fin = unif2[order(unif2[, 1], decreasing = FALSE), 2]
  unif.fin
}

# For each k, simulate multipliers from different distributions
n.boot            = 10000
multipliers.exp   = list()
multipliers.bern  = list()
multipliers.bound = list()
for (k in 1:K){
  set.seed    = 20170424 + k * 10
  n.obs.multi = k * M
  multipliers.exp[[k]] = matrix(rexp(n.boot * n.obs.multi, rate = 1),
                                ncol = n.obs.multi, nrow = n.boot)
  multipliers.bern[[k]] = matrix(2 * rbern(n.boot * n.obs.multi, 0.5),
                                 ncol = n.obs.multi, nrow = n.boot)
  multipliers.bound[[k]] = matrix(runif.mod(n.boot * n.obs.multi),
                                  ncol = n.obs.multi, nrow = n.boot)
}

# Funtion to compute bootstrapped penalized log-likelihood function
loglik.pen.MB = function(beta, lambda, X, Y, multipl){
  loglik1 = -(- (1/(2 * sd.eps^2) 
                 * t((Y.tmp - X.tmp %*% beta)^2) %*% multipl)
              - ((k * M) * t(pen.prime(beta, lambda, a)) %*% abs(beta)))
  loglik1
}

# Function to compute bootstrapped unpenalized log-likelihood function
loglik.MB = function(beta, lambda, X, Y, multipl){
  loglik1 = -( - (1/(2 * sd.eps^2) 
                  * t((Y.tmp - X.tmp %*% beta)^2) %*% multipl))
  loglik1
}

# Function to compute penalized log-likelihood function
loglik.pen = function(beta, lambda, X, Y){
  loglik1 = -(- (1/(2 * sd.eps^2) * t(Y.tmp - X.tmp %*% beta) %*% (Y.tmp - X.tmp %*% beta))
              - ((k * M) * t(pen.prime(beta, lambda, a)) %*% abs(beta)))
  loglik1
}

k           = 1
k1          = 0
k2          = 1
beta        = list()
lambda      = list()
sd.res      = list()
multipliers = multipliers.bound

# 1.step: Fit the model for assumed homogeneous interval I_t^(1)
X.tmp       = covariates[((k1 * M) + 1):(k2 * M), ]
Y.tmp       = BRP5[((k1 * M) + 1):(k2 * M)]
n.obs       = length(Y.tmp)
object.tmp  = Adapt.SCAD(X.tmp, Y.tmp, a, n.obs)
beta[[k]]   = object.tmp$beta
lambda[[k]] = object.tmp$lambda
sd.res[[k]] = sqrt(var(Y.tmp - X.tmp%*%beta[[k]]))
while (k < K){
  # 2.step: Fit the model over the next interval
  k2                = k2 + 1
  X.tmp             = covariates[((k1 * M) + 1):(k2 * M), ]
  Y.tmp             = BRP5[((k1 * M) + 1):(k2 * M)]
  n.obs             = length(Y.tmp)
  object.tmp        = Adapt.SCAD(X.tmp, Y.tmp, a, n.obs)
  beta[[(k + 1)]]   = object.tmp$beta
  lambda[[(k + 1)]] = object.tmp$lambda
  sd.res[[(k + 1)]] = sqrt(var(Y.tmp - X.tmp%*%beta[[(k + 1)]]))
  
  # 3.step: Find quantiles of LR with use of MB
  lik.ratio   = numeric(0)
  sd.eps      = sd.res[[(k + 1)]]

  set.seed    = 20170424 * k * 10
  multipliers = matrix(runif.mod(n.boot * n.obs),
                       ncol = n.obs, nrow = n.boot)

  for (l in 1:(n.boot/10)){
    multipl       = multipliers[l, ]
    boot.beta     = Adapt.SCAD.MB(X.tmp, Y.tmp, a, n.obs, lambda[[(k + 1)]], multipl)$beta
    lik.boot.beta = -loglik.pen.MB(boot.beta, lambda[[(k + 1)]], X.tmp, Y.tmp, multipl)
    lik.beta      = -loglik.pen.MB(beta[[(k + 1)]], lambda[[(k + 1)]], X.tmp, Y.tmp, multipl)
    lik.ratio     = c(lik.ratio, sqrt(2 * abs(lik.boot.beta - lik.beta)))
  }
  q50.med     = quantile(lik.ratio, probs = 0.5)
  q90         = quantile(lik.ratio, probs = 0.9)
  q95         = quantile(lik.ratio, probs = 0.95)
  q97.5       = quantile(lik.ratio, probs = 0.975)

  # 4.step: Evaluate the test statistic
  sd.eps = sd.res[[k]]
  LL1 = loglik.pen(beta[[k]], lambda[[k]], X.tmp, Y.tmp)
  sd.eps = sd.res[[(k + 1)]]
  LL2 = loglik.pen(beta[[(k + 1)]], lambda[[(k + 1)]], X.tmp, Y.tmp)
  t.stat = sqrt(2 * abs(LL1 - LL2))
 
  # 5.step: Test for homogeneity and evaluate current beta estimator
  if (t.stat <= q95){
    for (k3 in (k1 + 1):k){
      beta[[k3]]   = beta[[(k + 1)]]
      lambda[[k3]] = lambda[[(k + 1)]]
      sd.res[[k3]] = sd.res[[(k + 1)]]
    }
    k                 = k + 1
  } else {
    k1          = k2 - 1
    k           = k + 1
    X.tmp       = covariates[((k1 * M) + 1):(k2 * M), ]
    Y.tmp       = BRP5[((k1 * M) + 1):(k2 * M)]
    n.obs       = length(Y.tmp)
    object.tmp  = Adapt.SCAD(X.tmp, Y.tmp, a, n.obs)
    beta[[k]]   = object.tmp$beta
    lambda[[k]] = object.tmp$lambda
    sd.res[[k]] = sqrt(var(Y.tmp - X.tmp %*% beta[[k]]))
  } 
  print(k)
}

# Plot Y.fitted vs real data
for (k in 1:K){
  Y.fitted = c(Y.fitted, covariates[(((k - 1) * M) + 1):(k * M),] %*% beta[[k]])
}

plot(Y.fitted, type = "l", col = "red", ylim = c(-0.08, 0.08))

# What is actual number of nonzero coefficients for every interval?
act.set = numeric(0)
for (k in 1:K){
  act.set = c(act.set, sum(beta[[k]] != 0)) 
}


