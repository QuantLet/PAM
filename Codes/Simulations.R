# -------------------------------------------------------------------------------
# Simulation of structural changes in a model and an adaptive estimation
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
n.par      = 20              # no of parameters
n.sim      = 100             # no of simulations
seed1      = 20150206        # seed simulation X
seed2      = 20150602        # seed simulation epsilon
r          = 0.5             # Correlation parameter of X
sd.eps     = 1               # Standard deviation of epsilon
cp1        = 300
cp2        = 600
cp3        = 900

# True beta coefficients with change in t = {cp1, cp2, cp3}
tmp1.1  = c(1, 1.5, 1, 1, 2, -3, -1.5, 1, 0, 5, 0, 1)
tmp1.2  = rep(0, n.par - length(tmp1.1))
b1      = c(tmp1.1, tmp1.2)

tmp2.1  = c(1, 1.5, 0, 0, 2, -3, -1.5, 0, -4, 5)
tmp2.2  = rep(0, n.par - length(tmp2.1))
b2      = c(tmp2.1, tmp2.2)

tmp3.1  = c(2, 1)
tmp3.2  = rep(0, n.par - length(tmp3.1))
b3      = c(tmp3.1, tmp3.2)

tmp4.1  = c(1, 1.5, 0, 0, 0, 3, 0, 1.5)
tmp4.2  = rep(0, n.par - length(tmp4.1))
b4      = c(tmp4.1, tmp4.2)

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
  for (j in (cp1 + 1):cp2){
    Y.tmp = c(Y.tmp, b2 %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y2[[i]] = Y.tmp
}

# Computation of Y for t = cp2, ..., cp3
Y3    = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in (cp2 + 1):cp3){
    Y.tmp = c(Y.tmp, b3 %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y3[[i]] = Y.tmp
}

# Computation of Y after t = cp4
Y4    = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in (cp3 + 1):n.obs){
    Y.tmp = c(Y.tmp, b4 %*% X[[i]][j, ] + eps[[i]][j])
  }
  Y4[[i]] = Y.tmp
}

Y   = list()
for (i in 1:n.sim){
  Y[[i]] = c(Y1[[i]], Y2[[i]], Y3[[i]], Y4[[i]])
}

# ----------------------------------------------------------------------------------------
# Adaptive method for OLS estimators
# ----------------------------------------------------------------------------------------
# Set of critical values zeta, fitted betas and test statistics
beta.fct = function(t, s, k){
  X.tmp = X[[s]][t:(t + (k * M) - 1), ]
  Y.tmp = Y[[s]][t:(t + (k * M) - 1)]
  beta.fit  = solve(t(X.tmp) %*% X.tmp) %*% t(X.tmp) %*% Y.tmp
  beta.fit
}

test.stat = function(t, s, k){
  X.tmp   = as.matrix(X[[s]][t:(t + (k * M) - 1), ])
  Y.tmp   = as.matrix(Y[[s]][t:(t + (k * M) - 1)])
  beta.tmp1 = as.matrix(beta.fct(t, s, k))
  beta.tmp2 = as.matrix(beta.fct(t, s, k - 1))
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp1)
                  %*% (Y.tmp - X.tmp %*% beta.tmp1)))
  
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) 
               - (1 / (2 * sig^2) * t(Y.tmp - X.tmp %*% beta.tmp2)
                  %*% (Y.tmp - X.tmp %*% beta.tmp2)))
  
  test.stat = sqrt(abs(loglik1 - loglik2))
  test.stat
}

b.pts = list()
for (s in 1:100){
  t = 1
  k = 2
  m = 2
  theta = list()
  theta[[1]] = beta.fct(t, s, 1)
  b.pts.tmp = 0
  while (m <= K && t < 1200) {
    if (test.stat(t, s, k) <= zeta.riskbound.mle[k]) {
      l = m
      theta[[m]] = beta.fct(t, s, k)
      k = k + 1
      m = m + 1
    } else if ((l * M + 1) < (K * M)) {
      t = l * M + 1
      l = m
      b.pts.tmp = c(b.pts.tmp, t)
      k = 2
      theta[[m]] = beta.fct(t, s, 1)
      m = m + 1
    } else {
      m = m + 1
    }
  }
  b.pts[[s]] = b.pts.tmp
}

c2 = 0
for (i in 1:100){
  if (length(b.pts[[i]])==4) {
    print(i)
    c2 = c2 + 1
  }
}
# zeta.riskbound.mle funguje lepsie ako zeta.esterr.mle

plot(theta[[1]], type = "l", col = "black", ylim = c(-5, 5))
lines(theta[[2]], col = "black")
lines(theta[[3]], col = "black")
lines(theta[[4]], col = "black")
lines(theta[[5]], col = "black")
lines(theta[[6]], col = "black")
lines(theta[[7]], col = "blue")
lines(theta[[8]], col = "blue")
lines(theta[[9]], col = "blue")
lines(theta[[10]], col = "blue")
lines(theta[[11]], col = "blue")
lines(theta[[12]], col = "blue")
lines(theta[[13]], col = "red")
lines(theta[[14]], col = "red")
lines(theta[[15]], col = "red")
lines(theta[[16]], col = "red")
lines(theta[[17]], col = "red")
lines(theta[[18]], col = "red")
lines(theta[[19]], col = "green")
lines(theta[[20]], col = "green")
lines(theta[[21]], col = "green")
lines(theta[[22]], col = "green")
lines(theta[[23]], col = "green")
lines(theta[[24]], col = "green")
# ----------------------------------------------------------------------------------------
# Adaptive method for multistep SCAD estimators
# ----------------------------------------------------------------------------------------
# Set of critical values zeta, fitted betas and test statistics
adapt.fct = function(t, s, k){
  X.tmp = X[[s]][t:(t + k * M - 1), ]
  Y.tmp = Y[[s]][t:(t + k * M - 1)]
  adapt.tmp = Adapt.SCAD(X.tmp, Y.tmp, a, k * M)
  values        = list(adapt.tmp$beta, adapt.tmp$penalty, adapt.tmp$weights)
  names(values) = c("beta", "penalty", "weights")
  return(values)
}

test.stat = function(t, s, k){
  X.tmp   = as.matrix(X[[s]][t:(t + k * M - 1), ])
  Y.tmp   = as.matrix(Y[[s]][t:(t + k * M - 1)])
  tmp1    = adapt.fct(t, s, k)
  tmp2    = adapt.fct(t, s, k - 1)
  beta.tmp1 = tmp1$beta
  beta.tmp2 = tmp2$beta
  
  loglik1   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                             * t(Y.tmp - X.tmp %*% beta.tmp1)
                                             %*% (Y.tmp - X.tmp %*% beta.tmp1))
               - ((k * M) * t(tmp1$penalty)
                  %*% abs(beta.tmp1 * tmp1$weights)))
  loglik2   = (-(k * M)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
                                               * t(Y.tmp - X.tmp %*% beta.tmp2)
                                               %*% (Y.tmp - X.tmp %*% beta.tmp2))
               - ((k * M) * t(tmp2$penalty)
                  %*% abs(beta.tmp2 * tmp2$weights)))
  test.stat = sqrt(abs(loglik1 - loglik2))
  test.stat
}

b.pts = list()
for (s in 1:n.sim){
  t = 1
  k = 2
  m = 2
  theta = list()
  theta[[1]] = adapt.fct(t, s, 1)$beta
  b.pts.tmp = 0
  while (m <= K && t < (K * M)) {
    if (test.stat(t, s, k) <= zeta[k]) {
      l = m
      theta[[m]] = adapt.fct(t, s, k)$beta
      k = k + 1
      m = m + 1
    } else if ((l * M + 1) < 1200) {
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
  b.pts[[s]] = b.pts.tmp
}

b.pts.pre = b.pts

i1 = numeric(0)
i2 = numeric(0)
for (i in 1:1000){
  if (all(c(301, 601, 901) %in% b.pts[[i]])) {
    i2 = c(i2, i)
  }
  if (all(c(301, 601, 901) %in% b.pts[[i]]) && length(b.pts[[i]]) == 4){
    i1 = c(i1, i)
  }
}

length(i1)
length(i2)
c2
