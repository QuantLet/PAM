# -------------------------------------------------------------------------------
# Critical values calibration using Multiplier Bootstrap and MLE
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

# setwd("")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "optimx", "LaplacesDemon")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Simulation setup
n.obs      = 1000            # no of observations
n.par      = 20              # no of parameters
n.sim      = 1000               # no of simulations
seed1      = 20150206        # seed simulation X
seed2      = 20150602        # seed simulation epsilon
M          = 50              # increment of observations between successive subsamples
K          = 20              # number of subsamples
sd.eps     = 1               # st. dev. of the error term

# True beta coefficients (homogeneous for all t = 1, ..., n.obs)
tmp1  = rep(1, 20)
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
    beta.fit  = solve(t(X.tmp) %*% X.tmp) %*% t(X.tmp) %*% Y.tmp # vymenit za MLE
    # beta.fit  = as.vector(coef(optimx(rep(0.1, 5), loglikestim, method=c("BFGS"))))
    betas.tmp = cbind(betas.tmp, beta.fit)
  }
  betas[[s]] = betas.tmp 
}

# Take 1 simulated scenario and create lists of betas, X's and Y's
# s = 1

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

# Define bootstrapped log likelihood function
loglik = function(beta){
  X = X[[s]][1:(k * M), ]
  Y = Y[[s]][1:(k * M)]
  loglik1 = -( - (1/(2 * sd.eps^2)
               * t((Y - X %*% beta)^2) %*% multipliers[[k]][l, ]))
  as.numeric(loglik1)
}

# Define log likelihood function
loglik1 = function(beta){
  X = X[[s]][1:(k * M), ]
  Y = Y[[s]][1:(k * M)]
  loglik1 = -( - (1/(2 * sd.eps^2)
                  * t(Y - X %*% beta) %*% (Y - X %*% beta)))
  as.numeric(loglik1)
}

q90   = numeric(0)
q95   = numeric(0)
q97.5 = numeric(0)
q50.med   = numeric(0)
lik.ratio.tmp = numeric(0)
lik.ratio  = list()
multipliers = multipliers.norm
Sys.time()
for (k in 1:K){
  for (s in 1:n.sim){
    # boot.beta     = as.vector(coef(optimx(rep(1, n.par), loglik, method=c("BFGS"))))
    # lik.boot.beta = -loglik(boot.beta)
    #lik.beta      = -loglik(betas[[s]][, k])
    lik.boot.beta   = -loglik1(betas[[s]][, k])
    lik.beta        = -loglik1(b)
    lik.ratio.tmp   = c(lik.ratio.tmp, sqrt(2 * abs(lik.boot.beta - lik.beta)))
  }
  q50.med[k]     = quantile(lik.ratio.tmp, probs = 0.5)
  q90[k]         = quantile(lik.ratio.tmp, probs = 0.9)
  q95[k]         = quantile(lik.ratio.tmp, probs = 0.95)
  q97.5[k]       = quantile(lik.ratio.tmp, probs = 0.975)
  lik.ratio[[k]] = lik.ratio.tmp
  print(k)
}
Sys.time()

# Close cluster
stopCluster(cl)

zeta90   = q90    # C.v.'s defined by 90 % quantile
zeta95   = q95    # C.v.'s defined by 95 % quantile
zeta97.5 = q97.5  # C.v.'s defined by 97.5 % quantile

zeta90.mbmle = c(62.335280, 21.831584, 9.826057, 6.309536, 5.326342, 5.072161,
                 4.947249, 4.874284, 4.821878, 4.783719, 4.751660, 4.728847,
                 4.707288, 4.691725, 4.679660, 4.666294, 4.655433, 4.644915,
                 4.460106, 4.447540, 4.433778, 4.429444, 4.426873, 4.421670)

zeta95.mbmle = c(124.426184, 64.862442, 36.547382, 21.882525, 14.004713, 9.825995,
                 7.515807, 6.349502, 5.872919, 5.628082, 5.489680, 5.386063,
                 5.315502, 5.262278, 5.222190, 5.189325, 5.158804, 5.133013,
                 4.793418, 4.778509, 4.751524, 4.747599, 4.742204, 4.735851)

zeta975.mbmle = c(205.921336, 126.683578, 88.040199, 64.861718, 47.285695, 36.545890,
                  28.203231, 21.882335, 17.600535, 14.004404, 11.473446, 9.825964,
                  8.503063, 7.515776, 6.864645, 6.426817, 6.160208, 5.993371, 5.085935,
                  5.078349, 5.050726, 5.037851, 5.028184, 5.013502)

# Bez bootstrapove vysledky
q50.med.WMB = c(2.134521, 2.123462, 2.142640, 2.139456, 2.136819, 2.124787, 2.118479,
                2.113203, 2.110592, 2.107122, 2.104192, 2.103784, 2.104192, 2.105337,
                2.106733, 2.107884, 2.108641, 2.108947, 2.108806, 2.107515)

q90.WMB = c(3.034387, 3.053162, 3.056071, 3.042638, 3.060374, 3.061246, 3.051861,
            3.048830, 3.043923, 3.044635, 3.039544, 3.036385, 3.034257, 3.032448,
            3.032615, 3.032861, 3.033918, 3.034587, 3.036371, 3.036804)

q95.WMB = c(3.359526, 3.358533, 3.350934, 3.334379, 3.352901, 3.353074, 3.348185,
            3.338373, 3.338373, 3.338373, 3.340363, 3.337988, 3.333776, 3.334280,
            3.338080, 3.337988, 3.335117, 3.334970, 3.334970, 3.333726)

q97.5.WMB = c(3.618168, 3.597040, 3.597715, 3.588398, 3.588398, 3.583570, 3.582566,
              3.570430, 3.582566, 3.573030, 3.572293, 3.564956, 3.563634, 3.564956,
              3.562711, 3.562694, 3.562932, 3.560836, 3.560126, 3.561746)

# Normal bootstrap
q50.med.NMB = c(1.915715, 1.902974, 1.932203, 1.977069, 2.043084, 2.079298, 2.108219,
                2.124590, 2.136858, 2.143445, 2.153035, 2.164857, 2.170366, 2.172813,
                2.177307, 2.178891, 2.180456, 2.180278, 2.179779, 2.177185)

q90.NMB = c(3.041329, 2.944994, 2.940347, 2.971782, 3.074298, 3.154287, 3.192682,
            3.215052, 3.215413, 3.221387, 3.223632, 3.228461, 3.228634, 3.228876,
            3.231607, 3.242333, 3.233683, 3.233401, 3.228730, 3.224933)

q95.NMB = c(3.590353, 3.373618, 3.253012, 3.272038, 3.408464, 3.528514, 3.568569,
            3.581989, 3.580367, 3.585199, 3.575408, 3.579830, 3.579830, 3.580329,
            3.578033, 3.579830, 3.570589, 3.565451, 3.563475, 3.556526)

q97.5.NMB = c(4.074683, 3.790612, 3.692497, 3.626351, 3.747923, 3.866353, 3.888976,
              3.911209, 3.893995, 3.901638, 3.893749, 3.888385, 3.886641, 3.888284,
              3.888967, 3.888676, 3.875770, 3.875055, 3.871446, 3.864804)

# 2*Bernoulli bootstrap
q50.med.BMB = c(1.985151, 1.965416, 2.003293, 2.021584, 2.079625, 2.111503, 2.126868,
                2.145503, 2.160076, 2.168303, 2.173983, 2.180875, 2.182387, 2.187252,
                2.189630, 2.187216, 2.189496, 2.190008, 2.185872, 2.186870)

q90.BMB = c(2.886192, 2.875557, 2.918724, 2.941024, 3.020558, 3.093847, 3.125291,
            3.152887, 3.169079, 3.184951, 3.185920, 3.190922, 3.194190, 3.198784,
            3.199876, 3.202327, 3.208775, 3.207270, 3.202064, 3.201612)

q95.BMB = c(3.194889, 3.125052, 3.193247, 3.247424, 3.327154, 3.394933, 3.430391,
            3.449532, 3.468933, 3.478339, 3.482464, 3.486812, 3.488204, 3.495350,
            3.495400, 3.499534, 3.502671, 3.499327, 3.496374, 3.495227)

q97.5.BMB = c(3.436780, 3.345678, 3.410049, 3.454669, 3.609811, 3.697073, 3.745145,
              3.767786, 3.781051, 3.787090, 3.782957, 3.781051, 3.782957, 3.786452,
              3.786938, 3.789178, 3.792069, 3.787090, 3.787587, 3.785196)

# Exponential bootstrap
q50.med.EMB = c(1.698009, 1.747480, 1.797101, 1.850636, 1.921433, 1.968298, 2.004880,
                2.034350, 2.052124, 2.071072, 2.086350, 2.094786, 2.103506, 2.111070,
                2.113707, 2.118544, 2.123984, 2.126319, 2.127185, 2.129520)

q90.EMB = c(2.547952, 2.676471, 2.773629, 2.871313, 2.967674, 3.054859, 3.084846,
            3.132493, 3.140999, 3.151648, 3.162487, 3.164651, 3.164418, 3.179678,
            3.181276, 3.180071, 3.180714, 3.180714, 3.175093, 3.178362)

q95.EMB = c(2.920454, 3.002033, 3.083233, 3.183117, 3.318063, 3.408047, 3.440737,
            3.471478, 3.475068, 3.480855, 3.487812, 3.488501, 3.484810, 3.488501,
            3.496750, 3.488792, 3.490025, 3.489197, 3.484864, 3.488792)

q97.5.EMB = c(3.173327, 3.335200, 3.412157, 3.459616, 3.628592, 3.706298, 3.734808,
              3.766433, 3.767562, 3.761745, 3.771280, 3.772386, 3.763792, 3.767562,
              3.771280, 3.760191, 3.761745, 3.763792, 3.758779, 3.761745)

plot(q50.med.WMB, type = "l", col = "red", ylim = c(1.5,4.3))
lines(q90.WMB, col = "blue")
lines(q95.WMB, col = "black")
lines(q97.5.WMB, col = "green")
lines(q50.med.NMB, col = "red", lty = 2)
lines(q90.NMB, col = "blue", lty = 2)
lines(q95.NMB, col = "black", lty = 2)
lines(q97.5.NMB, col = "green", lty = 2)
lines(q50.med.EMB, col = "red", lty = 3)
lines(q90.EMB, col = "blue", lty = 3)
lines(q95.EMB, col = "black", lty = 3)
lines(q97.5.EMB, col = "green", lty = 3)
lines(q50.med.BMB, col = "red", lty = 4)
lines(q90.BMB, col = "blue", lty = 4)
lines(q95.BMB, col = "black", lty = 4)
lines(q97.5.BMB, col = "green", lty = 4)


lines(q50.med, col = "red", lty = 2)
lines(q90, col = "blue", lty = 2)
lines(q95, col = "black", lty = 2)
lines(q97.5, col = "green", lty = 2)

