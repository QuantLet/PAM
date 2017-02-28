# -------------------------------------------------------------------------------
# One-step weighted SCAD penalized regression
# -------------------------------------------------------------------------------

rm(list = ls(all = TRUE))
graphics.off()

libraries = c("MASS", "bbmle", "glmnet")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Define function for computing derivative of the SCAD penalty
pen.prime = function(beta, lambda, a){
  indicator  = ifelse(abs(beta) <= lambda, 1, 0) 
  tmp        = numeric(0)
  for (j in 1:n.par){
    tmp[j]   = max(a * lambda - abs(beta[j]), 0)
  }
  pen.value  = (lambda * indicator) + ((tmp/(a - 1)) * (1 - indicator))
  pen.value
}

# Define function for submatrix X_U
XU.fct = function(X, beta, lambda, a){
  pen   = pen.prime(beta, lambda, a)       # Evaluate penalty function
  if (min(pen) == 0){
    U.index = which(pen == 0)
    X.U     = as.matrix(X[, U.index])
  } else {
    X.U     = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  }
  X.U
}

# Define function for submatrix X_V
XV.fct = function(X, beta, lambda, a){
  weights = w.fct(beta)
  pen     = pen.prime(beta, lambda, a)    # Evaluate penalty function
  if (max(pen) == 0){
    X.V      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  } else {
    V.index  = which(pen != 0)
    pen.new  = pen[V.index]
    beta.new = beta[V.index]
    wght.new = weights[V.index]
    X.V      = as.matrix(X[, V.index])
    for (j in 1:length(V.index)){
      X.V[, j] = X.V[, j] * lambda / (pen.new[j] * abs(wght.new[j])) 
    }
  }
  X.V
}

# Define function for columns' order in  X_U and X_V
order.fct = function(X, beta, lambda, a){
  pen     = pen.prime(beta, lambda, a)       # Evaluate penalty function
  order   = seq(1, dim(X)[2], 1)
  U.index = which(pen == 0)
  V.index = which(pen != 0)
  U.order = order[U.index]
  V.order = order[V.index]
  order.new = c(U.order, V.order)
  
  order.new
}

# Define function for finding real X.V coefficients
V.coeff.fct = function(coeff, beta, lambda, a){
  weights   = w.fct(beta)
  pen       = pen.prime(beta, lambda, a) # Evaluate penalty function 
  V.index   = which(pen != 0)
  pen.new   = pen[V.index]
  wght.new  = weights[V.index]
  V.coeff   = numeric(0)
  for (j in 1:length(V.index)){
    V.coeff[j] = coeff[j] * lambda / (pen.new[j] * abs(wght.new[j])) # Transform coeff's back
  }
  V.coeff
}

# Define function determining weigths
w.fct = function(beta){
  wght = numeric(0)
  for (j in 1:length(beta)){
    wght[j] = ifelse(beta[j] == 0, n.obs ^ (1/4), 1/sqrt(abs(beta[j]))) # Compute weights cf. slide 2-6
  }
  return(wght)
}

# Simulation setup
n.obs       = 1000            # No of observations
n.par       = 100             # No of parameters
n.sim       = 100             # No of simulations
w           = 110             # Moving window size 
seed.X      = 2222            # Seed simulation for X
seed.eps    = 2333            # Seed simulation for epsilon
r           = 0.5             # Correlation parameter of X
sd.eps      = 1               # Standard deviation of epsilon
lambda.grid = seq(1/sqrt(n.obs), 30/sqrt(n.obs), 1/sqrt(n.obs))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD
s           = 1               # For now use only 1 simulated scenario

# Simulation of the error term
eps = list()
set.seed(seed.eps)
for (i in 1:n.sim){
eps[[i]] = rnorm(n.obs, mean = 0, sd = sd.eps)
}

# Simulation of true beta coefficients  
tmp1  = c(0.1,0.1,0.3,0.3,0.1,0.1,0.1)#rep(1, 5)
tmp2  = rep(0, (n.par - length(tmp1)))
b     = c(tmp1, tmp2)
beta0 = rep(0, n.par)

# Simulation of the design matrix
mu    = rep(0, n.par)
Sigma = matrix(0, nrow = n.par, ncol = n.par)
for (i in 1:n.par) {
  for (j in 1:n.par) {
    Sigma[i, j] = ifelse(i == j, 1, r^abs(i - j))
  }
}

X = list()
set.seed(seed.X)
for (i in 1:n.sim){
  X[[i]] = mvrnorm(n = n.obs, mu, Sigma)
}  

# Computation of Y
Y = list()
for (i in 1:n.sim){
  Y.tmp = numeric(0)
  for (j in 1:n.obs){
    Y.tmp = c(Y.tmp, b %*% X[[i]][j, ] + eps[[i]][j])
  }
Y[[i]] = Y.tmp
}
l = 1
# Define function to compute iterative weighted SCAD algorithm
Adapt.SCAD.old = function(X, Y, a){
  beta.tmp    = list()
  design.tmp  = list()
  weights.tmp = list()
  beta.prev   = list()
  steps       = numeric(0)
  for (l in 1:length(lambda.grid)){
    lambda    = lambda.grid[l]
    # 0. STEP: Find initial values of coefficients - OLS fit
    beta.fit  = solve(t(X[[s]]) %*% X[[s]]) %*% t(X[[s]]) %*% Y[[s]]
    mu0.hat   = X[[s]] %*% beta.fit # Initial fit
  
    # 1. STEP: Create working data
    # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
    # 1.a)
    Y.star = mu0.hat # Create response according to the inital fit
    # Y.star = Y[[s]] # Use response as it was observed/simulated
    X.star = X[[s]]
  
    # Repeat until convergence
    conv       = Inf
    step       = 0
    while (conv > 1 * 10 ^ -5) {
      # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
      # 1.b)
      X.U      = XU.fct(X.star, beta.fit, lambda, a)
      X.V      = XV.fct(X.star, beta.fit, lambda, a)
      ord      = order.fct(X.star, beta.fit, lambda, a)
      
      # 1.c)
      if (X.V[1, 1] != -Inf){
        if (X.U[1, 1] != -Inf){
          Y.2star      = Y.star - X.U %*% solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
          X.V2star     = X.V - X.U %*% solve(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
        } else {
          Y.2star      = Y.star
          X.V2star     = X.V
        }
    
        # 2. STEP: Lasso estimation using CV
        object      = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda.grid[l],
                             standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
    
        V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
 
        # 3. STEP: Coefficients associated with X_U and X_V
        if (X.U[1, 1] != -Inf){
          U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
        } else {
          U.coeff = numeric(0)
        }
    
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit, lambda, a)
        
        coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
        coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
        conv          = max(abs(beta.fit - coeff)) # Difference between 2 steps of iteration
        weights       = w.fct(beta.fit)
    
        act.set       = sum(coeff != 0)         # No of nonzero coefficients
        beta.prev.tmp = beta.fit
        beta.fit      = coeff                   # Define inital coefficients for next iteration step
        # X.star        = X.per                   # Define permutated design for next iteration step
      } else if (X.V[1, 1] == -Inf){
        # 2. STEP: Skip - no V coefficients
        # 3. STEP: Coefficients associated with X_U and X_V
        U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
        V.coeff = numeric(0)
    
        coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
        coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
        conv          = max(abs(beta.fit - coeff)) # Difference between 2 steps of iteration
        weights       = w.fct(beta.fit)
        
        act.set       = sum(coeff != 0)         # No of nonzero coefficients
        beta.prev.tmp = beta.fit
        beta.fit      = coeff                   # Define inital coefficients for next iteration step          
      }
      step = step + 1
    }
    beta.prev[[l]]   = beta.prev.tmp
    beta.tmp[[l]]    = beta.fit
    weights.tmp[[l]] = weights
    steps[l]         = step
  }
  loglik = numeric(0)
  for (l in 1:length(lambda.grid)){
    loglik[l] = (-n.obs/2 * log(2 * pi * 1) - (1/(2 * 1) 
                * t(Y[[s]] - X.star %*% beta.tmp[[l]])
                %*% (Y[[s]] - X.star %*% beta.tmp[[l]]))
                - (n.obs *t(pen.prime(beta.prev[[l]], lambda.grid[l], a))
                %*% abs(beta.tmp[[l]] * weights.tmp[[l]])))
  }
  index         = which.max(loglik)
  values        = list(beta.tmp[[index]], lambda.grid[index], steps[index])
  names(values) = c("beta", "lambda", "no.steps")
  return(values)
}

tmp.old = Adapt.SCAD.old(X, Y, a)
tmp.old$beta
tmp.old$lambda
lambda.grid
tmp.old$no.steps
