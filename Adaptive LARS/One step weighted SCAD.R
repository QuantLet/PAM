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
  indicator  = ifelse(beta <= lambda, 1, 0) 
  tmp        = numeric(0)
  for (j in 1:n.par){
    tmp[j]   = max(a * lambda - beta[j], 0)
  }
  pen.value  = lambda * (indicator 
                        + tmp/((a - 1) * lambda) * indicator) 
  pen.value
}

# Define function for submatrix X_U
XU.fct = function(X, beta, lambda, a){
  weights = w.fct(beta)
  pen  = pen.prime(beta * weights, lambda, a) # Evaluate penalty function with weights
  if (min(pen) == 0){
    U.index = which(pen == 0)
    X.U     = X[, U.index]
  } else {
    X.U     = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  }
  X.U
}

# Define function for submatrix X_V
XV.fct = function(X, beta, lambda, a){
  weights = w.fct(beta)
  pen     = pen.prime(beta * weights, lambda, a) # Evaluate penalty function with weights
  if (max(pen) == 0){
    X.V      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  } else {
    V.index  = which(pen != 0)
    pen.new  = pen[V.index]
    beta.new = beta[V.index]
    wgth.new = weights[V.index]
    X.V      = X[, V.index]
    for (j in 1:length(V.index)){
      X.V[, j] = X.V[, j] * lambda / (pen.new[j] * wgth.new[j]) 
    }
  }
  X.V
}

# Define function for finding real X.V coefficients
V.coeff.fct = function(coeff, beta, lambda, a){
  weights   = w.fct(beta)
  pen       = pen.prime(beta * weights, lambda, a) # Evaluate penalty function with weights
  V.index   = which(pen != 0)
  pen.new   = pen[V.index]
  wgth.new  = weights[V.index]
  coeff.new = coeff[V.index]
  V.coeff   = numeric(0)
  for (j in 1:length(V.index)){
    V.coeff[j] = coeff.new[j] * lambda / (pen.new[j] * wgth.new[j]) # Transform coeff's back
  }
  V.coeff
}

# Define function determining weigths
w.fct = function(beta){
  wght = numeric(0)
  for (j in 1:length(beta)){
    wght[j] = ifelse(beta[j] == 0, n.obs ^ (1/4), 1/sqrt(abs(beta[j]))) # Compute weights cf. slide 2-6
  }
  wght
}

# Function for permutation of inital beta fit (so it correspondents to X.U and X.V)
b.fct = function(beta, lambda, a){
  weights   = w.fct(beta)
  pen       = pen.prime(beta * weights, lambda, a) # Evaluate penalty function with weights
  V.index   = which(pen != 0)
  U.index   = which(pen == 0)
  beta.perm = c(beta[U.index], beta[V.index])
  beta.perm
}

# Function for permutation of columns of X
X.fct = function(X, beta, lambda, a){
  weights   = w.fct(beta)
  pen       = pen.prime(beta * weights, lambda, a) # Evaluate penalty function with weights
  V.index   = which(pen != 0)
  U.index   = which(pen == 0)
  X.perm    = cbind(X[, U.index], X[, V.index])
  X.perm
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
lambda.grid = seq(50/sqrt(n.obs), 100/sqrt(n.obs), 1/sqrt(n.obs))   # Define grid of lambdas
a           = 3.7             # Recommended value of parameter a for SCAD
s           = 1               # For now use only 1 simulated scenario

# Simulation of the error term
eps = list()
set.seed(seed.eps)
for (i in 1:n.sim){
eps[[i]] = rnorm(n.obs, mean = 0, sd = sd.eps)
}

# Simulation of true beta coefficients  
tmp1  = rep(1, 5)
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
Adapt.SCAD = function(X, Y, a){
  beta.tmp    = list()
  design.tmp  = list()
  weights.tmp = list()
  for (l in 1:length(lambda.grid)){
    lambda    <<- lambda.grid[l]
    # 0. STEP: Find initial values of coefficients - OLS fit
    beta.fit  = solve(t(X[[s]]) %*% X[[s]]) %*% t(X[[s]]) %*% Y[[s]]
    mu0.hat   = X[[s]] %*% beta.fit # Initial fit
  
    # 1. STEP: Create working data
    # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
    # 1.a)
    # Y.star = mu0.hat # Create response according to the inital fit
    Y.star = Y[[s]] # Use response as it was observed/simulated
    X.star = X[[s]]
  
    # Repeat until convergence
    conv       = Inf
    while (conv > 1 * 10 ^ -5) {
      # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
      # 1.b)
      X.U      = XU.fct(X.star, beta.fit, lambda, a)
      X.V      = XV.fct(X.star, beta.fit, lambda, a)
  
      # 1.c)
      if (X.V[1, 1] != -Inf){
        if (X.U[1, 1] != -Inf){
          Y.2star      = Y.star - X.U %*% solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
          X.V2star     = X.V - X.U %*% solve(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
        } else {
          Y.2star      = Y.star
          X.V2star     = X.V
        }
   
        X.V2starnorm = scale(X.V2star, TRUE, TRUE) # Normalize X_V matrix
        one          = as.vector(rep(1, n.obs))
        V.norm       = as.vector(sqrt(drop(one %*% (X.V2star^2))/(n.obs - 1))) # Compute the norm
    
        # 2. STEP: Lasso estimation with LARS using CV
        object     = glmnet(X.V2starnorm, Y.2star, family = "gaussian", alpha = 1, lambda = lambda.grid[l],
                            standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
    
        coefftmp   = as.vector(object$beta) # Extract coefficients from the fit
        V.coefftmp = coefftmp/V.norm        # Transform them back to the values before normalization
 
        # 3. STEP: Coefficients associated with X_U and X_V
        if (X.U[1, 1] != -Inf){
          U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coefftmp)) 
        } else {
          U.coeff = numeric(0)
        }
    
        V.coeff = V.coeff.fct(V.coefftmp, beta.fit, lambda, a)
    
        beta.fit.per  = b.fct(beta.fit, lambda, a)                # Permutated initial coefficients
        X.per         = X.fct(X.star, beta.fit, lambda, a)        # Permutated design matrix
        coeff         = c(U.coeff, V.coeff)            # Fitted coefficients
        conv          = max(abs(beta.fit.per - coeff)) # Difference between 2 steps of iteration
    
        act.set       = sum(coeff != 0)         # No of nonzero coefficients
    
        beta.fit      = coeff                   # Define inital coefficients for next iteration step
        X.star        = X.per                   # Define permutated design for next iteration step
      } else if (X.V[1, 1] == -Inf){
        # 2. STEP: Skip - no V coefficients
        # 3. STEP: Coefficients associated with X_U and X_V
        U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
        V.coeff = numeric(0)
    
        beta.fit.per  = b.fct(beta.fit, lambda, a)                # Permutated initial coefficients
        X.per         = X.fct(X.star, beta.fit, lambda, a)        # Permutated design matrix
        coeff         = c(U.coeff, V.coeff)            # Fitted coefficients
        conv          = max(abs(beta.fit.per - coeff)) # Difference between 2 steps of iteration
        weights       = w.fct(beta.fit.per)
        
        act.set       = sum(coeff != 0)         # No of nonzero coefficients
        beta.fit      = coeff                   # Define inital coefficients for next iteration step
        X.star        = X.per                   # Define permutated design for next iteration step        
      }
    }
    beta.tmp[[l]]    = beta.fit
    weights.tmp[[l]] = weights
    design.tmp[[l]]  = X.star
  }
  loglik = numeric(0)
  for (l in 1:length(lambda.grid)){
    loglik[l] = (n.obs/2 * log(2 * pi * 1) - 1/(2 * 1) 
                * t(Y[[s]] - design.tmp[[l]] %*% beta.tmp[[l]])
                %*% (Y[[s]] - design.tmp[[l]] %*% beta.tmp[[l]])
                - (n.obs * pen.prime(beta.tmp[[l]] * weights.tmp[[l]], lambda.grid[l], a)
                %*% abs(beta.tmp[[l]] * weights.tmp[[l]])))
  }
  index         = which.max(loglik)
  values        = list(beta.tmp[[index]], design.tmp[[index]], lambda.grid[index])
  names(values) = c("beta", "design", "lambda")
  return(values)
}

tmp = Adapt.SCAD(X, Y, a)

