# -------------------------------------------------------------------------------
# Multiple-step weighted SCAD penalized regression
# -------------------------------------------------------------------------------

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

# Define function for computing the SCAD penalty
pen.scad = function(beta, lambda, a){
  pen.value = numeric(0)
  for (j in 1:n.par){
    if (abs(beta[j]) <= lambda){
      pen.value[j] = lambda * abs(beta[j])
    } else if (lambda < abs(beta[j]) && abs(beta[j]) <= a * lambda){
      pen.value[j] = - ((abs(beta[j])^2 - 2 * a * lambda * abs(beta[j]) + lambda^2)/
                         (2 * (a - 1)))
    } else {
      pen.value[j] =  ((a + 1) * lambda^2)/2
    }
  }
  pen.value
}

# Define function for submatrix X_U
XU.fct = function(X, beta, lambda, a){
  pen   = pen.prime(beta, lambda, a)       # Evaluate penalty function
  if (min(pen) == 0){
    U.index  = which(pen == 0)
    X.U      = as.matrix(X[, U.index])
  } else {
    X.U      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  }
  X.U
}

# Define function for submatrix X_V
XV.fct = function(X, beta, lambda, a){
  weights = w.fct(beta)
  pen     = pen.prime(beta, lambda, a)    # Evaluate penalty function
  if (max(pen) == 0){
    X.V      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
    beta.new = numeric(0)
  } else {
    V.index  = which(pen != 0)
    pen.new  = pen[V.index]
    wght.new = weights[V.index]
    beta.new = beta[V.index]
    X.V      = as.matrix(X[, V.index])
    for (j in 1:length(V.index)){
      X.V[, j] = X.V[, j] * lambda/ (pen.new[j] * abs(wght.new[j])) 
    }
  }
  values = list(X.V, beta.new)
  names(values) = c("Xmatrix", "betapar")
  return(values)
}

# Define function for columns' order in  X_U and X_V
order.fct = function(X, beta, lambda, a){
  pen       = pen.prime(beta, lambda, a)       # Evaluate penalty function
  order     = seq(1, dim(X)[2], 1)
  U.index   = which(pen == 0)
  V.index   = which(pen != 0)
  U.order   = order[U.index]
  V.order   = order[V.index]
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


Adapt.SCAD = function(X, Y, a, n.obs){
  # X = X.tmp
  # Y = Y.tmp
  # n.obs = 12
  lambda.grid = c(c(1/(n.obs^seq(0.49, 0.01, -0.03)), seq(2,30,1)/(n.obs^0.01))/100, c(seq(2,30,1)/(n.obs^0.01)/10))     # Define grid of lambdas
  beta.tmp    = list()
  act.set     = list()
  weights.tmp = list()
  beta.prev   = list()
  
  # Normalize columns of X and Y
  # one     = rep(1, n.obs)
  # X.mean  = drop(one %*% X)/n.obs
  # X.cent  = scale(X, X.mean, FALSE)
  # X.norm  = sqrt(drop(one %*% (X.cent^2)))
  # X       = scale(X.cent, FALSE, X.norm)
  # 
  # Y.mean  = drop(one %*% Y)/n.obs
  # Y       = scale(Y, Y.mean, FALSE)
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  for (l in 1:length(lambda.grid)){
    # l = 1
    lambda1   = lambda.grid[l]
    
    # 0. STEP: Find initial values of coefficients - OLS fit
    beta.fit  = as.vector(glmnet(X, Y, family = "gaussian", alpha = 0, lambda = lambda1,
                       standardize = TRUE, intercept = FALSE)$beta)
    mu0.hat   = X %*% beta.fit # Initial fit
    
    # 1. STEP: Create working data
    # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
    # 1.a)
    Y.star = mu0.hat   # Create response according to the inital fit
    # Y.star = Y # Use response as it was observed/simulated
    X.star = as.matrix(X) 
    # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
    
      
    # 1.b)
    X.U      = XU.fct(X.star, beta.fit, lambda1, a)
    XV.tmp   = XV.fct(X.star, beta.fit, lambda1, a)
    X.V      = XV.tmp$Xmatrix
    V.beta   = XV.tmp$betapar
    ord      = order.fct(X.star, beta.fit, lambda1, a)
    
    # 1.c)
    if (X.V[1, 1] != -Inf){
      if (X.U[1, 1] != -Inf){
        Y.2star      = Y.star - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
        X.V2star     = X.V - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
      } else {
        Y.2star      = Y.star
        X.V2star     = X.V
      }
        
    # 2. STEP: Lasso estimation using CV
      if (dim(X.V)[2] == 1){
        V.coeff.tmp = ifelse(abs(V.beta) <= 2 * lambda1, 
                             sign(V.beta) * max((abs(V.beta) - lambda1), 0), 
                             ((a - 1) * V.beta - sign(V.beta) * a * lambda1)/(a - 2))
          
        V.coeff       = V.coeff.tmp
          
      } else {
        object      = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda1,
                             standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
          
        V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit, lambda1, a)
      }
        
      # 3. STEP: Coefficients associated with X_U and X_V
      if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
      } else {
          U.coeff = numeric(0)
      }
        
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      weights       = w.fct(beta.fit)
        
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit
      beta.fit      = coeff                   # Define inital coefficients for next iteration step
    } else if (X.V[1, 1] == -Inf){
      # 2. STEP: Skip - no V coefficients
      # 3. STEP: Coefficients associated with X_U and X_V
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
      V.coeff = numeric(0)
        
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      weights       = w.fct(beta.fit)
        
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit
      beta.fit      = coeff                   # Define inital coefficients for next iteration step          
    }
      
    act.set[[l]]     = act.set.tmp
    beta.prev[[l]]   = beta.prev.tmp
    beta.tmp[[l]]    = beta.fit
    weights.tmp[[l]] = weights
  }
    
  gcv    = numeric(0)
  for (l in 1:length(lambda.grid)){
    gcv[l] = n.obs^(-1) * (t(Y - X.star %*% beta.tmp[[l]])
                           %*% (Y - X.star %*% beta.tmp[[l]])) / (1 - (act.set[[l]]/n.obs))
  }
    
  index         = which.min(gcv)
  beta.fit.prev = beta.prev[[index]]
  beta.fit      = beta.tmp[[index]]
  lambda.fit    = lambda.grid[index]
  weights.fit   = weights.tmp[[index]]
  pen.fit       = pen.scad(beta.fit, lambda.fit, a)
  pen.der       = pen.prime(beta.fit, lambda.fit, a)
  
  values        = list(beta.fit, lambda.fit, act.set[[index]], weights.fit, pen.fit, pen.der)
  names(values) = c("beta", "lambda", "act.set", "weights", "penalty", "penaltyprime")
  return(values)
}


# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB = function(X, Y, a, n.obs, lambda.value, multipl){
  
  beta.tmp    = list()
  act.set     = list()
  weights.tmp = list()
  beta.prev   = list()
  
  Y = Y.tmp * sqrt(multipl)
  X = X.tmp * sqrt(multipl)

  # Normalize columns of X and Y
  # one     = rep(1, n.obs)
  # X.mean  = drop(one %*% X)/n.obs
  # X.cent  = scale(X, X.mean, FALSE)
  # X.norm  = sqrt(drop(one %*% (X.cent^2)))
  # X       = scale(X.cent, FALSE, X.norm)
  # 
  # Y.mean  = drop(one %*% Y)/n.obs
  # Y       = scale(Y, Y.mean, FALSE) 
  
  # 0. STEP: Find initial values of coefficients - OLS fit
  beta.fit  = as.vector(glmnet(X, Y, family = "gaussian", alpha = 0, lambda = lambda.value,
                               standardize = TRUE, intercept = FALSE)$beta)
  mu0.hat   = X %*% beta.fit # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star = as.matrix(X) 
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  lambda   = lambda.value
      
  # 1.b)
  X.U      = XU.fct(X.star, beta.fit, lambda, a)
  XV.tmp   = XV.fct(X.star, beta.fit, lambda, a)
  X.V      = XV.tmp$Xmatrix
  V.beta   = XV.tmp$betapar
  ord      = order.fct(X.star, beta.fit, lambda, a)
      
  # 1.c)
  if (X.V[1, 1] != -Inf){
    if (X.U[1, 1] != -Inf){
      Y.2star      = Y.star - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
      X.V2star     = X.V - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
    } else {
      Y.2star      = Y.star
      X.V2star     = X.V
    }
        
    # 2. STEP: Lasso estimation using CV
    if (dim(X.V)[2] == 1){
      V.coeff.tmp = ifelse(abs(V.beta) <= 2 * lambda, 
                           sign(V.beta) * max((abs(V.beta) - lambda), 0), 
                           ((a - 1) * V.beta - sign(V.beta) * a * lambda)/(a - 2))
          
      V.coeff       = V.coeff.tmp
          
    } else {
      object      = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda.value,
                           standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
          
      V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
      V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit, lambda, a)
    }
        
    # 3. STEP: Coefficients associated with X_U and X_V
    if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
    } else {
      U.coeff = numeric(0)
    }
        
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    weights.tmp   = w.fct(beta.fit)
        
    act.set       = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev     = beta.fit
    beta.fit      = coeff                   # Define inital coefficients for next iteration step
  } else if (X.V[1, 1] == -Inf){
    # 2. STEP: Skip - no V coefficients
    # 3. STEP: Coefficients associated with X_U and X_V
    U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
    V.coeff = numeric(0)
        
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    weights.tmp   = w.fct(beta.fit)
        
    act.set       = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev     = beta.fit
    beta.fit      = coeff                   # Define inital coefficients for next iteration step  
  }

  gcv    = n.obs^(-1) * (t(Y - X.star %*% beta.fit)
                         %*% (Y - X.star %*% beta.fit)) / (1 - (act.set/n.obs))
    
  lambda.fit    = lambda.value
  weights.fit   = weights.tmp
  pen.fit       = pen.scad(beta.fit, lambda.fit, a)
  pen.der       = pen.prime(beta.fit, lambda.fit, a)
  
  values        = list(beta.fit, lambda.fit, act.set, weights.fit, pen.fit, pen.der)
  names(values) = c("beta", "lambda", "act.set", "weights", "penalty", "penaltyprime")
  return(values)
}


