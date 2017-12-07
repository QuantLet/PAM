# -------------------------------------------------------------------------------
# One-step SCAD with BIC as a lambda choosing criterion
# -------------------------------------------------------------------------------

# Define function for computing derivative of the SCAD penalty
pen.prime = function(beta, lambda, a){
  indicator  = ifelse(abs(beta) <= lambda, 1, 0) 
  tmp        = numeric(0)
  for (jpp in 1:n.par){
    tmp[jpp]   = max(a * lambda - abs(beta[jpp]), 0)
  }
  pen.value  = (lambda * indicator) + ((tmp/(a - 1)) * (1 - indicator))
  pen.value
}

# Define function for computing the SCAD penalty
pen.scad = function(beta, lambda, a){
  pen.value = numeric(0)
  for (jps in 1:n.par){
    if (abs(beta[jps]) <= lambda){
      pen.value[jps] = lambda * abs(beta[jps])
    } else if (lambda < abs(beta[jps]) && abs(beta[jps]) <= a * lambda){
      pen.value[jps] = - ((abs(beta[jps])^2 - 2 * a * lambda * abs(beta[jps]) + lambda^2)/
                         (2 * (a - 1)))
    } else {
      pen.value[jps] =  ((a + 1) * lambda^2)/2
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
  pen     = pen.prime(beta, lambda, a)    # Evaluate penalty function
  if (max(pen) == 0){
    X.V      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
    beta.new = numeric(0)
  } else {
    V.index  = which(pen != 0)
    pen.new  = pen[V.index]
    beta.new = beta[V.index]
    X.V      = as.matrix(X[, V.index])
    for (jxv in 1:length(V.index)){
      X.V[, jxv] = X.V[, jxv] * lambda/ (pen.new[jxv])
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
  pen       = pen.prime(beta, lambda, a) # Evaluate penalty function 
  V.index   = which(pen != 0)
  pen.new   = pen[V.index]
  V.coeff   = numeric(0)
  for (jvc in 1:length(V.index)){
    V.coeff[jvc] = coeff[jvc] * lambda / (pen.new[jvc]) # Transform coeff's back
  }
  V.coeff
}

grid.fct = function(n.obs, Y, X){
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm  = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.maxl  = scale(X.cent, FALSE, X.norm)

  Y.mean  = drop(one %*% Y)/n.obs
  Y.maxl  = scale(Y, Y.mean, FALSE)
  
  max.lambda.tmp = numeric(0)
  for (j in 1:dim(X.maxl)[2]){
    max.lambda.tmp = c(max.lambda.tmp, Y.maxl %*% X.maxl[, j])
  }
  max.lambda = max.lambda.tmp[which.max(max.lambda.tmp)]
  
  lambda.grid = c(max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*1000)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*100)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03))
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/100)) # Define grid of lambdas
  lambda.grid
}

grid.fct.T = function(n.obs, Y, X){
  
  n.obs.L = length(Y)
  one     = rep(1, n.obs.L)
  X.mean  = drop(one %*% X)/n.obs.L
  X.cent  = scale(X, X.mean, FALSE)
  X.norm  = sqrt(drop(one %*% (X.cent^2)/n.obs.L))
  X.maxl  = scale(X.cent, FALSE, X.norm)
  
  Y.mean  = drop(one %*% Y)/n.obs.L
  Y.maxl  = scale(Y, Y.mean, FALSE)
  
  max.lambda.tmp = numeric(0)
  for (j in 1:dim(X.maxl)[2]){
    max.lambda.tmp = c(max.lambda.tmp, Y.maxl %*% X.maxl[, j])
  }
  max.lambda = max.lambda.tmp[which.max(max.lambda.tmp)]
  
  lambda.grid = c(max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*1000)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*100)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03))
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/100)) # Define grid of lambdas
  lambda.grid
}

Adapt.SCAD = function(X, Y, a, n.obs){

  # X = X.tmp.T
  # Y = Y.tmp.T
  # n.obs = n.obs.tmp.T
  beta.tmp    = list()
  act.set     = list()
  beta.prev   = list()
  
  lambda.grid = grid.fct(n.obs, Y, X)
 
  # Normalize design X (mean = 0, var = 1) and response Y (mean = 0)
  X.orig  = X
  Y.orig  = Y
 
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)

  X       = X.norm
  Y       = Y.norm
  
  # 0. STEP: Find initial values of coefficients - penalized fit
  
  beta.fit.ols  = as.vector(glmnet(X.norm, Y.norm, family = "gaussian", alpha = 0, lambda = 0,
                               standardize = TRUE, intercept = TRUE)$beta)
  
  beta.0        = beta.fit.ols/X.norm0

  # as.vector(glmnet(X.orig, Y.norm, family = "gaussian", alpha = 0, lambda = 0, # tieto dve sa rovnaju
  #                  standardize = TRUE, intercept = TRUE)$beta)
  # 
  # as.vector(glmnet(X.norm, Y.norm, family = "gaussian", alpha = 0, lambda = 0,
  #                  standardize = TRUE, intercept = TRUE)$beta)/X.norm0
  mu0.hat   = X %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star = as.matrix(X) 
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step

  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  for (llam in 1:length(lambda.grid)){

    lambda1   = lambda.grid[llam]
   
    # 1.b)
    X.U      = XU.fct(X.star, beta.fit.ols, lambda1, a)
    XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda1, a)
    X.V      = XV.tmp$Xmatrix
    V.beta   = XV.tmp$betapar
    ord      = order.fct(X.star, beta.fit.ols, lambda1, a)
    
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
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
      }
        
      # 3. STEP: Coefficients associated with X_U and X_V
      if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
      } else {
          U.coeff = numeric(0)
      }
        
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0           # Define inital coefficients for next iteration step
    } else if (X.V[1, 1] == -Inf){
      # 2. STEP: Skip - no V coefficients
      # 3. STEP: Coefficients associated with X_U and X_V
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
      V.coeff = numeric(0)
        
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0                   # Define inital coefficients for next iteration step          
    }
     
    act.set[[llam]]     = act.set.tmp
    beta.prev[[llam]]   = beta.prev.tmp
    beta.tmp[[llam]]    = beta.fit
  }
  
  bic    = numeric(0)
  for (lbic in 1:length(lambda.grid)){
    bic[lbic] = (log(n.obs) * act.set[[lbic]]/n.obs * max(1, sqrt(n.obs)/n.par)
                    + log((t(Y - X.cent %*% beta.tmp[[lbic]]) %*% (Y - X.cent %*% beta.tmp[[lbic]]))/n.obs))
  }
  index         = which.min(bic)
  beta.fit.prev = beta.prev[[index]]
  beta.fit      = beta.tmp[[index]]
  lambda.fit    = lambda.grid[index]

  values        = list(beta.0, beta.fit, lambda.fit, act.set[[index]], index, Y.mean, bic[which.min(bic)])
  names(values) = c("beta.0", "beta", "lambda", "act.set", "index", "intercept", "bic")
  return(values)
}

# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB = function(X, Y, a, n.obs, lambda.value, multipl){
  
  beta.tmp    = list()
  act.set     = list()
  beta.prev   = list()
  
  lambda.grid = grid.fct(n.obs, Y, X)
  
  Y = Y * sqrt(multipl)
  X = X * sqrt(multipl)

  # Normalize design X (mean = 0, var = 1) and response Y (mean = 0)
  X.orig  = X
  Y.orig  = Y
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)

  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)

  X       = X.norm
  Y       = Y.norm
  
  # 0. STEP: Find initial values of coefficients - penalized fit
  beta.fit.ols  = as.vector(glmnet(X, Y, family = "gaussian", alpha = 0, lambda = 0,
                               standardize = TRUE, intercept = TRUE)$beta)
  
  beta.0        = beta.fit.ols/X.norm0
  
  mu0.hat   = X %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star   = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star   = as.matrix(X) 
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  lambda   = lambda.value
  
  # 1.b)
  X.U      = XU.fct(X.star, beta.fit.ols, lambda, a)
  XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda, a)
  X.V      = XV.tmp$Xmatrix
  V.beta   = XV.tmp$betapar
  ord      = order.fct(X.star, beta.fit.ols, lambda, a)
      
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
      V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda, a)
    }
        
    # 3. STEP: Coefficients associated with X_U and X_V
    if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
    } else {
      U.coeff = numeric(0)
    }
        
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
    act.set       = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev     = beta.fit.ols
    beta.fit      = coeff/X.norm0                   # Define inital coefficients for next iteration step
  } else if (X.V[1, 1] == -Inf){
    # 2. STEP: Skip - no V coefficients
    # 3. STEP: Coefficients associated with X_U and X_V
    U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
    V.coeff = numeric(0)
        
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
        
    act.set       = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev     = beta.fit.ols
    beta.fit      = coeff/X.norm0                   # Define inital coefficients for next iteration step  
  }
    
  lambda.fit    = lambda.value
  
  bic    = (log(n.obs) * act.set/n.obs * max(1, sqrt(n.obs)/n.par)
            + log((t(Y - X.cent %*% beta.fit) %*% (Y - X.cent %*% beta.fit))/n.obs))
  
  values        = list(beta.0, beta.fit, lambda.fit, act.set, Y.mean, bic)
  names(values) = c("beta.0", "beta", "lambda", "act.set", "intercept", "bic")
  return(values)
}

# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB.varlambda = function(X, Y, a, n.obs, multipl){
  
  # X = X.tmp.L
  # Y = Y.tmp.L
  # n.obs = n.obs.tmp.L
  # multipl = multipl.L
  
  beta.tmp    = list()
  act.set     = list()
  beta.prev   = list()
  
  # Define grid of lambdas for data before multiplication by u_i's
  lambda.grid = grid.fct(n.obs, Y, X)
  
  # Multiply data by u_i's
  Y = Y * sqrt(multipl)
  X = X * sqrt(multipl)
  
  # Normalize design X (mean = 0, var = 1) and response Y (mean = 0)
  X.orig  = X
  Y.orig  = Y

  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm
  Y       = Y.norm
  
  # 0. STEP: Find initial values of coefficients - penalized fit
  
  beta.fit.ols  = as.vector(glmnet(X.norm, Y.norm, family = "gaussian", alpha = 0, lambda = 0,
                               standardize = TRUE, intercept = TRUE)$beta)
  
  beta.0    = beta.fit.ols/X.norm0
  
  # as.vector(glmnet(X.orig, Y.norm, family = "gaussian", alpha = 0, lambda = 0, # tieto dve sa rovnaju
  #                  standardize = TRUE, intercept = TRUE)$beta)
  # 
  # as.vector(glmnet(X.norm, Y.norm, family = "gaussian", alpha = 0, lambda = 0,
  #                  standardize = TRUE, intercept = TRUE)$beta)/X.norm0
  mu0.hat   = X %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star = as.matrix(X) 
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  for (llam in 1:length(lambda.grid)){
    
    lambda1   = lambda.grid[llam]
    
    # 1.b)
    X.U      = XU.fct(X.star, beta.fit.ols, lambda1, a)
    XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda1, a)
    X.V      = XV.tmp$Xmatrix
    V.beta   = XV.tmp$betapar
    ord      = order.fct(X.star, beta.fit.ols, lambda1, a)
    
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
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
      }
      
      # 3. STEP: Coefficients associated with X_U and X_V
      if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
      } else {
        U.coeff = numeric(0)
      }
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0           # Define inital coefficients for next iteration step
    } else if (X.V[1, 1] == -Inf){
      # 2. STEP: Skip - no V coefficients
      # 3. STEP: Coefficients associated with X_U and X_V
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
      V.coeff = numeric(0)
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0                   # Define inital coefficients for next iteration step          
    }
    
    act.set[[llam]]     = act.set.tmp
    beta.prev[[llam]]   = beta.prev.tmp
    beta.tmp[[llam]]    = beta.fit
  }
  
  bic    = numeric(0)
  for (lbic in 1:length(lambda.grid)){
    bic[lbic] = (log(n.obs) * act.set[[lbic]]/n.obs * max(1, sqrt(n.obs)/n.par) 
                 + log((t(Y - X.cent %*% beta.tmp[[lbic]]) %*% (Y - X.cent %*% beta.tmp[[lbic]]))/n.obs))
  }
  
  index         = which.min(bic)
  beta.fit.prev = beta.prev[[index]]
  beta.fit      = beta.tmp[[index]]
  lambda.fit    = lambda.grid[index]
  
  values        = list(beta.0, beta.fit, lambda.fit, act.set[[index]], index, Y.mean, bic[which.min(bic)])
  names(values) = c("beta.0", "beta", "lambda", "act.set", "index", "intercept", "bic")
  return(values)
}

# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB.shift = function(X.tmp.L, Y.tmp.L, X.tmp.R, Y.tmp.R, a, lambda.L, lambda.R, 
                               multipl.L, multipl.R, beta.L, beta.R){
  
  # lambda.L = lambda.tmp.L
  # lambda.R = lambda.tmp.R
  # beta.L = beta.tmp.L
  # beta.R = beta.tmp.R
  
  beta.tmp    = list()
  act.set     = list()
  beta.prev   = list()
  
  Y.L     = Y.tmp.L * sqrt(multipl.L)
  X.L     = X.tmp.L * sqrt(multipl.L)
  n.obs.L = length(Y.L)
  
  Y.R     = Y.tmp.R * sqrt(multipl.R)
  X.R     = X.tmp.R * sqrt(multipl.R)
  n.obs.R = length(Y.R)
  
  # Normalize design X.T (mean = 0, var = 1) and response Y.T (mean = 0)
  Y.orig.T  = c(Y.L, Y.R)
  X.orig.T  = rbind(X.L, X.R)
  
  n.obs.T   = n.obs.L + n.obs.R
  one.T     = rep(1, n.obs.T)
  X.mean.T  = drop(one.T %*% X.orig.T)/n.obs.T
  X.cent.T  = scale(X.orig.T, X.mean.T, FALSE)
  X.norm0.T = sqrt(drop(one.T %*% (X.cent.T^2)/n.obs.T))
  X.norm.T  = scale(X.cent.T, FALSE, X.norm0.T)
 
  # Y.mean.T  = drop(one.T %*% Y.orig.T)/n.obs.T
  # Y.norm.T  = scale(Y.orig.T, Y.mean.T, FALSE)
  
  X.T       = X.norm.T
  
  Y.T       = numeric(0)
  Y.T[1:n.obs.L] = Y.orig.T[1:n.obs.L]
  Y.T[(n.obs.L + 1):n.obs.T] = Y.orig.T[(n.obs.L + 1):n.obs.T] - (X.cent.T[(n.obs.L + 1):n.obs.T, ] %*% (beta.R - beta.L))
  
  Y.mean.T  = drop(one.T %*% Y.T)/n.obs.T
  Y.norm.T  = scale(Y.T, Y.mean.T, FALSE)
  
  beta.fit0 = as.vector(glmnet(X.cent.T, Y.norm.T, family = "gaussian", alpha = 0, lambda = 0,
                                   standardize = TRUE, intercept = TRUE)$beta)
  
  # Rescale beta.fit by X.norm0.T
  beta.fit.ols  = beta.fit0 * X.norm0.T
  
  beta.0    = beta.fit0
  
  mu0.hat   = X.T %*% beta.fit.ols # Initial fit

  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star    = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star    = X.T
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  
  index      = which(grid.fct(length(Y.L), Y.tmp.L, X.tmp.L) == lambda.L)
  lambda1    = grid.fct.T(length(Y.star), Y.tmp.L, X.tmp.L)[index]
    
  # 1.b)
  X.U        = XU.fct(X.star, beta.fit.ols, lambda1, a)
  XV.tmp     = XV.fct(X.star, beta.fit.ols, lambda1, a)
  X.V        = XV.tmp$Xmatrix
  V.beta     = XV.tmp$betapar
  ord        = order.fct(X.star, beta.fit.ols, lambda1, a)
    
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
      V.coeff.tmp   = ifelse(abs(V.beta) <= 2 * lambda1, 
                            sign(V.beta) * max((abs(V.beta) - lambda1), 0), 
                            ((a - 1) * V.beta - sign(V.beta) * a * lambda1)/(a - 2))
        
      V.coeff       = V.coeff.tmp
        
    } else {
      object        = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda1,
                           standardize = FALSE, intercept = FALSE) # Don't use LARS, lambda is fixed
        
      V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
      V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
    }
      
    # 3. STEP: Coefficients associated with X_U and X_V
    if (X.U[1, 1] != -Inf){
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
    } else {
      U.coeff = numeric(0)
    }
      
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
    act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev.tmp = beta.fit.ols
    beta.fit      = coeff/X.norm0.T         # Define inital coefficients for next iteration step
  } else if (X.V[1, 1] == -Inf){
    # 2. STEP: Skip - no V coefficients
    # 3. STEP: Coefficients associated with X_U and X_V
    U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
    V.coeff = numeric(0)
      
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
    act.set   = sum(coeff != 0)         # No of nonzero coefficients
    beta.prev = beta.fit.ols
    beta.fit  = coeff/X.norm0.T         # Define inital coefficients for next iteration step          
  }
    
  lambda.fit    = lambda1
  
  # bic    = (log(n.obs.T) * act.set/n.obs.T * max(1, sqrt(n.obs.T)/n.par)
  #           + log((t(Y.norm.T - X.cent.T %*% beta.fit) %*% (Y.norm.T - X.cent.T %*% beta.fit))/n.obs.T))

  values        = list(beta.0, beta.fit, lambda.fit, act.set)
  names(values) = c("beta.0", "beta", "lambda", "act.set")
  return(values)
}

# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB.shift.varlambda = function(X.tmp.L, Y.tmp.L, X.tmp.R, Y.tmp.R, a, 
                               multipl.L, multipl.R, beta.L, beta.R){
  
  beta.tmp    = list()
  act.set     = list()
  beta.prev   = list()
  
  Y.L     = Y.tmp.L * sqrt(multipl.L)
  X.L     = X.tmp.L * sqrt(multipl.L)
  n.obs.L = length(Y.L)
  
  Y.R     = Y.tmp.R * sqrt(multipl.R)
  X.R     = X.tmp.R * sqrt(multipl.R)
  n.obs.R = length(Y.R)
  
  # Normalize design X.T (mean = 0, var = 1) and response Y.T (mean = 0)
  Y.orig.T  = c(Y.L, Y.R)
  X.orig.T  = rbind(X.L, X.R)
  
  n.obs.T   = n.obs.L + n.obs.R
  one.T     = rep(1, n.obs.T)
  X.mean.T  = drop(one.T %*% X.orig.T)/n.obs.T
  X.cent.T  = scale(X.orig.T, X.mean.T, FALSE)
  X.norm0.T = sqrt(drop(one.T %*% (X.cent.T^2)/n.obs.T))
  X.norm.T  = scale(X.cent.T, FALSE, X.norm0.T)
  
  # Y.mean.T  = drop(one.T %*% Y.orig.T)/n.obs.T
  # Y.norm.T  = scale(Y.orig.T, Y.mean.T, FALSE)
  
  X.T       = X.norm.T
  
  Y.T       = numeric(0)
  Y.T[1:n.obs.L] = Y.orig.T[1:n.obs.L]
  Y.T[(n.obs.L + 1):n.obs.T] = Y.orig.T[(n.obs.L + 1):n.obs.T] - (X.cent.T[(n.obs.L + 1):n.obs.T, ] %*% (beta.R - beta.L))
  
  Y.mean.T  = drop(one.T %*% Y.T)/n.obs.T
  Y.norm.T  = scale(Y.T, Y.mean.T, FALSE)
  
  beta.fit0 = as.vector(glmnet(X.cent.T, Y.norm.T, family = "gaussian", alpha = 0, lambda = 0,
                               standardize = TRUE, intercept = TRUE)$beta)
  
  lambda.grid    = grid.fct(n.obs.T, Y.T, X.norm.T)
  
  # Rescale beta.fit by X.norm0.T
  beta.fit.ols  = beta.fit0 * X.norm0.T
  
  beta.0    = beta.fit0
  
  mu0.hat   = X.T %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star    = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star    = X.T
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  for (llam in 1:length(lambda.grid)){
    
    lambda1   = lambda.grid[llam]
    
    # 1.b)
    X.U      = XU.fct(X.star, beta.fit.ols, lambda1, a)
    XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda1, a)
    X.V      = XV.tmp$Xmatrix
    V.beta   = XV.tmp$betapar
    ord      = order.fct(X.star, beta.fit.ols, lambda1, a)
    
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
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
      }
      
      # 3. STEP: Coefficients associated with X_U and X_V
      if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
      } else {
        U.coeff = numeric(0)
      }
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0.T           # Define inital coefficients for next iteration step
    } else if (X.V[1, 1] == -Inf){
      # 2. STEP: Skip - no V coefficients
      # 3. STEP: Coefficients associated with X_U and X_V
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
      V.coeff = numeric(0)
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0)         # No of nonzero coefficients
      beta.prev.tmp = beta.fit.ols
      beta.fit      = coeff/X.norm0.T                   # Define inital coefficients for next iteration step          
    }
    
    act.set[[llam]]     = act.set.tmp
    beta.prev[[llam]]   = beta.prev.tmp
    beta.tmp[[llam]]    = beta.fit
  }
  
  bic    = numeric(0)
  for (lbic in 1:length(lambda.grid)){
    bic[lbic] = (log(n.obs.T) * act.set[[lbic]]/n.obs.T * max(1, sqrt(n.obs.T)/n.par)  
                 + log((t(Y.norm.T - X.cent.T %*% beta.tmp[[lbic]]) %*% (Y.norm.T - X.cent.T %*% beta.tmp[[lbic]]))/n.obs.T))
  }
  
  index         = which.min(bic)
  beta.fit.prev = beta.prev[[index]]
  beta.fit      = beta.tmp[[index]]
  lambda.fit    = lambda.grid[index]
  
  values        = list(beta.0, beta.fit, lambda.fit, act.set[[index]], index, Y.mean.T, bic[which.min(bic)])
  names(values) = c("beta.0", "beta", "lambda", "act.set", "index", "intercept", "bic")
  return(values)
}

