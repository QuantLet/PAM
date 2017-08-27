# ALASSO from lars with BIC stopping rule
Adapt.LASSO = function(X, Y, n.obs){ 
  
  # 0. STEP: Find initial values of coefficients - OLS fit
  beta.ols  = solve(t(X) %*% X) %*% t(X) %*% Y
  weights   = 1 / abs(beta.ols)
  
  # 1. STEP: Define corrected design matrix X.2star
  X.2star = X
  for (j in 1:length(beta.ols)){
    X.2star[, j] = X[, j] / weights[j]
  }
  
  # 2. STEP: Using LARS algorithm solve beta.2star with GCV
  object      = lars(X.2star, Y, type = "lasso", normalize = TRUE, intercept = FALSE)
  gcv         = n.obs^(-1) * as.vector(object$RSS) / (1 - (as.vector(object$df)/n.obs))
  step.gcv    = which.min(gcv)
  object.pred = predict.lars(object, X.2star, s = step.gcv, type = "coef", 
                             mode = "step")
  beta.2star  = object.pred$coefficients
  
  # 3. STEP: Find final coeffcients  
  beta.final   = beta.2star / weights
  act.set      = object.pred$s
  lambda       = object$lambda[step.gcv] 
 
  values = list(beta.final, act.set, lambda, weights)
  names(values) = c("beta.final", "act.set", "lambda", "weights")
  return(values)
}

Adapt.LASSO.MB = function(X, Y, n.obs, lambda.value, multipl){ 
  
  Y = Y.tmp * sqrt(multipl)
  X = X.tmp * sqrt(multipl)
  
  # 0. STEP: Find initial values of coefficients - OLS fit
  beta.ols  = solve(t(X) %*% X) %*% t(X) %*% Y
  weights   = 1 / abs(beta.ols)
  
  # 1. STEP: Define corrected design matrix X.2star
  X.2star = X
  for (j in 1:length(beta.ols)){
    X.2star[, j] = X[, j] / weights[j]
  }
  
  # 2. STEP: Solve for beta.2star with fixed lambda
  object = glmnet(X.2star, Y, family = "gaussian", alpha = 1, lambda = lambda.value,
                  standardize = TRUE, intercept = FALSE)
  beta.2star  = object$beta
  
  # 3. STEP: Find final coeffcients  
  beta.final   = beta.2star / weights
  
  values = list(beta.final, weights)
  names(values) = c("beta", "weights")
  return(values)
}