# -------------------------------------------------------------------------------
# Multiple-step weighted SCAD penalized regression
# -------------------------------------------------------------------------------

# rm(list = ls(all = TRUE))
# graphics.off()
# 
# setwd("")
#
# libraries = c("MASS", "bbmle", "glmnet")
# lapply(libraries, function(x) if (!(x %in% installed.packages())) {
#   install.packages(x)} )
# lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Define log-likelihood for LASSO (2nd step in the Adapt.SCAD algo)
loglik = function(beta){
  X.tmp1 = X.V2star
  Y.tmp1 = Y.2star
  beta   = ifelse(beta <= 10^(-5), 0, beta)
  loglik1 = -(-(n.obs)/2 * log(2 * pi * sd.eps^2) 
              - (1/(2 * sd.eps^2) 
                 * sum((Y.tmp1 - X.tmp1 %*% beta)^2))
              - n.obs * lambda.grid[l] * sum(abs(beta)))
  loglik1
}

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


Adapt.SCAD = function(X, Y, a, n.obs, max.steps){
  # X = X.tmp
  # Y = Y.tmp
  # n.obs = k*M
  beta.tmp    = list()
  act.set     = list()
  weights.tmp = list()
  beta.prev   = list()
  
  # 0. STEP: Find initial values of coefficients - OLS fit
  beta.fit  = solve(t(X) %*% X) %*% t(X) %*% Y
  mu0.hat   = X %*% beta.fit # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star = mu0.hat   # Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star = X 
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  
  # Repeat until convergence
  conv       = Inf
  step       = 0
  while (conv > 1 * 10 ^ -5) {
    # Compute the best couple (lambda, beta(lambda)) for the current step of the algorithm
    for (l in 1:length(lambda.grid)){
      # l = 1
      lambda   = lambda.grid[l]
      
      # 1.b)
      X.U      = XU.fct(X.star, beta.fit, lambda, a)
      XV.tmp   = XV.fct(X.star, beta.fit, lambda, a)
      X.V      = XV.tmp$Xmatrix
      V.beta   = XV.tmp$betapar
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
        if (dim(X.V)[2] == 1){
          V.coeff.tmp = ifelse(abs(V.beta) <= 2 * lambda, 
                               sign(V.beta) * max((abs(V.beta) - lambda), 0), 
                               ((a - 1) * V.beta - sign(V.beta) * a * lambda)/(a - 2))
          
          V.coeff       = V.coeff.tmp
          
        } else {
          object      = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda.grid[l],
                               standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
          
          V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
          V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit, lambda, a)
        }
        
        # 3. STEP: Coefficients associated with X_U and X_V
        if (X.U[1, 1] != -Inf){
          U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
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
        U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
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
    conv          = 10 ^ (-6) # max(abs(beta.fit.prev - beta.fit)) # Difference between 2 steps of iteration
    step          = step + 1
    if (step > max.steps) {break}
  }
  
  values        = list(beta.fit, lambda.grid[index], step, act.set[[index]], weights.fit, pen.fit, pen.der)
  names(values) = c("beta", "lambda", "no.steps", "act.set", "weights", "penalty", "penaltyprime")
  return(values)
}

# Define function to compute iterative weighted SCAD algorithm, Li&Zou(2008)
Adapt.SCAD.MB = function(X, Y, a, n.obs, max.steps, lambda.value, multipl){
  
  beta.tmp    = list()
  act.set     = list()
  weights.tmp = list()
  beta.prev   = list()
  
  Y = Y.tmp * sqrt(multipl)
  X = X.tmp * sqrt(multipl)
  # 0. STEP: Find initial values of coefficients - OLS fit
  beta.fit  = solve(t(X) %*% X) %*% t(X) %*% Y
  mu0.hat   = X %*% beta.fit # Initial fit
  
  # 1. STEP: Create working data
  # Skipping D_ii, in this case it should be 1 (not 2 as in the One-step SCAD paper)
  # 1.a)
  Y.star = mu0.hat # * sqrt(multipl)# Create response according to the inital fit
  # Y.star = Y # Use response as it was observed/simulated
  X.star = X # * sqrt(multipl)
  # Y.star   = X.star %*% beta.fit # Use this one to change Y in every step
  
  # Repeat until convergence
  conv       = Inf
  step       = 0
  while (conv > 1 * 10 ^ -5) {
    # Compute the best couple (lambda, beta(lambda)) for the current step of the algorithm
    for (l in 1:length(lambda.value)){
      # l = 1
      lambda   = lambda.value[l]
      
      # 1.b)
      X.U      = XU.fct(X.star, beta.fit, lambda, a)
      XV.tmp   = XV.fct(X.star, beta.fit, lambda, a)
      X.V      = XV.tmp$Xmatrix
      V.beta   = XV.tmp$betapar
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
        if (dim(X.V)[2] == 1){
          V.coeff.tmp = ifelse(abs(V.beta) <= 2 * lambda, 
                               sign(V.beta) * max((abs(V.beta) - lambda), 0), 
                               ((a - 1) * V.beta - sign(V.beta) * a * lambda)/(a - 2))
          
          V.coeff       = V.coeff.tmp
          
        } else {
          object      = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, lambda = lambda.grid[l],
                               standardize = TRUE, intercept = FALSE) # Don't use LARS, lambda is fixed
          
          V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
          V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit, lambda, a)
        }
        
        # 3. STEP: Coefficients associated with X_U and X_V
        if (X.U[1, 1] != -Inf){
          U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% (Y.star - X.V %*% V.coeff.tmp)) 
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
        U.coeff = as.vector(solve(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
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
    for (l in 1:length(lambda.value)){
      gcv[l] = n.obs^(-1) * (t(Y - X.star %*% beta.tmp[[l]])
                             %*% (Y - X.star %*% beta.tmp[[l]])) / (1 - (act.set[[l]]/n.obs))
    }
      
    # loglik  = numeric(0)  
    # for (l in 1:length(lambda.grid)){
    #   loglik[l] =  (-(n.obs)/2 * log(2 * pi * sd.eps^2) - (1/(2 * sd.eps^2) 
    #                                                        * t(Y.star - X.star %*% beta.tmp[[l]])
    #                                                        %*% (Y.star - X.star %*% beta.tmp[[l]]))
    #                 - ((n.obs) * sum(pen.prime(beta.tmp[[l]], lambda.grid[l], a))))
    # }
    
    # c = 999
    # loglik  = numeric(0)  
    # for (l in 1:length(lambda.grid)){
    #   loglik[l] =  (-(n.obs)/2 * log(2 * pi * sd.eps^2) * sum(multipliers[[k]][c, ]) - (1/(2 * sd.eps^2) 
    #                   * sum((Y.star - X.star %*% beta.tmp[[l]])^2 * multipliers[[k]][c, ]))
    #                 - ((n.obs) * sum(pen.prime(beta.tmp[[l]], lambda.grid[l], a))))
    # }
    
    index         = which.min(gcv)
    beta.fit.prev = beta.prev[[index]]
    beta.fit      = beta.tmp[[index]]
    lambda.fit    = lambda.value[index]
    weights.fit   = weights.tmp[[index]]
    pen.fit       = pen.scad(beta.fit, lambda.fit, a)
    pen.der       = pen.prime(beta.fit, lambda.fit, a)
    conv          = 10 ^ (-6) # max(abs(beta.fit.prev - beta.fit)) # Difference between 2 steps of iteration
    step          = step + 1
    if (step > max.steps) {break}
  }
  
  values        = list(beta.fit, lambda.value[index], step, act.set[[index]], weights.fit, pen.fit, pen.der)
  names(values) = c("beta", "lambda", "no.steps", "act.set", "weights", "penalty", "penaltyprime")
  return(values)
}

# beta.sim         = numeric(0)
# lambda.sim       = numeric(0)
# steps.sim        = numeric(0)
# act.set.sim      = numeric(0)
# beta.sim.1step   = numeric(0)
# lambda.sim.1step = numeric(0)
# steps.sim.1step  = numeric(0)
# act.set.sim.1step  = numeric(0)
# for (s in 1:n.sim){
#   tmp               = Adapt.SCAD(X[[s]], Y[[s]], a, n.obs)
#   tmp.1step         = SCAD.1step(X[[s]], Y[[s]], a)
#   beta.sim          = cbind(beta.sim, tmp$beta)
#   lambda.sim        = cbind(lambda.sim, tmp$lambda)
#   steps.sim         = cbind(steps.sim, tmp$no.steps)
#   act.set.sim       = cbind(act.set.sim, tmp$act.set)
#   beta.sim.1step    = cbind(beta.sim.1step, tmp.1step$beta)
#   lambda.sim.1step  = cbind(lambda.sim.1step, tmp.1step$lambda)
#   steps.sim.1step   = cbind(steps.sim.1step, tmp.1step$no.steps)
#   act.set.sim.1step = cbind(act.set.sim.1step, tmp.1step$act.set)
# }
# mean.beta.sim       = apply(beta.sim, 1, mean)
# mean.beta.sim.1step = apply(beta.sim, 1, mean)
# diff                = beta.sim - beta.sim.1step
# summary(diff)
# plot(mean.beta.sim, type = "l")
# lines(mean.beta.sim.1step, col = "red")
# beta.sim[,1]