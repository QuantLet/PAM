pen.prime = function(beta, lambda, a){
  indicator  = ifelse(abs(beta) <= lambda, 1, 0) 
  tmp        = numeric(0)
  for (j in 1:n.par){
    tmp[j]   = max(a * lambda - abs(beta[j]), 0)
  }
  pen.value  = (lambda * indicator) + ((tmp/(a - 1)) * (1 - indicator))
  pen.value
}

lambda = 2
a      = 3.7
beta   = 0.5
n.par  = 1

scad.fct = function(beta, lambda, a){
  scad = ifelse(abs(beta) <= lambda, lambda * abs(beta),
                ifelse(abs(beta) > lambda && abs(beta) <= a * lambda, 
                       -(abs(beta) ^ 2 - (2 * a * lambda * abs(beta)) + lambda ^ 2)/(2 * (a - 1)), 
                       ((a + 1) * (lambda ^ 2))/2))
}

lqa.fct =  function(beta, lambda, a, beta.new){
  lqa   = scad.fct(beta, lambda, a) + (1/2) * (pen.prime(beta, lambda, a)/abs(beta)) * (beta.new ^ 2 - beta ^ 2)
}

lla.fct = function(beta, lambda, a, beta.new){
  lla   = scad.fct(beta, lambda, a) + pen.prime(beta, lambda, a) * (abs(beta.new) - abs(beta))
}

pam.fct = function(beta, lambda, a, beta.new){
  pam   = scad.fct(beta, lambda, a) + pen.prime(beta, lambda, a) * (abs(beta.new) - abs(beta)) / sqrt(beta)
}


x = seq(-10,10,0.01)
scad.val = numeric(0)
lqa.val = numeric(0)
lla.val = numeric(0)
pam.val = numeric(0)
for (i in 1:length(x)){
  scad.val[i] = scad.fct(x[i], lambda, a)
  lqa.val[i] = lqa.fct(beta, lambda, a, x[i])
  lla.val[i] = lla.fct(beta, lambda, a, x[i])
  pam.val[i] = pam.fct(beta, lambda, a, x[i])
}

plot(scad.val, type = "l", ylim = c(0,20))
lines(lqa.val, lty = 2)
lines(lla.val, lty = 3)
lines(pam.val, col = "red")
