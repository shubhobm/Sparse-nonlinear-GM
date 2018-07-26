rm(list=ls())
setwd('c:/Study/My projects/Sparse-nonlinear-GM/Codes/')
source('HLGM.R')
source('Generator_new2.R')

library(glasso)
library(parallel)

##### Generate data
n = 100
p = 20
q = 20
set.seed(11192017)

## looping function
loopfun = function(rep){
  # set.seed(rep*11222017)
  # Omega0 = ThetaMat(p)[[1]]
  # E = mvrnorm(n, mu=matrix(rep(0,p),ncol=1), Sigma=solve(Omega0))
  # B0 = ThetaXYMat(p,q)
  # Omega.z = ThetaMat(q)[[1]]
  # Z0 = mvrnorm(n, mu=matrix(rep(0,q),ncol=1), Sigma=solve(Omega0))
  # X = Z0 %*% B0 + E

  # Generate data: p vars from a GGM
  # then a quadratic polynomial transformation with random coefs
  set.seed(rep*11222017)
  Omega0 = ThetaMat(p, SNR=3)[[1]]
  E = mvrnorm(n, mu=matrix(rep(0,p),ncol=1), Sigma=solve(Omega0))
  X = matrix(0,n,p)
  for(i in 1:p){
    iX = cbind(1,E[,i],E[,i]^2)
    X[,i] = iX %*% runif(3,-2,2)
  }
  
  # Apply HLGM
  Obj1 = hlgm1(X)
  iter1 = length(Obj1$B.list)
  z = lapply(1:iter1, function(x)
    eval.Omega(Obj1$Omega.list[[x]], Obj1$B.list[[x]], Omega0))
  round(matrix(unlist(z),ncol=5, byrow=T),2)
  
  # Take best model, apply HLGM to residuals
  which1 = which.min(Obj1$Obj[1:iter1])
  vars = Obj1$var.list[[which1]]
  B1 = Obj1$B.list[[which1]]
  Z1 =  with(vars, GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s))
  E1 = X - Z1 %*% B1
  Obj2 = hlgm(E1)
  iter2 = length(Obj2$B.list)
  z = lapply(1:iter2, function(x)
    eval.Omega(Obj2$Omega.list[[x]], Obj2$B.list[[x]], Omega0))
  round(matrix(unlist(z),ncol=5, byrow=T),2)
  
  # Take best model, apply HLGM to residuals
  which2 = which.min(Obj2$Obj[1:iter2])
  var2 = Obj2$var.list[[iter2]]
  B2 = Obj2$B.list[[iter2]]
  Z2 =  with(var2, GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s))
  E2 = E1 - Z2 %*% B2
  Obj3 = hlgm(E2)
  z = lapply(1:length(Obj3$B.list), function(x)
    eval.Omega(Obj3$Omega.list[[x]], Obj3$B.list[[x]], Omega0))
  round(matrix(unlist(z),ncol=5, byrow=T),2)

  ## get all models
  loopfun1 = function(m){
    jmle(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
         lambda = lambda.vec[m],
         gamma = sqrt(log(q)/n) * seq(1, 0.2, -0.2),
         init.option=1, tol=1e-3)
  }
  system.time(
    model.list <- mclapply(1:nlambda, loopfun1, mc.cores=min(detectCores(),nlambda))
  )
  # system.time(
  #   model.list <- lapply(1:nlambda, loopfun1)
  # )
  
  ## calculate BIC

  ## all model evaluations

}

eval.list = lapply(1:100, loopfun)
save(eval.list, file="eval_list.Rda")