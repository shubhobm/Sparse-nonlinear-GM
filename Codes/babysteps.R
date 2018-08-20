rm(list=ls())
setwd('c:/Study/My projects/Sparse-nonlinear-GM/Codes/')
setwd('c:/Study/Sparse-nonlinear-GM/Codes/')

source('HLGM.R')
source('HLGM-bin.R')
source('Generator_new2.R')

source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

analyze = function(Obj){
  iter = length(Obj$B.list)
  eval.list = lapply(1:iter, function(x)
    eval.Omega(Obj$Omega.list[[x]], Obj$B.list[[x]], Omega0))
  z = lapply(eval.list, function(x) x$Metrics)
  df = data.frame(cbind(round(matrix(unlist(z),ncol=5, byrow=T),2),
                          Obj$Obj))
  names(df) = c("TP","TN","TP1","TN1","Sp","Obj")
  df
}

trace.plot = function(Obj){
  par(mfrow=c(2,3))
  K = length(Obj$var.list[[1]][[1]])
  
  for(vi in 1:6){
    plot(as.numeric(lapply(Obj$var.list, function(x) x[[vi]][[1]])),
         type='l', ylim=c(0,1.5))
    for(k in 2:K){
      lines(as.numeric(lapply(Obj$var.list, function(x) x[[vi]][[k]])), lty=k)
    }
  }
  par(mfrow=c(1,1))
}

##### Generate data
n = 1e3
p = 20
q = 10
set.seed(11192017)

## looping function
loopfun = function(rep){
  # set.seed(rep*11222017)
  # Omega0 = ThetaMat(p, SNR=1)[[1]]
  # E = mvrnorm(n, mu=matrix(rep(0,p),ncol=1), Sigma=solve(Omega0))
  # B0 = ThetaXYMat(q,p)
  # Omega.z = ThetaMat(q)[[1]]
  # Z0 = mvrnorm(n, mu=matrix(rep(0,q),ncol=1), Sigma=solve(Omega.z))
  # X = Z0 %*% B0 + E

  # Generate data: p vars from a GGM
  # then a quadratic polynomial transformation with random coefs
  set.seed(rep*11222018)
  Omega0 = ThetaMat(p, SNR=1)[[1]]
  E = mvrnorm(n, mu=matrix(rep(0,p),ncol=1), Sigma=solve(Omega0))
  X = matrix(0,n,p)
  for(i in 1:p){
    iX = cbind(1,E[,i],E[,i]^2)
    X[,i] = iX %*% runif(3,-2,2)
  }

  # q=5
  system.time(Objq5 <- hlgm2.bin(X, q=5, K=3, momentum=0,maxit=50))
  a = analyze(Objq5)
  a
  plot(a$Obj, type='l')
  trace.plot(Objq5)
  
  # no momentum
  system.time(Objq10 <- hlgm2.bin(X, q=10, K=5, momentum=0,maxit=50))
  a = analyze(Objq10)
  a
  plot(a$Obj, type='l')
  trace.plot(Objq10)
  
  out200p30 = list(Objq5,Objq10)
  save(out200p30, file="outn200p30.Rda")
  
  # no momentum
  system.time(Objq20 <- hlgm2(X, q=20, momentum=0,maxit=20))
  a = analyze(Objq20)
  plot(a$Obj, type='l')
  trace.plot(Objq20)
  
  analyze(list200[[1]])
  analyze(list200[[2]])
  

  
  # q=20
  system.time(Objq20 <- hlgm2(X, q=20, momentum=0,maxit=30))
  iterq20 = length(Objq20$B.list)
  evalq20.list = lapply(1:iterq20, function(x)
    eval.Omega(Objq20$Omega.list[[x]], Objq20$B.list[[x]], Omega0))
  z = lapply(evalq20.list, function(x) x$Metrics)
  dfq20 = data.frame(cbind(round(matrix(unlist(z),ncol=5, byrow=T),2),
                        Objq20$Obj))
  names(dfq20) = c("TP","TN","TP1","TN1","Sp","Obj")
  dfq20
  
  # Take best model, apply HLGM to residuals
  set.seed(08012018)
  which1 = which.min(Obj11$Obj[1:iter11])
  vars = Obj11$var.list[[which1]]
  names(vars) = c("pi.m","mu.m","sig.m","pi.s","mu.s","sig.s")
  B1 = Obj11$B.list[[which1]]
  Adj1 = eval11.list[[which1]]$Adj
  Z1 =  with(vars, GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s))

  Obj2 = hlgm2(Z1, q=5, momentum=0)
  iter2 = length(Obj2$B.list)
  eval2.list = lapply(1:iter2, function(x)
    eval.Omega(Obj2$Omega.list[[x]], Obj2$B.list[[x]], Omega0,
               Obj11$Omega.list[[which1]]))
  z = lapply(eval2.list, function(x) x$Metrics)
  df = data.frame(cbind(round(matrix(unlist(z),ncol=5, byrow=T),2),
                        Obj2$Obj))
  names(df) = c("TP","TN","TP1","TN1","Sp","Obj")
  df
  
  # Take best model, apply HLGM to residuals
  which2 = which.min(Obj2$Obj[1:iter2])
  var2 = Obj2$var.list[[iter2]]
  B2 = Obj2$B.list[[iter2]]
  Adj2 = eval2.list[[which2]]$Adj
  Adj12 = Adj1 + Adj2
  sum(Adj12!=0 & Omega0!=0)/sum(Omega0!=0)
  sum(Adj12==0 & Omega0==0)/sum(Omega0==0)
  
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
