## ****************************************************#
## The Joint Multiple Multi-level Estimation (JMMLE) Method
## ****************************************************#

source('l1ML_Main.R')

## Clip a gradient
clip = function(x, ceil=1e3){
  
  xx = x
  for(i in 1:length(xx)){
    if(abs(x[i]) > ceil){
      xx[i] = sign(x[i])*ceil
    } 
  }
  xx
}

## Calculate objective function
Objfunc = function(X, Z, B, Omega, lambda, gamma,
                   pi.m, mu.m, sig.m, pi.s, mu.s, sig.s){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)
  
  E = X - Z %*% B
  Omega.skeleton = which(Omega!=0, arr.ind=T)
  Obj = 0
  
  # accumulate non-zero entries of the precision matrix
  for(j0 in 1:nrow(Omega.skeleton)){
    
    # calculate quantities
    j = Omega.skeleton[j0,1]
    j1 = Omega.skeleton[j0,2]
    Obj = Obj + Omega[j,j1]*sum(E[,j]*E[,j1])
  }
  Obj = Obj/n - log(det(Omega))
  
  # add penalties
  Obj = Obj + lambda*sum(abs(B)) +
    gamma*(sum(abs(Omega)) - sum(diag(Omega)))/2
  
  # add first KL divergence term
  Obj = Obj + q*sum(pi.m*((mu.m^2 + sig.m^2)/2 - log(sig.m))) +
    q*sum(pi.s*((mu.s^2 + exp(-2*sig.s))/2 + sig.s))
  
  # add second KL divergence term
  phi.m.sq = clip(sum(pi.m/sig.m^2), ceil=1e3)
  phi.s.sq = clip(sum(pi.s/sig.s^2), ceil=1e3)
  Obj = Obj +
    sum(pi.m*((sum(Z^2) - 2*sum(Z)*mu.m + q*(mu.m^2 + phi.m.sq))/sig.m^2/2 +
                q*log(sig.m)) +
          pi.s*((sum(Z^2) - 2*sum(Z)*mu.s + q*(mu.s^2 + phi.s.sq))/sig.s^2/2 +
                  q*log(sig.s)))
  - q*log(phi.m.sq)/2 - q*log(phi.s.sq)/2
  Obj
}

## Evaluate outputs from a model object
eval.Omega = function(Omega, B, Omega0, Adj=NULL){
  # adjacency matrix
  Adj.matrix = Omega
  # add edges to adj matrix from 2nd level connections
  for(j in 1:q){
    j.nodes = which(B[,j]!=0)
    if(length(j.nodes)>1){
      for(j1 in j.nodes){
        for(j2 in j.nodes){
          if(j1 != j2){
            Adj.matrix[j1,j2] = 1
            Adj.matrix[j2,j1] = 1
          }
        }
      }
    }
  }
  Adj.matrix[which(Adj.matrix!=0, arr.ind=T)] = 1
  
  # if Adj matrix is supplied then add
  if(!is.null(Adj)){
    Adj.matrix = Adj.matrix + Adj
  }
  
  return(list(Adj=Adj.matrix,
              Metrics=c((sum(Omega!=0 & Omega0!=0)-p)/(sum(Omega0!=0)-p),
                        sum(Omega==0 & Omega0==0)/sum(Omega0==0),
                        (sum(Adj.matrix!=0 & Omega0!=0)-p)/(sum(Omega0!=0)-p),
                        sum(Adj.matrix==0 & Omega0==0)/sum(Omega0==0),
                        sum(Adj.matrix!=0))))
}

## Main algorithm
hlgm = function(X, q, B_init=NULL, Theta_init=NULL, init.option=1,
                lambda=NULL, gamma=NULL, K=3, eta=1e-4,
                refit.B=TRUE, tol=1e-3, grad.max=1e3,
                maxit=20, eps=1e-6, max.epoch=1e2,VERBOSE=TRUE){
  
  #****************************************************#
  # Define and initialize some quantities
  #****************************************************#
  n = nrow(X)
  p = ncol(X)
  
  ## default values of arguments if they are NULL
  if(is.null(lambda)){
    lambda = sqrt(log(p)/n) * seq(2, 1.25, -.25)
  }
  if(is.null(gamma)){
    gamma = sqrt(log(q)/n) * seq(2, 1.25, -.25)
  }
  nl = length(lambda)
  ng = length(gamma)
  
  #****************************************************#
  # Initialization of iterates
  #****************************************************#
  
  ## initialize variational parameters and generate pseudo-data
  set.seed(07232018)
  pi.m = rep(1/K,K)
  mu.m = runif(K,-1,1)
  sig.m = runif(K,0,1)
  pi.s = pi.m
  mu.s = runif(K,0,1)
  sig.s = runif(K,0,0.1)
  Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)

  ## Initialize B, if initial B and Omega not supplied
  if(init.option==1){
    ## generate Z from initialized variational params,
    ## then fit l1ML to initialize B and Omega
    cat("Initializing B and Omega\n")
    B.array = array(0,c(p,q,nl*ng))
    Omega.array = array(0, c(q,q,nl*ng))
    
    # # randomly initialize until at least one non-zero B is found
    
    bic.vec = rep(0, nl*ng)
    n.tune = 0
    pb = txtProgressBar(0,nl*ng)
    for(a in 1:nl){
      for(b in 1:ng){
        n.tune = n.tune+1
        init.Obj = l1ML_Main(X, Z, lambda=lambda[a], rho=gamma[b],
                             StabilizeTheta=F, VERBOSE=F)
        # cat("lambda=",lambda[a],"gamma=",gamma[b],
        #     sum(abs(init.Obj$B.est)),"done!\n")
        
        # calculate BIC
        B.array[,,n.tune] = init.Obj$B.est
        Omega.array[,,n.tune] = init.Obj$Theta.est
        bic.vec[n.tune] = with(init.Obj, -log(det(Theta.est)) +
                                 sum(diag(crossprod(X - Z %*% B.est) %*% Theta.est))/n +
                                 log(n)/n * ((sum(Theta.est != 0)-q)/2 +
                                               sum(B.est != 0)))
        setTxtProgressBar(pb,n.tune)
      }
    }
    close(pb)
    
    B_init = B.array[,,which.min(bic.vec)]
    # random initialization if all coefs are 0 in B matrix
    if(sum(abs(B_init))==0){
      for(i in 1:p){
        for(j in 1:q){
          B_init[i,j] = rbinom(1,1,1/p)*sample(c(-1,1),1)*runif(1,0.5,1);
        }
      }
    }
    Omega_init = Omega.array[,,which.min(bic.vec) ]
  }
  
  #****************************************************#
  # Alternating algorithm
  #****************************************************#
  
  # initialize
  iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
  updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
  B_new = B_init
  Omega_new = Omega_init
  Normdiff = rep(0, maxit)
  Objval = rep(0, maxit)
  B.list = list()
  Omega.list = list()
  var.list = list()
  
  ## start with the alternating procedure
  cat('-----\n')
  while(!CONVERGE){
    iter = iter + 1
    B_old = B_new
    Omega_old = Omega_new
    Omega.skeleton = which(Omega_old!=0, arr.ind=T)
    cat("Iteration ", iter, ":\n")
    
    ## Update variational parameters using SGD
    cat("Updating variational parameters\n")
    epoch = 0
    CONVERGE.var=FALSE
    while(!CONVERGE.var){
      
      epoch = epoch+1
      eta.epoch = eta*exp(-epoch+1)
      perm = gtools::permute(1:n) # randomly order samples
      pi.m.old = pi.m
      pi.s.old = pi.s
      mu.m.old = mu.m
      mu.s.old = mu.s
      sig.m.old = sig.m
      sig.s.old = sig.s
      
      ## terms corresp to squared error loss
      for(i in perm){ # SGD: go through all samples randomly
        xi = X[i,]
        for(j0 in 1:nrow(Omega.skeleton)){
          
          # calculate quantities
          j = Omega.skeleton[j0,1]
          j1 = Omega.skeleton[j0,2]
          xj.1.bj1 = xi[j]*sum(B_old[j1,])
          bj.t.bj1 = sum(B_old[j,]*B_old[j1,])
          mult1 = -(2*Omega_old[j,j1]*xj.1.bj1)/n
          # mult2 = exp(clip(2*(mu.s.old + sig.s.old^2), ceil=3))*bj.t.bj1
          mult2 = exp(-2*(mu.s + sig.s^2))*bj.t.bj1
          # mult2 = exp(2*(mu.s.old + sig.s.old^2))*(bj.t.bj1 +
          #                                            q/2/nrow(Omega.skeleton))
          # calculatre gradients
          grad.pi.m = mult1*mu.m + (mu.m^2+sig.m^2)*bj.t.bj1
          grad.mu.m = pi.m*(mult1 + 2*mu.m*bj.t.bj1)
          grad.sig.m = 2*sig.m*pi.m*bj.t.bj1
          grad.pi.s = mult2
          grad.mu.s = -2*mult2*pi.s
          grad.sig.s = -4*mult2*pi.s*sig.s
          
          # update
          pi.m = pi.m - eta.epoch*clip(grad.pi.m)
          mu.m = mu.m - eta.epoch*clip(grad.mu.m)
          sig.m = sig.m - eta.epoch*clip(grad.sig.m)
          pi.s = pi.s - eta.epoch*clip(grad.pi.s)
          mu.s = mu.s - eta.epoch*clip(grad.mu.s)
          sig.s = sig.s - eta.epoch*clip(grad.sig.s)
          
          # normalize parameters
          pi.m = abs(pi.m)/sum(abs(pi.m))
          pi.s = abs(pi.s)/sum(abs(pi.s))
          sig.m = abs(sig.m)
          sig.s = abs(sig.s)
        }
      }
      # cat('\t',sig.s,'\n----------\n')
      
      ## add terms corresp to first KL divergence
      mult3 = exp(-2*(mu.s.old + sig.s.old^2))/2
      grad1.pi.m = (mu.m.old^2+sig.m.old^2)/2
      grad1.mu.m = pi.m.old*mu.m.old
      grad1.sig.m = pi.m.old*sig.m.old
      grad1.pi.s = mult3 + mu.s.old
      grad1.mu.s = -2*pi.s.old*mult3 + 1
      grad1.sig.s = -4*pi.s.old*sig.s.old*mult3
      
      # update
      pi.m = pi.m - eta.epoch*clip(grad1.pi.m)
      mu.m = mu.m - eta.epoch*clip(grad1.mu.m)
      sig.m = sig.m - eta.epoch*clip(grad1.sig.m)
      pi.s = pi.s - eta.epoch*clip(grad1.pi.s)
      mu.s = mu.s - eta.epoch*clip(grad1.mu.s)
      sig.s = sig.s - eta.epoch*clip(grad1.sig.s)
      
      ## add terms corresp. to second KL divergence
      # update parameters of r(.)
      phi.m.sq = clip(sum(pi.m/sig.m^2), ceil=1e3)
      phi.s.sq = clip(sum(pi.s/sig.s^2), ceil=1e3)
      
      # calculate gradients
      sumZ2 = sum(Z^2)/n
      sumZ = sum(Z)/n
      grad2.pi.m = (sumZ2 - 2*sumZ*mu.m.old +
                      q*(mu.m.old^2 + phi.m.sq))/2/sig.m.old^2
      grad2.mu.m = pi.m.old*(q*mu.m.old - sumZ)/2/sig.m.old^2
      grad2.sig.m = -pi.m.old*(sumZ2 - 2*sumZ*mu.m.old +
                                 q*(mu.m.old^2 + phi.m.sq))/sig.m.old^3
      grad2.pi.s = (sumZ2 - 2*sumZ*mu.s.old +
                      q*(mu.s.old^2 + phi.s.sq))/2/sig.s.old^2
      grad2.mu.s = pi.s.old*(q*mu.s.old - sumZ)/2/sig.s.old^2
      grad2.sig.s = -pi.s.old*(sumZ2 - 2*sumZ*mu.s.old +
                                 q*(mu.s.old^2 + phi.s.sq))/sig.s.old^3
      # update
      pi.m = pi.m - eta.epoch*clip(grad2.pi.m)
      mu.m = mu.m - eta.epoch*clip(grad2.mu.m)
      sig.m = sig.m - eta.epoch*clip(grad2.sig.m)
      pi.s = pi.s - eta.epoch*clip(grad2.pi.s)
      mu.s = mu.s - eta.epoch*clip(grad2.mu.s)
      sig.s = sig.s - eta.epoch*clip(grad2.sig.s)
      
      # normalize parameters
      pi.m = abs(pi.m)/sum(abs(pi.m))
      pi.s = abs(pi.s)/sum(abs(pi.s))
      sig.m = abs(sig.m)
      sig.s = abs(sig.s)
      
      # check convergence
      var.diff = sqrt(sum((pi.m.old-pi.m)^2)/sum(pi.m^2)) +
        sqrt(sum((pi.s.old-pi.s)^2)/sum(pi.s^2)) +
        sqrt(sum((mu.m.old-mu.m)^2)/sum(mu.m^2)) +
        sqrt(sum((mu.s.old-mu.s)^2)/sum(mu.s^2)) +
        sqrt(sum((sig.m.old-sig.m)^2)/sum(sig.m^2)) +
        sqrt(sum((sig.s.old-sig.s)^2)/sum(sig.s^2))
      cat("Norm_diff =",round(var.diff,4),sig.s,'\n-----\n')
      CONVERGE.var = (var.diff<tol)
      if(epoch==max.epoch){
        cat("Max iterations reached.",'\n')
        break;
      }
    }
    
    ## generate pseudo-data Z from variational distributions
    Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    ## Update generative parameters B, Omega
    cat("Updating B and Omega\n")
    B.array = array(0,c(p,q,nl*ng))
    Omega.array = array(0, c(q,q,nl*ng))
    bic.vec = rep(0, nl*ng)
    n.tune = 0
    pb = txtProgressBar(0, nl*ng)
    for(a in 1:nl){
      for(b in 1:ng){
        n.tune = n.tune+1
        ML.Obj = l1ML_Main(X, Z, lambda=lambda[a], rho=gamma[b],
                           StabilizeTheta=F, VERBOSE=F)
        # cat("lambda=",lambda[a],"gamma=",gamma[b],
        #     sum(abs(ML.Obj$B.est)),"done!\n")
        
        # calculate BIC
        B.array[,,n.tune] = ML.Obj$B.est
        Omega.array[,,n.tune] = ML.Obj$Theta.est
        bic.vec[n.tune] = with(ML.Obj, -log(det(Theta.est)) +
                                 sum(diag(crossprod(X - Z %*% B.est) %*% Theta.est))/n +
                                 log(n)/n * ((sum(Theta.est != 0)-q)/2 +
                                               sum(B.est != 0)))
        setTxtProgressBar(pb,n.tune)
      }
    }
    close(pb)
    
    B_new = B.array[,,which.min(bic.vec)]
    Omega_new = Omega.array[,,which.min(bic.vec)]
    
    # store outputs
    B.list[[iter]] = B_new
    Omega.list[[iter]] = Omega_new
    var.list[[iter]] = list(pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    # check convergence
    best.ind = which(matrix(1:(nl*ng),nl,ng)==which.min(bic.vec), arr.ind=T)
    best.lam = lambda[best.ind[1]]
    best.gam = gamma[best.ind[2]]
    Objval[iter] = Objfunc(X, Z, B_new, Omega_new,
                           best.lam, best.gam,
                           pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    Normdiff[iter] = sqrt(sum((B_new - B_old)^2)/sum(B_new^2)) +
      sqrt(sum((Omega_new - Omega_old)^2)/sum(Omega_new^2))
    
    # Outputs
    cat('Obj_val',Objval[iter],'\n-----\n')
    CONVERGE = (min(Normdiff[iter])<tol)
    
    if (iter == maxit){ # break if max iterations reached
      cat("Max iterations reached.",'\n')
      break;
    }
    if(sum(abs(B_new!=0))==0){ # break if no connections to next layer
      cat("No more connections to latent layer.",'\n')
      break;
    }
  }
  
  if(CONVERGE){
    cat("Converged after",iter,"iterations.\n")
  }
  
  ## return
  list(B.list=B.list, Omega.list=Omega.list,
       var.list=var.list, Obj=Objval[1:iter])
}

hlgm1 = function(X, q, B_init=NULL, Theta_init=NULL, init.option=1,
                lambda=NULL, gamma=NULL, K=3, eta=1e-4,
                refit.B=TRUE, tol=1e-3, grad.max=1e3,
                maxit=20, eps=1e-6, max.epoch=1e2,VERBOSE=TRUE){
  
  #****************************************************#
  # Define and initialize some quantities
  #****************************************************#
  n = nrow(X)
  p = ncol(X)
  
  ## default values of arguments if they are NULL
  if(is.null(lambda)){
    lambda = sqrt(log(p)/n) * seq(2, 1.25, -.25)
  }
  if(is.null(gamma)){
    gamma = sqrt(log(q)/n) * seq(2, 1.25, -.25)
  }
  nl = length(lambda)
  ng = length(gamma)
  
  #****************************************************#
  # Initialization of iterates
  #****************************************************#
  
  ## initialize variational parameters and generate pseudo-data
  set.seed(07232018)
  pi.m = rep(1/K,K)
  mu.m = runif(K,-1,1)
  sig.m = runif(K,0,1)
  pi.s = pi.m
  mu.s = runif(K,0,1)
  sig.s = runif(K,0,1)
  Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
  
  ## Initialize B, if initial B and Omega not supplied
  if(init.option==1){
    ## generate Z from initialized variational params,
    ## then fit l1ML to initialize B and Omega
    cat("Initializing B and Omega\n")
    B.array = array(0,c(p,q,nl*ng))
    Omega.array = array(0, c(p,p,nl*ng))
    
    # # randomly initialize until at least one non-zero B is found
    
    bic.vec = rep(0, nl*ng)
    n.tune = 0
    pb = txtProgressBar(0,nl*ng)
    for(a in 1:nl){
      for(b in 1:ng){
        n.tune = n.tune+1
        init.Obj = l1ML_Main(X, Z, lambda=lambda[a], rho=gamma[b],
                             StabilizeTheta=F, VERBOSE=F)
        # cat("lambda=",lambda[a],"gamma=",gamma[b],
        #     sum(abs(init.Obj$B.est)),"done!\n")
        
        # calculate BIC
        B.array[,,n.tune] = init.Obj$B.est
        Omega.array[,,n.tune] = init.Obj$Theta.est
        bic.vec[n.tune] = with(init.Obj, -log(det(Theta.est)) +
                                 sum(diag(crossprod(X - Z %*% B.est) %*% Theta.est))/n +
                                 log(n)/n * ((sum(Theta.est != 0)-q)/2 +
                                               sum(B.est != 0)))
        setTxtProgressBar(pb,n.tune)
      }
    }
    close(pb)
    
    B_init = B.array[,,which.min(bic.vec)]
    # random initialization if all coefs are 0 in B matrix
    if(sum(abs(B_init))==0){
      for(i in 1:q){
        for(j in 1:q){
          B_init[i,p] = rbinom(1,1,1/q)*sample(c(-1,1),1)*runif(1,0.5,1);
        }
      }
    }
    Omega_init = Omega.array[,,which.min(bic.vec) ]
  }
  
  #****************************************************#
  # Alternating algorithm
  #****************************************************#
  
  # initialize
  iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
  updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
  B_new = B_init
  Omega_new = Omega_init
  Normdiff = rep(0, maxit)
  Objval = rep(0, maxit)
  B.list = list()
  Omega.list = list()
  var.list = list()
  
  ## start with the alternating procedure
  cat('-----\n')
  while(!CONVERGE){
    iter = iter + 1
    B_old = B_new
    Omega_old = Omega_new
    Omega.skeleton = which(Omega_old!=0, arr.ind=T)
    cat("Iteration ", iter, ":\n")
    
    ## Update variational parameters using SGD
    cat("Updating variational parameters\n")
    epoch = 0
    CONVERGE.var=FALSE
    while(!CONVERGE.var){
      
      epoch = epoch+1
      eta.epoch = eta*exp(-epoch+1)
      perm = gtools::permute(1:n) # randomly order samples
      pi.m.old = pi.m
      pi.s.old = pi.s
      mu.m.old = mu.m
      mu.s.old = mu.s
      sig.m.old = sig.m
      sig.s.old = sig.s
      
      ## terms corresp to squared error loss
      for(i in perm){ # SGD: go through all samples randomly
        xi = X[i,]
        
        grad.pi.m = 0
        grad.mu.m = 0
        grad.sig.m = 0
        grad.pi.s = 0
        grad.mu.s = 0
        grad.sig.s = 0
        mult2 = exp(-2*(mu.s + sig.s^2))
        
        for(j0 in 1:nrow(Omega.skeleton)){
          
          # calculate quantities
          j = Omega.skeleton[j0,1]
          j1 = Omega.skeleton[j0,2]
          xj.1.bj1 = xi[j]*sum(B_old[j1,])
          bj.t.bj1 = sum(B_old[j,]*B_old[j1,])
          mult1 = -(2*Omega_old[j,j1]*xj.1.bj1)/n
          # mult2 = exp(clip(2*(mu.s.old + sig.s.old^2), ceil=3))*bj.t.bj1
          # mult2 = exp(2*(mu.s.old + sig.s.old^2))*(bj.t.bj1 +
          #                                            q/2/nrow(Omega.skeleton))
          # calculatre gradients
          grad.pi.m = grad.pi.m + mult1*mu.m + (mu.m^2+sig.m^2)*bj.t.bj1
          grad.mu.m = grad.mu.m + pi.m*(mult1 + 2*mu.m*bj.t.bj1)
          grad.sig.m = grad.sig.m + 2*sig.m*pi.m*bj.t.bj1
          grad.pi.s = grad.pi.s + mult2*bj.t.bj1
          grad.mu.s = grad.mu.s - 2*pi.s*mult2*bj.t.bj1
          grad.sig.s = grad.sig.s - 4*pi.s*sig.s*mult2*bj.t.bj1
        }
        
        ## terms corresp to first KL divergence
        mult3 = mult2/2*q/n
        grad1.pi.m = (mu.m^2+sig.m^2)/2*q/n
        grad1.mu.m = pi.m*mu.m*q/n
        grad1.sig.m = pi.m*sig.m*q/n
        grad1.pi.s = mult3 + mu.s*q/n
        grad1.mu.s = -2*pi.s.old*mult3 + 1*q/n
        grad1.sig.s = -4*pi.s.old*sig.s.old*mult3
        
        ## terms corresp. to second KL divergence
        # update parameters of r(.)
        phi.m.sq = clip(sum(pi.m/sig.m^2), ceil=1e3)
        phi.s.sq = clip(sum(pi.s/sig.s^2), ceil=1e3)
        
        # calculate gradients
        sumZi2 = sum(Z[i,]^2); sumZi = sum(Z[i,])
        grad2.pi.m = (sumZi2 - 2*sumZi*mu.m +
                        q*(mu.m^2 + phi.m.sq))/2/n/sig.m^2 + q*log(sig.m)
        grad2.mu.m = pi.m*(q*mu.m - sumZi)/n/sig.m^2
        grad2.sig.m = -pi.m*(sumZi2 - 2*sumZi*mu.m +
                               q*(mu.m^2 + phi.m.sq))/n/sig.m^3 + clip(q/sig.m)
        grad2.pi.s = (sumZi2 - 2*sumZi*mu.s +
                        q*(mu.s^2 + phi.s.sq))/2/n/sig.s^2 + q*log(sig.s)
        grad2.mu.s = pi.s*(q*mu.s - sumZi)/n/sig.s^2
        grad2.sig.s = -pi.s*(sumZi2 - 2*sumZi*mu.s +
                               q*(mu.s^2 + phi.s.sq))/n/sig.s^3 + clip(q/sig.s)
        
        # update
        pi.m = pi.m - eta.epoch*clip(grad.pi.m+grad1.pi.m+grad2.pi.m)
        mu.m = mu.m - eta.epoch*clip(grad.mu.m+grad2.mu.m+grad2.mu.m)
        sig.m = sig.m - eta.epoch*clip(grad.sig.m+grad1.sig.m+grad2.sig.m)
        pi.s = pi.s - eta.epoch*clip(grad.pi.s+grad1.pi.s+grad2.pi.s)
        mu.s = mu.s - eta.epoch*clip(grad.mu.s+grad1.mu.s+grad2.mu.s)
        sig.s = sig.s - eta.epoch*clip(grad.sig.s+grad1.sig.s+grad2.sig.s)

        # normalize parameters
        pi.m = abs(pi.m)/sum(abs(pi.m))
        pi.s = abs(pi.s)/sum(abs(pi.s))
        sig.m = abs(sig.m)
        sig.s = abs(sig.s)
      }
      # cat('\t',sig.s,'\n----------\n')

      # check convergence
      var.diff = sqrt(sum((pi.m.old-pi.m)^2)/sum(pi.m^2)) +
        sqrt(sum((pi.s.old-pi.s)^2)/sum(pi.s^2)) +
        sqrt(sum((mu.m.old-mu.m)^2)/sum(mu.m^2)) +
        sqrt(sum((mu.s.old-mu.s)^2)/sum(mu.s^2)) +
        sqrt(sum((sig.m.old-sig.m)^2)/sum(sig.m^2)) +
        sqrt(sum((sig.s.old-sig.s)^2)/sum(sig.s^2))
      # cat("Norm_diff =",round(var.diff,4),sig.s,'\n-----\n')
      CONVERGE.var = (var.diff<tol)
      if(epoch==max.epoch){
        cat("Max iterations reached.",'\n')
        break;
      }
    }
    
    ## generate pseudo-data Z from variational distributions
    Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    ## Update generative parameters B, Omega
    cat("Updating B and Omega\n")
    B.array = array(0,c(p,q,nl*ng))
    Omega.array = array(0, c(q,q,nl*ng))
    bic.vec = rep(0, nl*ng)
    n.tune = 0
    pb = txtProgressBar(0, nl*ng)
    for(a in 1:nl){
      for(b in 1:ng){
        n.tune = n.tune+1
        ML.Obj = l1ML_Main(X, Z, lambda=lambda[a], rho=gamma[b],
                           StabilizeTheta=F, VERBOSE=F)
        # cat("lambda=",lambda[a],"gamma=",gamma[b],
        #     sum(abs(ML.Obj$B.est)),"done!\n")
        
        # calculate BIC
        B.array[,,n.tune] = ML.Obj$B.est
        Omega.array[,,n.tune] = ML.Obj$Theta.est
        bic.vec[n.tune] = with(ML.Obj, -log(det(Theta.est)) +
                                 sum(diag(crossprod(X - Z %*% B.est) %*% Theta.est))/n +
                                 log(n)/n * ((sum(Theta.est != 0)-q)/2 +
                                               sum(B.est != 0)))
        setTxtProgressBar(pb,n.tune)
      }
    }
    close(pb)
    
    B_new = B.array[,,which.min(bic.vec)]
    Omega_new = Omega.array[,,which.min(bic.vec)]
    
    # store outputs
    B.list[[iter]] = B_new
    Omega.list[[iter]] = Omega_new
    var.list[[iter]] = list(pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    # check convergence
    best.ind = which(matrix(1:(nl*ng),nl,ng)==which.min(bic.vec), arr.ind=T)
    best.lam = lambda[best.ind[1]]
    best.gam = gamma[best.ind[2]]
    Objval[iter] = Objfunc(X, Z, B_new, Omega_new,
                           best.lam, best.gam,
                           pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    Normdiff[iter] = sqrt(sum((B_new - B_old)^2)/sum(B_new^2)) +
      sqrt(sum((Omega_new - Omega_old)^2)/sum(Omega_new^2))
    
    # Outputs
    cat('Obj_val',Objval[iter],'\n-----\n')
    CONVERGE = (min(Normdiff[iter])<tol)
    
    if (iter == maxit){ # break if max iterations reached
      cat("Max iterations reached.",'\n')
      break;
    }
    if(sum(abs(B_new!=0))==0){ # break if no connections to next layer
      cat("No more connections to latent layer.",'\n')
      break;
    }
  }
  
  if(CONVERGE){
    cat("Converged after",iter,"iterations.\n")
  }
  
  ## return
  list(B.list=B.list, Omega.list=Omega.list,
       var.list=var.list, Obj=Objval[1:iter])
}

# HLGM with JMMLE
hlgm2 = function(X, q, B_init=NULL, Theta_init=NULL, init.option=1,
                 lambda=NULL, gamma=NULL, K=3, eta=1e-4,
                 refit.B=TRUE, momentum=1, tol=1e-3, grad.max=1e3,
                 maxit=20, eps=1e-6, max.epoch=1e2,VERBOSE=TRUE){
  
  #****************************************************#
  # Define and initialize some quantities
  #****************************************************#
  n = nrow(X)
  p = ncol(X)
  
  ## default values of arguments if they are NULL
  if(is.null(lambda)){
    lambda = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
  }
  if(is.null(gamma)){
    gamma = sqrt(log(q)/n) * seq(1.1, 0.1, -0.2)
  }
  nl = length(lambda)
  ng = length(gamma)
  
  #****************************************************#
  # Initialization of iterates
  #****************************************************#
  
  ## initialize variational parameters and generate pseudo-data
  set.seed(07232018)
  pi.m = rep(1/K,K)
  mu.m = runif(K,-1,1)
  sig.m = runif(K,0,1)
  pi.s = pi.m
  mu.s = runif(K,0,1)
  sig.s = runif(K,0,.1)
  Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
  
  ## Initialize B, if initial B and Omega not supplied
  if(init.option==1){
    ## generate Z from initialized variational params,
    ## then fit l1ML to initialize B and Omega
    cat("Initializing B and Omega\n")

    ## tune JMLE model
    pb = txtProgressBar(0, nl)
    
    ## get all models
    B.group.array = array(matrix(1:(p*q), nrow=p, byrow=T), dim=c(q,p,1))
    Theta.groups=rep(list(matrix(1,p-1,1)),p)
    loopfun1 = function(m){
      Obj = jmmle.1step(Y.list=list(X), X.list=list(Z),
                        B.group.array=B.group.array,
                        Theta.groups=Theta.groups,
                        lambda = lambda[m],
                        gamma = gamma,
                        init.option=1, tol=1e-3, VERBOSE=F)
      setTxtProgressBar(pb,m)
      Obj
    }
    # system.time(
      model.list <- lapply(1:nl, loopfun1)
    # )
    close(pb)
    
    ## calculate HBIC
    hbic.vec = rep(NA, nl)
    for(m in 1:nl){
      jmle.model = model.list[[m]]
      
      if(class(jmle.model)=="list"){ ## if no error in training the model
        Theta.k = jmle.model$Theta_refit$Theta[[1]]
        for(j in 1:q)
        {
          Theta.k[j,j] = 0
        }
        SSE = sum(diag(crossprod((X - Z %*% jmle.model$B.refit[,,1]) %*%
                                   (diag(1,p) - Theta.k))))/n
        hbic.pen = log(log(n))*log(p*(p-1)/2)/n * sum(Theta.k != 0)/2 +
          log(log(n))*log(p*q)/n * sum(jmle.model$B.refit[,,1] != 0)
        hbic.vec[m] = SSE + hbic.pen 
      }
    }
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    B_init = jmmle.model$B.refit[,,1]
    # random initialization if all coefs are 0 in B matrix
    if(sum(abs(B_init))==0){
      for(i in 1:q){
        for(j in 1:p){
          B_init[i,j] = rbinom(1,1,1/q)*sample(c(-1,1),1)*runif(1,0.5,1);
        }
      }
    }
    Omega_init = jmmle.model$Theta_refit$Omega[[1]]
    Theta_init = jmmle.model$Theta_refit$Theta[[1]]
  }
  
  #****************************************************#
  # Alternating algorithm
  #****************************************************#
  
  # initialize
  iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
  updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
  B_new = B_init
  Omega_new = Omega_init
  Theta_new = Theta_init
  B.array = array(0, dim=c(q,p,1))
  Theta.array = array(0, dim=c(p,p,1))
  
  Normdiff = rep(0, maxit)
  Objval = rep(0, maxit)
  B.list = list()
  Omega.list = list()
  var.list = list()
  v = rep(list(0),6) # momentum
  
  ## start with the alternating procedure
  cat('-----\n')
  while(!CONVERGE){
    iter = iter + 1
    B_old = B_new
    Omega_old = Omega_new
    Theta_old = Theta_new
    B.array[,,1] = B_old
    Theta.array[,,1] = Theta_old
    Omega.skeleton = which(Omega_old!=0, arr.ind=T)
    cat("Iteration ", iter, ":\n")
    
    ## E step: Update variational parameters using SGD
    cat("Updating variational parameters\n")
    epoch = 0
    CONVERGE.var=FALSE
    while(!CONVERGE.var){
      
      epoch = epoch+1
      eta.epoch = eta*exp(-epoch+1)
      perm = gtools::permute(1:n) # randomly order samples
      pi.m.old = pi.m
      pi.s.old = pi.s
      mu.m.old = mu.m
      mu.s.old = mu.s
      sig.m.old = sig.m
      sig.s.old = sig.s
      
      ## terms corresp to squared error loss
      for(i in perm){ # SGD: go through all samples randomly
        xi = X[i,]
        
        grad.pi.m = 0
        grad.mu.m = 0
        grad.sig.m = 0
        grad.pi.s = 0
        grad.mu.s = 0
        grad.sig.s = 0
        mult2 = exp(-2*(mu.s + sig.s^2))
        
        for(j0 in 1:nrow(Omega.skeleton)){
          
          # calculate quantities
          j = Omega.skeleton[j0,1]
          j1 = Omega.skeleton[j0,2]
          xj.1.bj1 = xi[j]*sum(B_old[,j1])
          bj.t.bj1 = sum(B_old[,j]*B_old[,j1])
          mult1 = -(2*Omega_old[j,j1]*xj.1.bj1)/n
          # mult2 = exp(clip(2*(mu.s.old + sig.s.old^2), ceil=3))*bj.t.bj1
          # mult2 = exp(2*(mu.s.old + sig.s.old^2))*(bj.t.bj1 +
          #                                            q/2/nrow(Omega.skeleton))
          # calculatre gradients
          grad.pi.m = grad.pi.m + mult1*mu.m + (mu.m^2+sig.m^2)*bj.t.bj1
          grad.mu.m = grad.mu.m + pi.m*(mult1 + 2*mu.m*bj.t.bj1)
          grad.sig.m = grad.sig.m + 2*sig.m*pi.m*bj.t.bj1
          grad.pi.s = grad.pi.s + mult2*bj.t.bj1
          grad.mu.s = grad.mu.s - 2*pi.s*mult2*bj.t.bj1
          grad.sig.s = grad.sig.s - 4*pi.s*sig.s*mult2*bj.t.bj1
        }
        
        ## terms corresp to first KL divergence
        mult3 = mult2/2*q/n
        grad1.pi.m = (mu.m^2+sig.m^2)/2*q/n
        grad1.mu.m = pi.m*mu.m*q/n
        grad1.sig.m = pi.m*sig.m*q/n
        grad1.pi.s = mult3 + mu.s*q/n
        grad1.mu.s = -2*pi.s.old*mult3 + 1*q/n
        grad1.sig.s = -4*pi.s.old*sig.s.old*mult3
        
        ## terms corresp. to second KL divergence
        # update parameters of r(.)
        phi.m.sq = clip(sum(pi.m/sig.m^2), ceil=1e3)
        phi.s.sq = clip(sum(pi.s/sig.s^2), ceil=1e3)
        
        # calculate gradients
        sumZi2 = sum(Z[i,]^2); sumZi = sum(Z[i,])
        grad2.pi.m = (sumZi2 - 2*sumZi*mu.m +
                        q*(mu.m^2 + phi.m.sq))/2/n/sig.m^2 + q/n*log(sig.m)
        grad2.mu.m = pi.m*(q*mu.m - sumZi)/n/sig.m^2
        grad2.sig.m = -pi.m*(sumZi2 - 2*sumZi*mu.m +
                               q*(mu.m^2 + phi.m.sq))/n/sig.m^3 + clip(q/n/sig.m)
        grad2.pi.s = (sumZi2 - 2*sumZi*mu.s +
                        q*(mu.s^2 + phi.s.sq))/2/n/sig.s^2 + q/n*log(sig.s)
        grad2.mu.s = pi.s*(q*mu.s - sumZi)/n/sig.s^2
        grad2.sig.s = -pi.s*(sumZi2 - 2*sumZi*mu.s +
                               q*(mu.s^2 + phi.s.sq))/n/sig.s^3 + clip(q/n/sig.s)
        
        # update
        if(momentum==1){
          # SGD with momentum
          gam.v = .5#ifelse(iter<5, 0.5, 0.9)
          v = lapply(v, function(x) gam.v*x)
          v[[1]] = v[[1]] + eta.epoch*clip(grad.pi.m+grad1.pi.m+grad2.pi.m)
          v[[2]] = v[[2]] + eta.epoch*clip(grad.mu.m+grad2.mu.m+grad2.mu.m)
          v[[3]] = v[[3]] + eta.epoch*clip(grad.sig.m+grad1.sig.m+grad2.sig.m)
          v[[4]] = v[[4]] + eta.epoch*clip(grad.pi.s+grad1.pi.s+grad2.pi.s)
          v[[5]] = v[[5]] + eta.epoch*clip(grad.mu.s+grad1.mu.s+grad2.mu.s)
          v[[6]] = v[[6]] + eta.epoch*clip(grad.sig.s+grad1.sig.s+grad2.sig.s)
          
          pi.m = pi.m - v[[1]]
          mu.m = mu.m - v[[2]]
          sig.m = sig.m - v[[3]]
          pi.s = pi.s - v[[4]]
          mu.s = mu.s - v[[5]]
          sig.s = sig.s - v[[6]]
        } else{
          # vanilla SGD
          pi.m = pi.m - eta.epoch*clip(grad.pi.m+grad1.pi.m+grad2.pi.m)
          mu.m = mu.m - eta.epoch*clip(grad.mu.m+grad2.mu.m+grad2.mu.m)
          sig.m = sig.m - eta.epoch*clip(grad.sig.m+grad1.sig.m+grad2.sig.m)
          pi.s = pi.s - eta.epoch*clip(grad.pi.s+grad1.pi.s+grad2.pi.s)
          mu.s = mu.s - eta.epoch*clip(grad.mu.s+grad1.mu.s+grad2.mu.s)
          sig.s = sig.s - eta.epoch*clip(grad.sig.s+grad1.sig.s+grad2.sig.s)
        }
        # normalize parameters
        pi.m = abs(pi.m)/sum(abs(pi.m))
        pi.s = abs(pi.s)/sum(abs(pi.s))
        sig.m = abs(sig.m)
        sig.s = abs(sig.s)
      }
      # cat('\t',sig.s,'\n----------\n')
      
      # check convergence
      var.diff = sqrt(sum((pi.m.old-pi.m)^2)/sum(pi.m^2)) +
        sqrt(sum((pi.s.old-pi.s)^2)/sum(pi.s^2)) +
        sqrt(sum((mu.m.old-mu.m)^2)/sum(mu.m^2)) +
        sqrt(sum((mu.s.old-mu.s)^2)/sum(mu.s^2)) +
        sqrt(sum((sig.m.old-sig.m)^2)/sum(sig.m^2)) +
        sqrt(sum((sig.s.old-sig.s)^2)/sum(sig.s^2))
      # cat("Norm_diff =",round(var.diff,4),sig.s,'\n-----\n')
      CONVERGE.var = (var.diff<tol)
      if(epoch==max.epoch){
        cat("Max iterations reached.",'\n')
        break;
      }
    }
    
    # ## M step: Update generative parameters B, Omega
    # m1 = sum(pi.m*mu.m)
    # m2 = sum(pi.m*(mu.m^2+sig.m^2))
    # m3 = sum(pi.s*exp(-2*(mu.s + sig.s^2)))
    # X1 = X*m1/(m2+m3)
    # Z1 = matrix(1, n, q)*sqrt(m2+m3)
    
    # Generate pseudo-data
    Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    cat("Updating B and Omega\n")
    ## tune JMLE model
    pb = txtProgressBar(0, nl)
    
    ## get all models
    loopfun1 = function(m){
      Obj = jmmle.1step(Y.list=list(X), X.list=list(Z),
                        B.group.array=B.group.array,
                        Theta.groups=Theta.groups,
                        B_init.array=B.array,
                        Theta_init.array=Theta.array,
                        init.gamma=jmle.model$best.gamma,
                        lambda=lambda[m], gamma=gamma,
                        init.option=1, tol=1e-3, VERBOSE=F)
      setTxtProgressBar(pb,m)
      Obj
    }
    # system.time(
      model.list <- lapply(1:nl, loopfun1)
    # )
    close(pb)
    
    ## calculate HBIC
    hbic.vec = rep(NA, nl)
    for(m in 1:nl){
      jmle.model = model.list[[m]]
      
      if(class(jmle.model)=="list"){ ## if no error in training the model
        Theta.k = jmle.model$Theta_refit$Theta[[1]]
        for(j in 1:q)
        {
          Theta.k[j,j] = 0
        }
        SSE = sum(diag(crossprod((X - Z %*% jmle.model$B.refit[,,1]) %*%
                                   (diag(1,p) - Theta.k))))/n
        hbic.pen = log(log(n))*log(p*(p-1)/2)/n * sum(Theta.k != 0)/2 +
          log(log(n))*log(p*q)/n * sum(jmle.model$B.refit[,,1] != 0)
        hbic.vec[m] = SSE + hbic.pen 
      }
    }
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    B_new = jmmle.model$B.refit[,,1]
    Omega_new = jmmle.model$Theta_refit$Omega[[1]]
    Theta_new = jmmle.model$Theta_refit$Theta[[1]]
    
    # store outputs
    B.list[[iter]] = B_new
    Omega.list[[iter]] = Omega_new
    var.list[[iter]] = list(pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    # check convergence
    best.lam = lambda[which.min(hbic.vec)]
    best.gam = jmmle.model$best.gamma
    Objval[iter] = Objfunc(X, Z, B_new, Omega_new,
                           best.lam, best.gam,
                           pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    Normdiff[iter] = sqrt(sum((B_new - B_old)^2)/sum(B_new^2)) +
      sqrt(sum((Omega_new - Omega_old)^2)/sum(Omega_new^2))
    
    # Outputs
    cat('Obj_val',Objval[iter])
    CONVERGE = (Normdiff[iter]<tol)
    
    # change momentum variable if obj function didn't decrease
    # if(iter==1){
    #   cat("Momentum is",momentum)
    # } else{
    #   flag.obj = (Objval[iter]<Objval[iter-1])
    #   if(!flag.obj){
    #     momentum = 1-momentum
    #     cat("Obj function didn't go down. Momentum changed to",
    #         momentum,".")
    #   }
    # }
    cat("\n-----\n")
    
    if (iter == maxit){ # break if max iterations reached
      cat("Max iterations reached.",'\n')
      break;
    }
    if(sum(abs(B_new!=0))==0){ # break if no connections to next layer
      cat("No more connections to latent layer.",'\n')
      break;
    }
  }
  
  if(CONVERGE){
    cat("Converged after",iter,"iterations.\n")
  }
  
  ## return
  list(B.list=B.list, Omega.list=Omega.list,
       var.list=var.list, Obj=Objval[1:iter])
}

# HLGM from second level onwards: estimate only diagonals of residual matrix
hlgm.diag = function(X, q, B_init=NULL, init.option=1,
                 lambda=NULL, K=3, eta=1e-4,
                 refit.B=TRUE, momentum=1, tol=1e-3, grad.max=1e3,
                 maxit=20, eps=1e-6, max.epoch=1e2,VERBOSE=TRUE){
  
  #****************************************************#
  # Define and initialize some quantities
  #****************************************************#
  n = nrow(X)
  p = ncol(X)
  
  ## default values of arguments if they are NULL
  if(is.null(lambda)){
    lambda = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
  }
  nl = length(lambda)

  #****************************************************#
  # Initialization of iterates
  #****************************************************#
  
  ## initialize variational parameters and generate pseudo-data
  set.seed(07232018)
  pi.m = rep(1/K,K)
  mu.m = runif(K,-1,1)
  sig.m = runif(K,0,1)
  pi.s = pi.m
  mu.s = runif(K,0,1)
  sig.s = runif(K,0,.1)
  Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
  
  ## Initialize B, if initial B and Omega not supplied
  if(init.option==1){
    ## generate Z from initialized variational params,
    ## then fit l1ML to initialize B and Omega
    cat("Initializing B and Omega\n")
    
    ## tune JMLE model
    pb = txtProgressBar(0, nl)
    
    ## get all models
    loopfun1 = function(m){
      Obj = l1LS_Main(X, Z, skeleton.hat = array(1, c(q,p)),
                      lambda=lambda[m])
      setTxtProgressBar(pb,m)
      Obj
    }
    # system.time(
    model.list <- lapply(1:nl, loopfun1)
    # )
    close(pb)
    
    ## calculate HBIC
    hbic.vec = rep(NA, nl)
    for(m in 1:nl){
      jmle.model = model.list[[m]]
      
      if(class(jmle.model)=="list"){ ## if no error in training the model
        B.hat = jmle.model$B0
        Omega.diag = 1/diag(var(X - Z %*% B.hat))
        for(j in 1:q)
        SSE = sum(diag(crossprod((X - Z %*% B.hat) %*% Omega.diag)))/n
        hbic.pen = log(log(n))*log(p*q)/n * sum(B.hat != 0)
        hbic.vec[m] = SSE + hbic.pen 
      }
    }
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    B_init = jmmle.model$B0
    # random initialization if all coefs are 0 in B matrix
    if(sum(abs(B_init))==0){
      for(i in 1:q){
        for(j in 1:p){
          B_init[i,j] = rbinom(1,1,1/q)*sample(c(-1,1),1)*runif(1,0.5,1);
        }
      }
    }
    Omega_init = 1/diag(var(X - Z %*% B_init))
  }
  
  #****************************************************#
  # Alternating algorithm
  #****************************************************#
  
  # initialize
  iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
  updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
  B_new = B_init
  Omega_new = Omega_init

  Normdiff = rep(0, maxit)
  Objval = rep(0, maxit)
  B.list = list()
  Omega.list = list()
  var.list = list()
  v = rep(list(0),6) # momentum
  
  ## start with the alternating procedure
  cat('-----\n')
  while(!CONVERGE){
    iter = iter + 1
    B_old = B_new
    Omega_old = Omega_new
    cat("Iteration ", iter, ":\n")
    
    ## E step: Update variational parameters using SGD
    cat("Updating variational parameters\n")
    epoch = 0
    CONVERGE.var=FALSE
    while(!CONVERGE.var){
      
      epoch = epoch+1
      eta.epoch = eta*exp(-epoch+1)
      perm = gtools::permute(1:n) # randomly order samples
      pi.m.old = pi.m
      pi.s.old = pi.s
      mu.m.old = mu.m
      mu.s.old = mu.s
      sig.m.old = sig.m
      sig.s.old = sig.s
      
      ## terms corresp to squared error loss
      for(i in perm){ # SGD: go through all samples randomly
        xi = X[i,]
        
        grad.pi.m = 0
        grad.mu.m = 0
        grad.sig.m = 0
        grad.pi.s = 0
        grad.mu.s = 0
        grad.sig.s = 0
        mult2 = exp(-2*(mu.s + sig.s^2))
        
        for(j in 1:q){
          
          # calculate quantities
          xj.1.bj1 = xi[j]*sum(B_old[,j])
          bj.t.bj1 = sum(B_old[,j]*B_old[j,j])
          mult1 = -(2*Omega_old[j]*xj.1.bj1)/n
          # mult2 = exp(clip(2*(mu.s.old + sig.s.old^2), ceil=3))*bj.t.bj1
          # mult2 = exp(2*(mu.s.old + sig.s.old^2))*(bj.t.bj1 +
          #                                            q/2/nrow(Omega.skeleton))
          # calculatre gradients
          grad.pi.m = grad.pi.m + mult1*mu.m + (mu.m^2+sig.m^2)*bj.t.bj1
          grad.mu.m = grad.mu.m + pi.m*(mult1 + 2*mu.m*bj.t.bj1)
          grad.sig.m = grad.sig.m + 2*sig.m*pi.m*bj.t.bj1
          grad.pi.s = grad.pi.s + mult2*bj.t.bj1
          grad.mu.s = grad.mu.s - 2*pi.s*mult2*bj.t.bj1
          grad.sig.s = grad.sig.s - 4*pi.s*sig.s*mult2*bj.t.bj1
        }
        
        ## terms corresp to first KL divergence
        mult3 = mult2/2*q/n
        grad1.pi.m = (mu.m^2+sig.m^2)/2*q/n
        grad1.mu.m = pi.m*mu.m*q/n
        grad1.sig.m = pi.m*sig.m*q/n
        grad1.pi.s = mult3 + mu.s*q/n
        grad1.mu.s = -2*pi.s.old*mult3 + 1*q/n
        grad1.sig.s = -4*pi.s.old*sig.s.old*mult3
        
        ## terms corresp. to second KL divergence
        # update parameters of r(.)
        phi.m.sq = clip(sum(pi.m/sig.m^2), ceil=1e3)
        phi.s.sq = clip(sum(pi.s/sig.s^2), ceil=1e3)
        
        # calculate gradients
        sumZi2 = sum(Z[i,]^2); sumZi = sum(Z[i,])
        grad2.pi.m = (sumZi2 - 2*sumZi*mu.m +
                        q*(mu.m^2 + phi.m.sq))/2/n/sig.m^2 + q/n*log(sig.m)
        grad2.mu.m = pi.m*(q*mu.m - sumZi)/n/sig.m^2
        grad2.sig.m = -pi.m*((sumZi2 - 2*sumZi*mu.m +
                                q*(mu.m^2 + phi.m.sq))/n/sig.m^3 + clip(q/n/sig.m))
        grad2.pi.s = (sumZi2 - 2*sumZi*mu.s +
                        q*(mu.s^2 + phi.s.sq))/2/n/sig.s^2 + q/n*log(sig.s)
        grad2.mu.s = pi.s*(q*mu.s - sumZi)/n/sig.s^2
        grad2.sig.s = -pi.s*((sumZi2 - 2*sumZi*mu.s +
                                q*(mu.s^2 + phi.s.sq))/n/sig.s^3 + clip(q/n/sig.s))
        
        # update
        if(momentum==1){
          # SGD with momentum
          gam.v = .5#ifelse(iter<5, 0.5, 0.9)
          v = lapply(v, function(x) gam.v*x)
          v[[1]] = v[[1]] + eta.epoch*clip(grad.pi.m+grad1.pi.m+grad2.pi.m)
          v[[2]] = v[[2]] + eta.epoch*clip(grad.mu.m+grad2.mu.m+grad2.mu.m)
          v[[3]] = v[[3]] + eta.epoch*clip(grad.sig.m+grad1.sig.m+grad2.sig.m)
          v[[4]] = v[[4]] + eta.epoch*clip(grad.pi.s+grad1.pi.s+grad2.pi.s)
          v[[5]] = v[[5]] + eta.epoch*clip(grad.mu.s+grad1.mu.s+grad2.mu.s)
          v[[6]] = v[[6]] + eta.epoch*clip(grad.sig.s+grad1.sig.s+grad2.sig.s)
          
          pi.m = pi.m - v[[1]]
          mu.m = mu.m - v[[2]]
          sig.m = sig.m - v[[3]]
          pi.s = pi.s - v[[4]]
          mu.s = mu.s - v[[5]]
          sig.s = sig.s - v[[6]]
        } else{
          # vanilla SGD
          pi.m = pi.m - eta.epoch*clip(grad.pi.m+grad1.pi.m+grad2.pi.m)
          mu.m = mu.m - eta.epoch*clip(grad.mu.m+grad2.mu.m+grad2.mu.m)
          sig.m = sig.m - eta.epoch*clip(grad.sig.m+grad1.sig.m+grad2.sig.m)
          pi.s = pi.s - eta.epoch*clip(grad.pi.s+grad1.pi.s+grad2.pi.s)
          mu.s = mu.s - eta.epoch*clip(grad.mu.s+grad1.mu.s+grad2.mu.s)
          sig.s = sig.s - eta.epoch*clip(grad.sig.s+grad1.sig.s+grad2.sig.s)
        }
        # normalize parameters
        pi.m = abs(pi.m)/sum(abs(pi.m))
        pi.s = abs(pi.s)/sum(abs(pi.s))
        sig.m = abs(sig.m)
        sig.s = abs(sig.s)
      }
      # cat('\t',sig.s,'\n----------\n')
      
      # check convergence
      var.diff = sqrt(sum((pi.m.old-pi.m)^2)/sum(pi.m^2)) +
        sqrt(sum((pi.s.old-pi.s)^2)/sum(pi.s^2)) +
        sqrt(sum((mu.m.old-mu.m)^2)/sum(mu.m^2)) +
        sqrt(sum((mu.s.old-mu.s)^2)/sum(mu.s^2)) +
        sqrt(sum((sig.m.old-sig.m)^2)/sum(sig.m^2)) +
        sqrt(sum((sig.s.old-sig.s)^2)/sum(sig.s^2))
      # cat("Norm_diff =",round(var.diff,4),sig.s,'\n-----\n')
      CONVERGE.var = (var.diff<tol)
      if(epoch==max.epoch){
        cat("Max iterations reached.",'\n')
        break;
      }
    }
    
    # ## M step: Update generative parameters B, Omega
    # m1 = sum(pi.m*mu.m)
    # m2 = sum(pi.m*(mu.m^2+sig.m^2))
    # m3 = sum(pi.s*exp(-2*(mu.s + sig.s^2)))
    # X1 = X*m1/(m2+m3)
    # Z1 = matrix(1, n, q)*sqrt(m2+m3)
    
    # Generate pseudo-data
    Z = GenerateZ(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    cat("Updating B and Omega\n")
    ## tune JMLE model
    pb = txtProgressBar(0, nl)
    
    ## get all models
    loopfun1 = function(m){
      Obj = l1LS_Main(X, Z, skeleton.hat = array(1, c(q,p)),
                      lambda=lambda[m])
      setTxtProgressBar(pb,m)
      Obj
    }
    # system.time(
    model.list <- lapply(1:nl, loopfun1)
    # )
    close(pb)
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    B_new = jmmle.model$B0
    Omega_new = 1/diag(var(X - Z %*% B_new))
    
    # store outputs
    B.list[[iter]] = B_new
    Omega.list[[iter]] = Omega_new
    var.list[[iter]] = list(pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    
    # check convergence
    best.lam = lambda[which.min(hbic.vec)]
    best.gam = jmmle.model$best.gamma
    Objval[iter] = Objfunc(X, Z, B_new, diag(Omega_new),
                           best.lam, gamma=0,
                           pi.m, mu.m, sig.m, pi.s, mu.s, sig.s)
    Normdiff[iter] = sqrt(sum((B_new - B_old)^2)/sum(B_new^2)) +
      sqrt(sum((Omega_new - Omega_old)^2)/sum(Omega_new^2))
    
    # Outputs
    cat('Obj_val',Objval[iter],".")
    CONVERGE = (Normdiff[iter]<tol)
    
    # change momentum variable if obj function didn't decrease
    if(iter==1){
      cat("Momentum is",momentum)
    } else{
      flag.obj = (Objval[iter]<Objval[iter-1])
      if(!flag.obj){
        momentum = 1-momentum
        cat("Obj function didn't go down. Momentum changed to",
            momentum,".")
      }
    }
    cat("\n-----\n")
    
    if (iter == maxit){ # break if max iterations reached
      cat("Max iterations reached.",'\n')
      break;
    }
    if(sum(abs(B_new!=0))==0){ # break if no connections to next layer
      cat("No more connections to latent layer.",'\n')
      break;
    }
  }
  
  if(CONVERGE){
    cat("Converged after",iter,"iterations.\n")
  }
  
  ## return
  list(B.list=B.list, Omega.list=Omega.list,
       var.list=var.list, Obj=Objval[1:iter])
}

## misc
# if (iter == 1){
#   Norm_diff = Normfunc[1]
# }
# else{
#   Norm_diff = Normfunc[iter] - Normfunc[iter-1]
#   Obj_diff = (Objfunc[iter] - Objfunc[iter-1])/Objfunc[iter-1]
# }
# Obj_diff = Objval[iter+1]/Objval[iter] - 1

# convergence criterion value
# nd = abs(Normdiff[iter])
# if(iter>1){
#   nd = c(nd, abs(rev(diff(Normdiff))[1]))
# }
# if(iter>2){
#   nd = c(nd, abs(rev(diff(Normdiff,2))[1]))
# }

#****************************************************#
# EOF
#****************************************************#