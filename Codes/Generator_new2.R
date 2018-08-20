# ThetaMat() generates the precision matrix
ThetaMat = function(q,type="random",CN=NULL,
                    Theta1=NULL,Theta2=NULL, SNR=1, prob=NULL){
	# by default, generates a random graph with sparsity level being 5/q;
	# if conditional number is not specified, the default will be the dimension of the matrix (as in Cai)
	# alternatively, type could be "band", with the off-diagonals provided
	
	Theta = array(0,c(q,q))
	if(is.null(prob)){
	  prob = 5/q
	}
	
	if (type=="random"){
		diag(Theta) = 0;
		if (is.null(CN))
			CN = q;
		# for off-diagonals:
		for (i in 1:(q-1)){
			for (j in (i+1):q){
				Theta[j,i] = SNR*rbinom(1,1,prob)*sample(c(-1,1),1)*runif(1,0.5,1);
				Theta[i,j] = Theta[j,i];
			}
		}
		# now bump-up the diagonal to get the desired condition number:
		Theta = Theta + diag(q)*(0.001+abs(min(eigen(Theta)$values))); 
		egval = eigen(Theta)$values; 
		CN_Theta = max(egval)/min(egval)
		
		while(CN_Theta>CN){
			Theta = Theta + 0.01*diag(q);
			egval = eigen(Theta)$values;
			CN_Theta = max(egval)/min(egval);
		}
		return(list(Theta=Theta,CN_Theta=CN_Theta))
	}
	if (type=="band"){
		diag(Theta) = 1;
		if (is.null(Theta2)){ # tridiagonal, inverse matrix of AR(1)
			if (abs(Theta1)>=1){
				stop("Wrong Input Value for Theta1, Theta needs to be diagonally dominant.\n")
			}
			for (i in 1:(q-1)){
				for (j in (i+1):q){
					Theta[i,j] = ifelse(j-i==1,Theta1,0);
					Theta[j,i] = Theta[i,j];
				}
			}	
		}
		else{	# two bands 
			if ((abs(Theta1)+abs(Theta2))>1)
				stop("Wrong Input Value for Off-Diagonals, Theta needs to be diagonally dominant.\n")
			for (i in 1:(q-1)){
				for (j in (i+1):q){
					Theta[i,j] = ifelse(j-i==1,sample(c(-1,1),1)*Theta1,ifelse(j-i==2,sample(c(-1,1),1)*Theta2,0));
					Theta[j,i] = Theta[i,j];
				}
			}
		}
		egval = eigen(Theta)$values;
		CN_Theta = max(egval)/min(egval);
		return(list(Theta=Theta,CN_Theta=CN_Theta))
	}
}

# ThetaXYMat() generates the regression matrix
ThetaXYMat = function(p,q,sparsity=NULL,SNR=NULL){
	# p regressors, q regressions
	
	if (is.null(sparsity))
		sparsity = 5/p;
	
	ThetaXY = array(0,c(p,q));
	if (is.null(SNR)){
		for (i in 1:p)
			for (j in 1:q)
				ThetaXY[i,j] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,0.8,1);
	}
	else{
		for (i in 1:p)
			for (j in 1:q)
				ThetaXY[i,j] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,SNR-0.05,SNR+0.05);
	}
	return(ThetaXY)
}

# SigmaMatGiant() generates the big covariance matrix
SigmaMatGiant = function(ThetaX,ThetaXY,ThetaY){
	ThetaYX = t(ThetaXY);
	ThetaGiant = rbind(cbind(ThetaX,ThetaXY),cbind(ThetaYX,ThetaY));
	
	eigmin = min(eigen(ThetaGiant)$values);
	if (eigmin<0)
		ThetaGiant = ThetaGiant + (0.2+abs(eigmin))*diag(nrow(ThetaGiant))
	eig.info = eigen(ThetaGiant);
	CN = max(abs(eig.info$values))/min(abs(eig.info$values))
	SigmaGiant = solve(ThetaGiant)
	return(list(SigmaGiant=SigmaGiant,ThetaCN = CN))
}

gen.2layer = function(n,SigmaGiant,p1,p2){
	XY = mvrnorm(n,mu=rep(0,(p1+p2)),SigmaGiant);
	X = XY[,1:p1]; X = scale(X,center=T,scale=F);
	Y = XY[,(p1+1):(p1+p2)]; Y = scale(Y,center=T,scale=F);
	return(list(XY=XY,X=X,Y=Y))
}

# GenerateZ() generates data in the latent layer given variational parameters
GenerateZ = function(n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s){
  # Generate M and S
  M = matrix(0, n, q)
  S = M
  K = length(pi.m)
  # for(i in 1:n){
  #   for(j in 1:q){
  #     ran1 = runif(1,0,1)
  #     for(k in 1:K){
  #       if(ran1 <= cumsum(pi.m)[k]){
  #         M[i,j] = rnorm(1,mu.m[k],sig.m[k])
  #         break;
  #       }
  #     }
  #     
  #     ran2 = runif(1,0,1)
  #     for(k in 1:K){
  #       if(ran2 <= cumsum(pi.s[k])){
  #         S[i,j] = rnorm(1,mu.s[k],sig.s[k])
  #         break;
  #       }
  #     }
  #   }
  # }
  for(i in 1:n){
    for(j in 1:q){
      M[i,j] = sum(pi.m*rnorm(rep(1,K),mu.m,sig.m))
      S[i,j] = sum(pi.s*rnorm(rep(1,K),mu.s,sig.s))
    }
  }
  
  # generate Z
  M + matrix(rnorm(n*q),n,q)*exp(-S)
}

GenerateZ.bin = function(Z.bin, n, q, pi.m, mu.m, sig.m, pi.s, mu.s, sig.s){
  # Generate M and S
  M = matrix(0, n, q)
  S = M
  K = length(pi.m)
  
  for(i in 1:n){
    for(j in 1:q){
      bin = Z.bin[i,j]
      M[i,j] = sum(pi.m[,bin]*rnorm(rep(1,K),mu.m[,bin],sig.m[,bin]))
      S[i,j] = sum(pi.s[,bin]*rnorm(rep(1,K),mu.s[,bin],sig.s[,bin]))
    }
  }
  
  # generate Z
  M + matrix(rnorm(n*q),n,q)*exp(-S)
}










