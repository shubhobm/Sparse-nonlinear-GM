require(huge)
require(gdata)

StabilitySelection = function(X,rho.seq,nbootstrap=20,VERBOSE=TRUE){
	# Stability selection in the spirit of Meinshausen & Buhlman
	# Specifically for glasso
	
	# X is a n*p matrix; rho.seq is a sequence of tuning parameters
	# here, we don't do the randomized reweighting
	# the return value is a score (of length p(p-1)/2) for all possible edges 
	if (VERBOSE){
		cat("Stability Selection In Progress with", 2*nbootstrap, "bootstrapped samples \n")
	}
	n = nrow(X)
	p = ncol(X)
	halfsize = as.integer(n/2)
	freq = matrix(0,length(rho.seq),p*(p-1)/2)
	
	for (i in 1:nbootstrap){
		perm = sample(n)
		i1 = perm[1:halfsize]
		i2 = perm[(halfsize+1):n]
			
		glasso.est = huge(X[i1,],lambda=rho.seq,method="glasso",verbose=FALSE)$icov
		for (j in 1:length(rho.seq)){
			freq[j,] = freq[j,] + upperTriangle(abs(sign(as.matrix(glasso.est[[j]]))))
		}
		
		glasso.est = huge(X[i2,],lambda=rho.seq,method="glasso",verbose=FALSE)$icov
		for (j in 1:length(rho.seq)){
			freq[j,] = freq[j,] + upperTriangle(abs(sign(as.matrix(glasso.est[[j]]))))
		}
	}
	freq = freq/(2*nbootstrap)
	
	result = array(0,c(p,p))
	upperTriangle(result) = apply(freq,2,max)
	result = result + t(result)
	diag(result) = 1
	
	return(result)
}
