library(gdata)
loglik = function(Y,X,Theta,B)
{ 
	# this function calculates the loglikelihood 
	n = nrow(X)
	p1 = ncol(X)
	p2 = ncol(Y)
	
	expterm = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))
	l = (n/2)*log(det(Theta))-expterm/2
	return(l)
}


Obj = function(Y,X,Theta,lambda,rho,B)
{
	# this function calculates the objective function
	n = nrow(X)
	p1 = ncol(X)
	p2 = ncol(Y)

	dataterm = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n
	
	B.pen = lambda*sum(abs(B))
	Theta.pen = rho*(sum(abs(Theta))-sum(diag(Theta)))
	
	Objfunc = dataterm - log(det(Theta)) + B.pen + Theta.pen
	
	return(Objfunc)
}

BICfunc = function(Y,X,Theta,B)
{	#BIC = -log(det(Theta)) + tr(S%*%Theta) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n + (log(n)/n)*(nonzeros in upperTriangle(Theta) + nonzeros in B)
	#And we prefer smaller BIC value
	n = nrow(X)
	PENs = log(n)/n*(sum(abs(upperTriangle(Theta))>1e-6) + sum(abs(B)>1e-6))
	FIT = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n - log(det(Theta)) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n 
	BIC = FIT + PENs
	return(BIC)
}



		
		
				

