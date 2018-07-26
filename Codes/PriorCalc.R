source("lasso_inference.r")


PriorCalc = function(Ymat,Xmat,alpha=0.1,adjust=TRUE,correction="BH"){
	# Ymat are the nodes in layer 2, n*p2 (or n*q)
	# Xmat are the nodes in layer 1, n*p1
	# lambda is the tuning parameter for lasso
	
	p1 = ncol(Xmat)
	p2 = ncol(Ymat)
	n = nrow(Xmat)
	
	p.val = array(0,c(p1,p2))
	output.list_prior = foreach(j = 1:p2)%dopar%{
		f = SSLasso(Xmat,Ymat[,j],verbose=FALSE,intercept=FALSE)
		p.val_j = f$pvals
		p.val_j # the p-values of the jth regression
	}
	
	# collect the result
	for (j in 1:p2){
		p.val[,j] = output.list_prior[[j]] 
	}
	
	if (adjust)
		p.val = array(p.adjust(as.vector(p.val),method=correction),c(p1,p2))
	
	Skeleton.hat = (p.val<= alpha)
		
	return(Skeleton.hat)
}

# Note here, Nodes.hat are arranged as an p1 * p2 matrix, with column denoting existing variables corresponding to a specific node in Layer 2
