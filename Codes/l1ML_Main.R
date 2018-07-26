source("PriorCalc.R")
source("l1LS_Main.R")
source("Obj.R")
source("StabilitySelection.R")

require(huge)
require(glasso)

# --------------------------------------------
# Args:
#	Y -- response variables in layer 2.
#	X -- predictor variables in layer 1.
#	initializer -- choose between Lasso, Ridge and CAPME, default is Lasso.
# 	lambda, rho -- tuning parameter for Lasso/Ridge, graphical Lasso; if NULL, sqrt(log(p)/n).
#	Screening -- if we would like to invoke the screening procedure, default is TRUE. 
#	correction -- the method we use for multiple testing correction, choose between "BH" and "bonferroni". Default is "BH". 
#	alpha -- control level for the multiple testing procedure, a value between 0 and 1. Default is 0.1. 
#	StabilizeTheta -- TRUE/FALSE, if we do stability selection on Theta. Default is TRUE.
#	rho.seq -- a sequence of parameter for stability selection, if NULL, seq(0.5,3,by=0.5)*sqrt(log(p2)/n).
#	B0,Theta0 -- initial value for B/Theta (mainly for debug use), default is NULL.
#	Trace = TRUE/FALSE; if we trace the change of B and Theta for every iteration, default is FALSE. 

l1ML_Main = function(Y,X,initializer="Lasso",lambda=NULL,rho=NULL,Screening=TRUE,correction="BH",alpha=0.1,StabilizeTheta=TRUE,rho.seq=NULL,B0=NULL,Theta0=NULL,VERBOSE=TRUE,TRACE=FALSE){
	n = nrow(X);
	if (n!=nrow(Y)){
		stop("The numbers of samples in Layer 1 and Layer 2 differ!")
	}
	p1 = ncol(X);
	p2 = ncol(Y);
	
	## if initial value for B and Theta are provided
	## mainly for debug/algorithm testing
	if (!is.null(B0) && !is.null(Theta0)){
		if (VERBOSE){
			cat("Initial value for B and Theta detected.\n")
		}
		Theta.initial = Theta0; B.initial=B0;
		# we restrict the subsequent alternating procedure on the support of the initial value
		skeletlon.hat = ifelse(abs(B.initial)>1e-6,1,0);
		Existing.edges = diag(1:p1) %*% skeleton.hat;
	}
	else{ # if either one of the initial value is not provided, we need to initialize
		
		## first we calculate the support of the regressions based on J-M
		if (Screening){		
			if (VERBOSE)
				cat("J-M in progress with ",as.character(correction),"correction.\n");
			skeleton.hat = PriorCalc(Ymat=Y,Xmat=X,alpha=alpha,correction=correction);
			Existing.edges = diag(1:p1) %*% skeleton.hat;

		}
		else {# no screening
			if (VERBOSE)
				cat("No screening.\n")
			skeleton.hat = array(1,c(p1,p2)); # assuming all directed edges possibly exist
			Existing.edges = diag(1:p1) %*% skeleton.hat;
		}
			
		
		if (initializer == "Lasso" | initializer == "Ridge"){
			if (VERBOSE){
				cat("Obtaining initial estimate from ", as.character(initializer), " ...",'\n')
			}	
			LS = l1LS_Main(Y,X,skeleton.hat=skeleton.hat,lambda=0.02,initializer=initializer);
			if (is.null(B0)){
				if (VERBOSE)
					cat("Using ", as.character(initializer), " estimates as the initializer for B",'\n')
				B.initial = LS$B0;	
			}
			else{
				if (VERBOSE)
					cat("Initial value of B is detected\n")
				B.initial = B0;
			}
			if (is.null(Theta0)){
				if (VERBOSE)
					cat("Use glasso on residuals as initializer for Theta\n")
				Theta.initial = LS$Theta0;
			}
			else{
				if (VERBOSE)
					cat("Initial value for Theta is detected\n")
				Theta.initial = Theta0;
			}
		}
		if (initializer=="CAPME"){  ## this is for debug use 
			library(capme)
			if (VERBOSE)
				cat("CAPME cross validating ... ")
			capmecv = cv.capme(fold=3,loss="likelihood",x=X,y=Y)
			if (VERBOSE)
				cat("DONE. Obtaining initial estimate from CAPME...",'\n')
			lambda.capme = capmecv$lambdaopt
			tau.capme = capmecv$tauopt
			capme.result = capme(x=X,y=Y,lambda=lambda.capme,tau=tau.capme)
			
			if (is.null(B0)){
				if (VERBOSE)
					cat("Using CAPME estimates as the initial value for B",'\n')
				B.initial = capme.result$Gammalist[[1]]
			}
			else if(!is.null(B0)){
				if (VERBOSE)
					cat("Initial value of B is detectd",'\n')
				B.initial = B0
			}
	
			if (is.null(Theta0)){
				if (VERBOSE)
					cat("Using CAPME estimates as the initial value for Theta",'\n')
				Theta.initial = capme.result$Omegalist[[1]]
				## Here we might need to threshold CAPME estimator a bit
				Theta.initial = ifelse(abs(Theta.initial)>1e-3,Theta.initial,0);
			}
			else if (!is.null(Theta0)){
				if (VERBOSE)
					cat("Initial value of Theta is detectd",'\n')
				Theta.initial = Theta0
			}
		}
	}# till this step, we finish the initialization before the alternating procedure
	
	# calculate the value of the objective function evaluated at the initial values
	Obj_initial = Obj(Y,X,Theta.initial,lambda,rho,B.initial)
	
	Objfunc = c(); iter = 0; CONVERGE=FALSE; update.counter=0;
	updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
	
	## initialization of the residual matrix
	Rmat = Y - X %*% B.initial;
	R = array(0,c(n,p2));  	
	for (j in 1:p2){
		for (k in (1:p2)[-j]){
			R[,j] = R[,j] + Theta.initial[j,k]*Rmat[,k];
		}  
		R[,j] = R[,j]/(2*Theta.initial[j,j]);  # we get r_j at this step
	}
	
	## record the result from every iteration
	## for algorithm testing/debugging purpose only
	if (TRACE){
		TRACE.est =  vector("list")
	}
	
	## Here we start with the alternating procedure
	B_new = B.initial; Theta_new = Theta.initial;
	while(!CONVERGE){
		iter = iter + 1;
		B_old = B_new ; Theta_old = Theta_new; 
		
		# Updating B
		output.list_B = foreach (j=1:p2)%dopar% {
			# the update of one column in B
			B_j = rep(0,p1);
			if (length(which(Existing.edges[,j]!=0))==0)
				B_j = B_j;
			if (length(which(Existing.edges[,j]!=0))==1){
				coef1 = lm((Y[,j]+R[,j])~X[,which(Existing.edges[,j]!=0)]+0)$coef;
				## this is the least square solution, but we need to apply a soft-thresholding
				coef1_st = max(0,coef1-lambda/(2*Theta_old[j,j]))*sign(coef1);
				B_j[which(Existing.edges[,j]!=0)] = coef1_st;
			}				
			if (length(which(Existing.edges[,j]!=0))>=2){
				temp = glmnet(X[,which(Existing.edges[,j]!=0)],(Y[,j]+R[,j]),intercept=FALSE);
				B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda/(2*Theta_old[j,j]),type="coefficients")[-1];
			}
			B_j;
		}
		for (j in 1:p2){
			B_new[,j] = output.list_B[[j]];
		}
		  		
		# Updating Theta
		if (iter >=10 | norm(B_new-B_old,"F")<0.1 | updateTheta == TRUE){
			updateTheta = TRUE; # once we start updating Theta, we just start from now on;
			ResMat = Y - X%*%B_new; 
			Theta_new = as.matrix(huge(ResMat,rho,method="glasso",verbose=FALSE)$icov[[1]]);
			update.counter = update.counter + 1;
		}
		
		# Updating r_i's to appear in the regression of the next iteration
		Rmat = Y - X %*% B_new;
		R = array(0,c(n,p2));
		# update R_i's
		for (i in 1:p2){
			for (k in (1:p2)[-i]){
				R[,i] = R[,i]+ Theta_new[i,k]*Rmat[,k];
			}
			R[,i] = R[,i]/(2*Theta_new[i,i]);
		}
		
		# record the value of the objective function
		Objfunc[iter] = Obj(Y,X,Theta_new,lambda,rho,B_new);
		if (TRACE){
			TRACE.est[[iter]] = vector("list",2)
			TRACE.est[[iter]][[1]] = B_new; TRACE.est[[iter]][[2]] = Theta_new; 
		}
		
		## check convergence
		if (iter == 1){
			Obj_diff = Objfunc[1]-Obj_initial;
		}
		else{
			Obj_diff = Objfunc[iter] - Objfunc[iter-1];
		}
		CONVERGE = (abs(Obj_diff)<1e-4);

		if (iter == 100){
			cat("Exceeds the maximum number of iterations allowed",'\n');
			break;
		}
		if (VERBOSE){
			cat(paste(c('iter=','Obj_diff=','lambda=','rho='),c(iter,round(Obj_diff,4),round(lambda,3),round(rho,3)),sep=" "),'\n');
		}
			
	
		## if converge, we start the refitting procedure
		if (CONVERGE){ 
			if (VERBOSE)
				cat("Theta was updated ",update.counter,"times",'\n')			

			if (VERBOSE)
				cat("Start refitting B with OLS ...\n")
			
			B.refit = array(0,c(p1,p2));
			output.list_refit = foreach (j = 1:p2) %dopar% {
				B.refit_j = rep(0,p1);
				if (length(which(abs(B_new[,j])>1e-6))>0){
					B.refit_j[which(abs(B_new[,j])>1e-6)]=lm((Y[,j]+R[,j])~X[,which(abs(B_new[,j])>1e-6)]+0)$coef;
				}
				B.refit_j;
			}	#finish refitting all the Bs
			for (j in 1:p2){
				B.refit[,j] = output.list_refit[[j]]
			}
			
			if (VERBOSE)
				cat("Start re-estimating Theta ...\n")
			ResMat.refit = Y - X %*% B.refit
			
			if (StabilizeTheta){
				require(glasso)
				if (is.null(rho.seq)){
					rho.seq = seq(from=0.5,to=3.0,by=0.5)*sqrt(log(p2)/n)
				}
				## the actual number of bootstrap is 2 times the given value ## 
				ProbMat = StabilitySelection(ResMat.refit,rho.seq,nbootstrap=25,VERBOSE=VERBOSE)
				WeightMat = 1 - ProbMat
				Theta.refit = glasso(var(ResMat.refit),rho*WeightMat,penalize.diagonal=FALSE)$wi
			}
			else{
				Theta.refit = as.matrix(huge(ResMat.refit,rho,method="glasso",verbose=FALSE)$icov[[1]])
			}
			
		}
	}
	B.est = B.refit; Theta.est = Theta.refit;
	BICvalue = BICfunc(Y,X,Theta.est,B.est);
	
	if (TRACE){
		return(list(B.screening=skeleton.hat, B.est=B.est,Theta.est=Theta.est,B.initial=B.initial,Theta.initial=Theta.initial,Obj.initial=Obj_initial,Objvals=Objfunc,BIC=BICvalue,BTheta.trace=TRACE.est))
	}
	if (!TRACE){
		return(list(B.screening=skeleton.hat, B.est=B.est,Theta.est=Theta.est,B.initial=B.initial,Theta.initial=Theta.initial,Obj.initial=Obj_initial,Objvals=Objfunc,BIC=BICvalue))
	}
}

# Values:
#	B.est, Theta,est -- the final estimates of the regression matrix and the Theta matrix;
#	B.initial, Theta.initial -- the initial value of B matrix and Theta matrix;
#	Obj.initial -- the initial value of the objective function;
#	Objvals -- the value change in the objective functions;
# 	BIC -- the BIC value given the final estimates;
#	BTheta.trace -- the value of B and Theta during each iteration, only available when TRACE=TRUE;
# -------------------------------------------
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

