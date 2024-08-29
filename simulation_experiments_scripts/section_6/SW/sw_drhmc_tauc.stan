functions{
	
	/*
	
	  Tri-diagonal Cholesky routines
	
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = [diagFirstLastElem,diagElem,...,diagElem,diagFirstLastElem]
		* first sup/sub-diagonal elements are all equal to offDiagElem
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	vector CIP_TriDiagChol_const1n(int n, real diagFirstLastElem, real diagElem, real offDiagElem){
		vector[2*n] L;
		real LlogDet;
		// first iteration
		L[1] = sqrt(diagFirstLastElem);
		LlogDet = log(L[1]);
		// iteration 2:n-1
		for ( t in 2:(n-1)){
			L[n+t-1] = offDiagElem/L[t-1];
			L[t] = sqrt(diagElem - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		// last iteration
		L[2*n-1] = offDiagElem/L[n-1];
		L[n] = sqrt(diagFirstLastElem - pow(L[2*n-1],2));
		LlogDet += log(L[n]);
		// done Cholesky
		
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = diagElem (diagElem has length n)
		* first sup/sub-diagonal elements are all equal to offDiagElem
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	
	
	vector CIP_TriDiagChol_diag_constod(vector diagElem, real offDiagElem){
		int n = rows(diagElem);
		vector[2*n] L;
		real LlogDet;
		L[1] = sqrt(diagElem[1]);
		LlogDet = log(L[1]);
		for ( t in 2:n){
			L[n+t-1] = offDiagElem/L[t-1];
			L[t] = sqrt(diagElem[t] - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = diagElem (diagElem has length n)
		* first sup/sub-diagonal elements are in offDiagElem ( offDiagElem has length n-1)
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	
	
	vector CIP_TriDiagChol(vector diagElem, vector offDiagElem){
		int n = rows(diagElem);
		vector[2*n] L;
		real LlogDet;
		L[1] = sqrt(diagElem[1]);
		LlogDet = log(L[1]);
		for ( t in 2:n){
			L[n+t-1] = offDiagElem[t-1]/L[t-1];
			L[t] = sqrt(diagElem[t] - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Solves L^T x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LT_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[n] = b[n]/L[n];
		// remaining solves
		for ( tt in 1:(n-1)){
			x[n-tt] = (b[n-tt] - x[n-tt+1]*L[2*n-tt])/L[n-tt];
		}
		return(x);
	}
	
	/*
		Solves L x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_L_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[1] = b[1]/L[1];
		// remaining solves
		for ( i in 2:n){
			x[i] = (b[i] - x[i-1]*L[n+i-1])/L[i];
		}
		return(x);
	}
	
	/*
		Solves L L^T x = G x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LLT_solve(vector L, vector b){
		return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b)));
	}

}

data{
	int T;
	// real y[T];
	array[T] real y;
	real alpha;
	real beta;
	real lambdamean;
	real lambdaprec;
	// real zmean[T-1];
	// real xmean[T];
	array[T - 1] real zmean;
	array[T] real xmean;
}



parameters{
	real lambda_bar;
	vector[T-1] z_bar;
	vector[T] x_bar;
	vector[T] tau_bar;
}

transformed parameters{
	real lambda;
	real sigma;
	real lprec;
	vector[T-1] z;
	vector[T] x;
	vector[T] tau;
	
	vector[2*(T-1)] Lz;
	vector[2*T] Lx;
	vector[T-1] expmz;
	vector[T] G4diag;
	vector[2*T] Ltau;
	vector[T] Pyy;
	
	// block 1
	lambda =  lambdamean + lambda_bar/sqrt(lambdaprec);
	sigma = exp(-0.5*lambda);
	lprec = exp(lambda);
	
	// block 2
	Lz = CIP_TriDiagChol_const1n(T-1,lprec+0.5,2.0*lprec+0.5,-lprec);
	z = CIP_TriDiagChol_LT_solve(Lz,z_bar);
	for(t in 1:(T-1)){
		z[t] += zmean[t];
	}
	
	// block 3
	Lx = CIP_TriDiagChol_const1n(T,lprec+0.5,2.0*lprec+0.5,-lprec);
	x = CIP_TriDiagChol_LT_solve(Lx,x_bar);
	for(t in 1:T){
		x[t] += xmean[t];
	}
	
	// block 4
	expmz = exp(-z);
	G4diag = exp(-x);
	for(t in 1:T){
		Pyy[t] = G4diag[t]*y[t];
	}
	
	G4diag[1] += expmz[1];
	for(t in 2:(T-1)){
		G4diag[t] += expmz[t-1]+expmz[t];
	}
	G4diag[T] += expmz[T-1];
	
	
	
	Ltau = CIP_TriDiagChol(G4diag,-expmz);
	tau = CIP_TriDiagChol_LT_solve(Ltau,tau_bar) + CIP_TriDiagChol_LLT_solve(Ltau,Pyy);
	
	
	
	
			 
}


model{
	target += alpha*lambda - beta*exp(lambda);
	
	for( t in 2:(T-1)) {
		target += normal_lpdf(z[t] | z[t-1], sigma);
	}
	
	for( t in 2:T){
		target += normal_lpdf(x[t] | x[t-1], sigma);
	}
	
	for( t in 2:T){
		target += normal_lpdf(tau[t] | tau[t-1],exp(0.5*z[t-1]));
	}
	
	for( t in 1:T){
		target += normal_lpdf(y[t] | tau[t], exp(0.5*x[t]));
	}
	
	target += -Lz[2*(T-1)] -Lx[2*T] -Ltau[2*T];
	
}




