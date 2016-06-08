

combination.rule<- function(ret, gamma=1, superset=1:7, subset=NULL, detailed.output=FALSE, RHO.grid.size= 100, Kmax.init= 500, tail.cut.exp= 20){

	TT<- dim(ret)[1]
	N<- dim(ret)[2]
	D<- length(superset)
	D.set<- (1:7)[superset]

	get.atilde<- function(TT, N){
		z<- c( 1/sqrt(2*pi) * beta(TT/2-2, 3/2), 1/(TT-3)/(TT-5) ) ## first m=1, then m=2
		for(i in 1:(N-1)) z<- z* ( 1 + 1/(TT-i-3))
		list( a1tilde=z[1] , a2tilde=z[2])
	}
	atilde<- get.atilde(N=N, TT=TT)
	a1tilde<- atilde$a1tilde #
	a2tilde<- atilde$a2tilde #



	get.aRR<- function(TT, RHO.grid.size= 200, Kmax.init= 10000, tail.cut.exp= 30, TRACE=FALSE){ ## 30
	#Kmax.init #=# initial Kmax.init, not that 10000, size of the 'infinite' sum at beginning, if this is to small Kmax.init will be increased by the algorithm, if it is larger than required it will be reduced (using tail.cut.exp).
	#RHO.grid.size #=# creating rho grid with equidistant knots and distance 1/RHO.grid.size
	#tail.cut.exp	#=# exp, digits bound to cut series...
	#TRACE -- if TRUE print progress
		Kmax<- rep(Kmax.init, 3)
		RHO<- seq(0, 1, by=1/RHO.grid.size)
		RHO<- RHO[-c(1, length(RHO))]
		AFX<- array(, dim=c(length(RHO),3)  )
		## get.AFX
		get.AFX<- function(i.m, rho, MUtfu.fun, MX.fun){
			check<- TRUE
			while(check){
				K<- 0:Kmax[i.m]
				MUtfu<-  MUtfu.fun(K, rho)
				wm.WUtfu<- which.max(MUtfu)
				if( wm.WUtfu <= Kmax[i.m]) { ## if so, then everything okay, might reduce Kmax[i.m] [caution K starts at 0]
						bb<-MUtfu[wm.WUtfu] - MUtfu[wm.WUtfu:Kmax[i.m]] > tail.cut.exp
						if(!any(bb)){ ## then all==FALSE
							# double Kmax
							Kmax[i.m]<- Kmax[i.m]*2				#
						} else {
							Kmax[i.m]<- wm.WUtfu + which.max(bb)
							check<- FALSE
						}
				} else { # double Kmax[i.m]
					Kmax[i.m]<- Kmax[i.m]*2				#
				}
			}# check
			if(TRACE)print(Kmax[i.m])
			MX<-  MX.fun(K, rho)
			sum(exp(MUtfu)*MX)  ## CORRECT
		}# get.AFX
		## MUftu and MX
		MUTFU<- list()
		MXL<- list()
		MUTFU[[1]]<- function(K, rho) 2*K*log(rho) + (TT-5)/2*log(1-rho^2) - lbeta(K+1,(TT-1)/2)
		MXL[[1]]<- function(K, rho) 2/ ( (TT+ 2*K-3 )^2 * (TT+ 2*K-1 )   )	
		MUTFU[[2]]<- function(K, rho) 2*K*log(rho) + (TT-4)/2*log(1-rho^2) + lgamma( (2*K+TT-2)/2 ) - lgamma((TT-1)/2) - lgamma(K+1) 
		MXL[[2]]<- function(K, rho) 1/sqrt(2)/(TT+2*K-3)	
		MUTFU[[3]]<- function(K, rho)  2*K*log(rho) + (TT-3)/2*log(1-rho^2) + 2*lgamma( (2*K+TT-2)/2 ) - lgamma((TT-1)/2) - lgamma(K+1) - lgamma((2*K+TT-1)/2)  
		MXL[[3]]<- function(K, rho) 1/2	
		for(i.r in length(RHO):1 ){	## computation loop about RHO
			if(TRACE) print(i.r)
			rho<- RHO[i.r]
			for(i.m in 1:3)	AFX[i.r,i.m]<-get.AFX(i.m, rho=rho, MUtfu.fun=MUTFU[[i.m]], MX.fun=MXL[[i.m]])
		}
		fixed<- array(, dim=c( 3,2) )	## fixed moments
		## 22,11,12 | (0,1)
		fixed[1,1]<- 1/(TT-3)^2	# E(e_i^2)E(e_j^2)
		fixed[1,2]<- 1/(TT-3)/(TT-5) ## E(e_i^4)
		fixed[2,1]<- beta(TT/2-1, .5) / sqrt(2*pi) * 1/(TT-3)  # E(e_i)E(e_j^2)
		fixed[2,2]<- beta(TT/2-2, 3/2)/sqrt(pi)/sqrt(2) 
		fixed[3,1]<- beta(TT/2-1, .5)^2 / (2*pi)  # E(e_i)E(e_j)
		fixed[3,2]<- 1/(TT-3)  # E(e_i^2)
		AFX1<- c(fixed[1,2], rev(AFX[,1]) ,fixed[1,1], AFX[,1], fixed[1,2])
		AFX2<- c(fixed[2,2], rev(AFX[,2]) ,fixed[2,1], AFX[,2], fixed[2,2])
		AFX3<- c(fixed[3,2], rev(AFX[,3]) ,fixed[3,1], AFX[,3], fixed[3,2])
		RHOe<- c(-1, -rev(RHO), 0, RHO, 1)
		atilde22<- approxfun(RHOe, AFX1)
		atilde12<- approxfun(RHOe, AFX2)
		atilde11<- approxfun(RHOe, AFX3)
		list(atilde22=atilde22, atilde12=atilde12, atilde11=atilde11)
	}# end get.aRR

	aRR<- get.aRR(TT, RHO.grid.size, Kmax.init, tail.cut.exp, FALSE)

	a22tilde<- aRR[[1]] #approxfun(RHO, AFX22) ##S^{-2}, S^{-2}
	a12tilde<- aRR[[2]] #approxfun(RHO, AFX12) ##S^{-1}, S^{-2}
	a11tilde<- aRR[[3]] #approxfun(RHO, AFX11) ##S^{-1}, S^{-1}

tr<- function(z) sum(diag(z))
U<- function(w, m, S) as.numeric( t(w) %*% m - gamma/2 *t(w) %*% S %*% w )
ones<- rep.int(1, N)
I<- diag(N)

#a1,a2,...
a1<- TT/(TT-N-2)
a2<- N/TT
a3<- TT*(TT-2)/(TT-N-1)/(TT-N-4)
a4<- TT/(TT-3)
a5<- TT^2 * a2tilde
a6<- sqrt(TT)/sqrt(2) * beta(TT/2-1, .5) / sqrt(pi)
a7<- TT^(3/2) * a1tilde
## options... 2^D-1 ...
comb<- list()
cnames<- list()
k<- 1
for(i in 1:D) {
	tmp<- as.matrix(combn(D,i))
	for(j in 1:dim(tmp)[2])	{
		comb[[k]]<- tmp[, j]
		cnames[[k]]<- paste(D.set[comb[[k]]], sep="", collapse="")
		k<-k+1
	}#j
}#i
Dcomb<- length(comb)
cnames<- unlist(cnames)
names(comb)<- cnames

#### true A's
get.Abvec<- function(mu, Sigma, adjust=0){
## adjust== 0 normal plug-in,  adjust==1 , sharphe-ratio squared plug-in, adjust==2 psi-squared adj. as KZ 3fund.
	Sigma.inv<- solve(Sigma)
	sinv2<- 1/diag(Sigma)
	Sinv2 <- diag(sinv2)
	sinv1<- sqrt(sinv2)
	Sinv1<- diag(sinv1)

	RR<- cov2cor(Sigma)
	Upsilon22<- TT^2 * Sigma * a22tilde(RR)
	Upsilon11<- TT * Sigma * a11tilde(RR)
	Upsilon12<- TT^(3/2) * Sigma * a12tilde(RR)

#### ADJUSTED squared sharpe ratio [adjust==1]
	sharpe.squared<- as.numeric( mu %*% Sigma.inv %*% mu )	
	ibeta<- function(x,a,b) pbeta(x,a,b) * beta(a,b) ## incomplete beta, see kanzhou2007
	sharpe.squared.adj<- ((TT-N-2)*sharpe.squared - N)/TT + 2*(sharpe.squared^(N/2)*(1+ sharpe.squared)^(-(TT-2)/2))/TT/ibeta(sharpe.squared/(1+sharpe.squared), N/2, (TT-N)/2)
## 
#### ADJUSTED psi squared [adjust==2]
	mug<- t(ones) %*% Sigma.inv %*% mu / as.numeric( t(ones) %*% Sigma.inv %*% ones )
	psi.squared<- as.numeric( sharpe.squared - mug * t(ones) %*% Sigma.inv %*% mu )
	psi.squared.adj<- ((TT-N-1)*psi.squared - (N-1))/TT + 2*(psi.squared)^((N-1)/2)*(1+psi.squared)^(-(TT-2)/2)/( TT * ibeta(psi.squared/(1+psi.squared),(N-1)/2,(TT-N+1)/2))
##################
	b<- numeric(8)
#	b[1]<- a1 * sharpe.squared.adj #t(mu) %*% Sigma.inv %*% mu
	if(adjust==1){
		b[1]<- a1 * sharpe.squared.adj #t(mu) %*% Sigma.inv %*% mu
	} else if(adjust==2) {
		b[1]<- a1 * (psi.squared.adj  + t(ones) %*% Sigma.inv %*% ones * mug^2 ) #t(mu) %*% Sigma.inv %*% mu
	} else {
		b[1]<- a1 * sharpe.squared #t(mu) %*% Sigma.inv %*% mu
	}
	b[2]<- a1 * t(ones) %*% Sigma.inv %*% mu ## the same for all elements
	b[3]<-  t(mu) %*% mu
	b[4]<-  t(ones) %*% mu
	b[5]<- a4 * t(mu) %*% Sinv2 %*% mu
	b[6]<- a4 * t(sinv2) %*% mu
	b[7]<- a6 * t(sinv1) %*% mu
	b[8]<- a1 * ones %*% mu * t(ones) %*% Sigma.inv %*% mu
	b<- b[superset]
	# A matrix:
#########
	A<- array(, dim=c(8,8))
	if(adjust==1){
		A[1,1]<- a1*a3*(sharpe.squared.adj + a2)
		A[1,2]<- a1*a3* t(mu) %*% Sigma.inv %*% ones
		A[2,2]<- a1*a3* t(ones) %*% Sigma.inv %*% ones
	} else if(adjust==2) {
		A[1,1]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones * (  (psi.squared.adj + a2)/ as.numeric(t(ones) %*% Sigma.inv %*% ones)  +  mug^2 ) #t(mu) %*% Sigma.inv %*% mu
		A[1,2]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones * mug
		A[2,2]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones 
	} else {
		A[1,1]<- a1*a3*(sharpe.squared + a2)
		A[1,2]<- a1*a3* t(mu) %*% Sigma.inv %*% ones
		A[2,2]<- a1*a3* t(ones) %*% Sigma.inv %*% ones
	}
	A[1,3]<- a1*t(mu) %*% mu + a1*tr(Sigma)/TT
	A[1,4]<- a1* t(mu) %*% ones
	A[1,5]<- a5* (t(mu) %*% Sinv2 %*% mu + a2)
	A[1,6]<- a5* t(mu) %*% sinv2
	A[1,7]<- a7* t(mu) %*% sinv1
	A[1,8]<- a1*a3* ( ones %*% mu * t(mu) %*% Sigma.inv %*% ones +a2)
	A[2,3]<- a1 * t(ones) %*% mu
	A[2,4]<- a1*N
	A[2,5]<- a5 * t(mu) %*% sinv2
	A[2,6]<- a5 * t(ones) %*% sinv2
	A[2,7]<- a7* t(ones) %*% sinv1
	A[2,8]<- a1 * a3 * t(ones) %*% mu * t(ones) %*% Sigma.inv %*% ones 
	A[3,3]<- t(mu) %*% Sigma %*% mu + tr(Sigma %*% Sigma)/TT
	A[3,4]<- t(mu) %*% Sigma %*% ones
	A[3,5]<- a4*t(mu) %*% Sigma %*% Sinv2 %*% mu + a4*tr(Sigma %*% Sinv2 %*% Sigma)/TT
	A[3,6]<- a4 * t(mu) %*% Sigma %*% sinv2
	A[3,7]<- a6 * t(mu) %*% Sigma %*% sinv1
	A[3,8]<- a1 * ( t(ones) %*% mu )^2 + a1* t(ones) %*% Sigma %*% ones /TT
	A[4,4]<- t(ones) %*% Sigma %*% ones
	A[4,5]<- a4* t(ones) %*% Sigma %*% Sinv2 %*% mu
	A[4,6]<- a4* t(ones) %*% Sigma %*% sinv2 
	A[4,7]<- a6 * t(ones) %*% Sigma %*% sinv1
	A[4,8]<- a1 * t(ones) %*% mu * N 
	A[5,5]<- t(mu) %*% Sinv2 %*% Upsilon22 %*% Sinv2 %*% mu + tr( Sinv2 %*% Upsilon22 %*% Sinv2 %*% Sigma)/TT
	A[5,6]<- t(mu) %*% Sinv2 %*% Upsilon22 %*% sinv2
	A[5,7]<- t(mu) %*% Sinv2 %*% Upsilon12 %*% sinv1 
	A[5,8]<- a5 *(  t(ones) %*% mu * t(mu) %*% sinv2 + t(ones) %*% Sigma %*% sinv2 /TT )
	A[6,6]<- t(sinv2) %*% Upsilon22 %*% sinv2 
	A[6,7]<- t(sinv2) %*% Upsilon12 %*% sinv1 
	A[6,8]<- a5 *  t(ones) %*% mu * t(ones) %*% sinv2 
	A[7,7]<- t(sinv1) %*% Upsilon11 %*% sinv1 
	A[7,8]<- a7 *  t(ones) %*% mu * t(ones) %*% sinv1 
	A[8,8]<- a1*a3*(  (t(ones) %*% mu)^2 + t(ones) %*% Sigma %*% ones / TT ) * t(ones) %*% Sigma.inv %*% ones
	## upper
	A[lower.tri(A)]<-t(A)[lower.tri(A)]
	A<- A[superset, superset]
	## w.vec
	vec<- array(, dim=c(8,N) )
	vec[1,]<- Sigma.inv %*% mu
	vec[2,]<- Sigma.inv %*% ones
	vec[3,]<- mu
	vec[4,]<- ones
	vec[5,]<- Sinv2 %*% mu
	vec[6,]<- sinv2
	vec[7,]<- sinv1
	vec[8,]<- as.numeric(t(ones) %*% mu) * Sigma.inv %*% ones
	vec<- vec[superset,]
	return(list(A,b,vec))
}	#get.Abcvec




if(is.null(subset)) {
	A.index<- 1:Dcomb
	A.names<- cnames 
}	else {
	A.index<- which(apply( sapply( subset, grepl, x=cnames) , 1, all))
	A.names<- cnames[A.index] ## subset
}
A.len<- length(A.names)




adj.index<- grep("1",A.names )
A.names.adj<- gsub("1", "1'",A.names[adj.index])
adj.len<- length(adj.index)

adj2.index<- grep("1",A.names )
A.names.adj2<- gsub("1", "1''",A.names[adj2.index])
adj2.len<- length(adj2.index)
AW.len<- length(A.names)+ adj.len + adj2.len



## weights
W.hat<- array(, dim=c(AW.len, N))
#DELTA.hat<- array(, dim=c(AW.len,D))

	DELTA.hat<- array(, dim=c( A.len, D))	## on est
	DELTA.adj.hat<- array(, dim=c( adj.len, D))	## on est
	DELTA.adj2.hat<- array(, dim=c( adj2.len, D))	## on est
	##benchmark
	mu.hat<- apply(ret, 2, mean)	## exess return
	Sigma.hat<- cov(ret) * (TT-1)/TT
	Abvec.hat<- get.Abvec(mu.hat, Sigma.hat, 0)
	A.hat<- Abvec.hat[[1]]
	b.hat<- Abvec.hat[[2]]
	vec.hat<- Abvec.hat[[3]]
	scales.hat<-  1/apply(vec.hat,1, sum)
	Scales.hat<- diag(scales.hat)
	vec.hat.scaled.hat<- vec.hat*scales.hat
	## adj
	Abvec.adj.hat<- get.Abvec(mu.hat, Sigma.hat, 1)
	A.adj.hat<- Abvec.adj.hat[[1]]
	b.adj.hat<- Abvec.adj.hat[[2]]
	vec.adj.hat<- Abvec.adj.hat[[3]]
	scales.adj.hat<-  1/apply(vec.adj.hat,1, sum)
	Scales.adj.hat<- diag(scales.adj.hat)
	vec.adj.hat.scaled.hat<- vec.adj.hat*scales.adj.hat
	## adj
	Abvec.adj2.hat<- get.Abvec(mu.hat, Sigma.hat, 2)
	A.adj2.hat<- Abvec.adj2.hat[[1]]
	b.adj2.hat<- Abvec.adj2.hat[[2]]
	vec.adj2.hat<- Abvec.adj2.hat[[3]]
	scales.adj2.hat<-  1/apply(vec.adj2.hat,1, sum)
	Scales.adj2.hat<- diag(scales.adj2.hat)
	vec.adj2.hat.scaled.hat<- vec.adj2.hat*scales.adj2.hat
	for(i.D in 1:A.len){
		delta.hat<- solve(A.hat[comb[[A.index[i.D]]] ,comb[[A.index[i.D]]]] %*% Scales.hat[comb[[A.index[i.D]]] ,comb[[A.index[i.D]]]], b.hat[comb[[A.index[i.D]]]])/gamma
		DELTA.hat[ i.D,comb[[A.index[i.D]]]]<- delta.hat ## scales estimated
		W.hat[i.D,]<- delta.hat %*% vec.hat.scaled.hat[comb[[A.index[i.D]]],]
	}	
## adjusted
	for(i.D in seq_along(adj.index) ){
		delta.adj.hat<- solve(A.adj.hat[comb[[A.index[adj.index[i.D]]]] ,comb[[A.index[adj.index[i.D]]]]] %*% Scales.adj.hat[comb[[A.index[adj.index[i.D]]]] ,comb[[A.index[adj.index[i.D]]]]], b.adj.hat[comb[[A.index[adj.index[i.D]]]]])/gamma
		DELTA.adj.hat[i.D,comb[[A.index[adj.index[i.D]]]]]<- delta.adj.hat ## scales estimated
		W.hat[A.len + i.D,]<- delta.adj.hat %*% vec.adj.hat.scaled.hat[comb[[A.index[adj.index[i.D]]]],]
	}#dev.off()
	for(i.D in seq_along(adj2.index) ){
		delta.adj2.hat<- solve(A.adj2.hat[comb[[A.index[adj2.index[i.D]]]] ,comb[[A.index[adj2.index[i.D]]]]] %*% Scales.adj2.hat[comb[[A.index[adj2.index[i.D]]]] ,comb[[A.index[adj2.index[i.D]]]]], b.adj2.hat[comb[[A.index[adj2.index[i.D]]]]])/gamma
		DELTA.adj2.hat[i.D,comb[[A.index[adj2.index[i.D]]]]]<- delta.adj2.hat ## scales estimated
		W.hat[A.len +adj.len+ i.D,]<- delta.adj2.hat %*% vec.adj2.hat.scaled.hat[comb[[A.index[adj2.index[i.D]]]],]
	}#dev.off()
	## bench marks
	all.names<- c(A.names, A.names.adj, A.names.adj2)
	dimnames(W.hat)<- list( all.names ,dimnames(ret)[[2]]  )
	delta.hat<- rbind(DELTA.hat, DELTA.adj.hat, DELTA.adj2.hat)
	dimnames(delta.hat)<- list( all.names ,D.set)
	V<- vec.hat.scaled.hat
	dimnames(V)<- list( D.set ,dimnames(ret)[[2]]  )
	if(!detailed.output) return(W.hat[order(all.names), ]) else return( list(w=W.hat[order(all.names), ], delta=delta.hat[order(all.names), ],  V=t(V) ) )
}#8fund



combination.rule.restriction<- function(ret, HC, h0, rule, gamma=1, detailed.output=FALSE, RHO.grid.size= 100, Kmax.init= 500, tail.cut.exp= 20){
	
	superset<- rule
	subset<- rule

		TT<- dim(ret)[1]
	N<- dim(ret)[2]
	D<- length(superset)
	D.set<- (1:7)[superset]

	get.atilde<- function(TT, N){
		z<- c( 1/sqrt(2*pi) * beta(TT/2-2, 3/2), 1/(TT-3)/(TT-5) ) ## first m=1, then m=2
		for(i in 1:(N-1)) z<- z* ( 1 + 1/(TT-i-3))
		list( a1tilde=z[1] , a2tilde=z[2])
	}
	atilde<- get.atilde(N=N, TT=TT)
	a1tilde<- atilde$a1tilde #
	a2tilde<- atilde$a2tilde #



	get.aRR<- function(TT, RHO.grid.size= 200, Kmax.init= 10000, tail.cut.exp= 30, TRACE=FALSE){ ## 30
	#Kmax.init #=# initial Kmax.init, not that 10000, size of the 'infinite' sum at beginning, if this is to small Kmax.init will be increased by the algorithm, if it is larger than required it will be reduced (using tail.cut.exp).
	#RHO.grid.size #=# creating rho grid with equidistant knots and distance 1/RHO.grid.size
	#tail.cut.exp	#=# exp, digits bound to cut series...
	#TRACE -- if TRUE print progress
		Kmax<- rep(Kmax.init, 3)
		RHO<- seq(0, 1, by=1/RHO.grid.size)
		RHO<- RHO[-c(1, length(RHO))]
		AFX<- array(, dim=c(length(RHO),3)  )
		## get.AFX
		get.AFX<- function(i.m, rho, MUtfu.fun, MX.fun){
			check<- TRUE
			while(check){
				K<- 0:Kmax[i.m]
				MUtfu<-  MUtfu.fun(K, rho)
				wm.WUtfu<- which.max(MUtfu)
				if( wm.WUtfu <= Kmax[i.m]) { ## if so, then everything okay, might reduce Kmax[i.m] [caution K starts at 0]
						bb<-MUtfu[wm.WUtfu] - MUtfu[wm.WUtfu:Kmax[i.m]] > tail.cut.exp
						if(!any(bb)){ ## then all==FALSE
							# double Kmax
							Kmax[i.m]<- Kmax[i.m]*2				#
						} else {
							Kmax[i.m]<- wm.WUtfu + which.max(bb)
							check<- FALSE
						}
				} else { # double Kmax[i.m]
					Kmax[i.m]<- Kmax[i.m]*2				#
				}
			}# check
			if(TRACE)print(Kmax[i.m])
			MX<-  MX.fun(K, rho)
			sum(exp(MUtfu)*MX)  ## CORRECT
		}# get.AFX
		## MUftu and MX
		MUTFU<- list()
		MXL<- list()
		MUTFU[[1]]<- function(K, rho) 2*K*log(rho) + (TT-5)/2*log(1-rho^2) - lbeta(K+1,(TT-1)/2)
		MXL[[1]]<- function(K, rho) 2/ ( (TT+ 2*K-3 )^2 * (TT+ 2*K-1 )   )	
		MUTFU[[2]]<- function(K, rho) 2*K*log(rho) + (TT-4)/2*log(1-rho^2) + lgamma( (2*K+TT-2)/2 ) - lgamma((TT-1)/2) - lgamma(K+1) 
		MXL[[2]]<- function(K, rho) 1/sqrt(2)/(TT+2*K-3)	
		MUTFU[[3]]<- function(K, rho)  2*K*log(rho) + (TT-3)/2*log(1-rho^2) + 2*lgamma( (2*K+TT-2)/2 ) - lgamma((TT-1)/2) - lgamma(K+1) - lgamma((2*K+TT-1)/2)  
		MXL[[3]]<- function(K, rho) 1/2	
		for(i.r in length(RHO):1 ){	## computation loop about RHO
			if(TRACE) print(i.r)
			rho<- RHO[i.r]
			for(i.m in 1:3)	AFX[i.r,i.m]<-get.AFX(i.m, rho=rho, MUtfu.fun=MUTFU[[i.m]], MX.fun=MXL[[i.m]])
		}
		fixed<- array(, dim=c( 3,2) )	## fixed moments
		## 22,11,12 | (0,1)
		fixed[1,1]<- 1/(TT-3)^2	# E(e_i^2)E(e_j^2)
		fixed[1,2]<- 1/(TT-3)/(TT-5) ## E(e_i^4)
		fixed[2,1]<- beta(TT/2-1, .5) / sqrt(2*pi) * 1/(TT-3)  # E(e_i)E(e_j^2)
		fixed[2,2]<- beta(TT/2-2, 3/2)/sqrt(pi)/sqrt(2) 
		fixed[3,1]<- beta(TT/2-1, .5)^2 / (2*pi)  # E(e_i)E(e_j)
		fixed[3,2]<- 1/(TT-3)  # E(e_i^2)
		AFX1<- c(fixed[1,2], rev(AFX[,1]) ,fixed[1,1], AFX[,1], fixed[1,2])
		AFX2<- c(fixed[2,2], rev(AFX[,2]) ,fixed[2,1], AFX[,2], fixed[2,2])
		AFX3<- c(fixed[3,2], rev(AFX[,3]) ,fixed[3,1], AFX[,3], fixed[3,2])
		RHOe<- c(-1, -rev(RHO), 0, RHO, 1)
		atilde22<- approxfun(RHOe, AFX1)
		atilde12<- approxfun(RHOe, AFX2)
		atilde11<- approxfun(RHOe, AFX3)
		list(atilde22=atilde22, atilde12=atilde12, atilde11=atilde11)
	}# end get.aRR

	aRR<- get.aRR(TT, RHO.grid.size, Kmax.init, tail.cut.exp, FALSE)

	a22tilde<- aRR[[1]] #approxfun(RHO, AFX22) ##sinv2, sinv2
	a12tilde<- aRR[[2]] #approxfun(RHO, AFX12) ##sinv1, sinv2
	a11tilde<- aRR[[3]] #approxfun(RHO, AFX11) ##sinv1, sinv1

tr<- function(z) sum(diag(z))
U<- function(w, m, S) as.numeric( t(w) %*% m - gamma/2 *t(w) %*% S %*% w )
ones<- rep.int(1, N)
I<- diag(N)

#a1,a2,...
a1<- TT/(TT-N-2)
a2<- N/TT
a3<- TT*(TT-2)/(TT-N-1)/(TT-N-4)
a4<- TT/(TT-3)
a5<- TT^2 * a2tilde
a6<- sqrt(TT)/sqrt(2) * beta(TT/2-1, .5) / sqrt(pi)
a7<- TT^(3/2) * a1tilde
## options... 2^D-1 ...
comb<- list()
cnames<- list()
k<- 1
for(i in 1:D) {
	tmp<- as.matrix(combn(D,i))
	for(j in 1:dim(tmp)[2])	{
		comb[[k]]<- tmp[, j]
		cnames[[k]]<- paste(D.set[comb[[k]]], sep="", collapse="")
		k<-k+1
	}#j
}#i
Dcomb<- length(comb)
cnames<- unlist(cnames)
names(comb)<- cnames

#### true A's
get.Abvec<- function(mu, Sigma, adjust=0){
## adjust== 0 normal plug-in,  adjust==1 , sharphe-ratio squared plug-in, adjust==2 psi-squared adj. as KZ 3fund.
	Sigma.inv<- solve(Sigma)
	sinv2<- 1/diag(Sigma)
	Sinv2<- diag(sinv2)
	sinv1<- sqrt(sinv2)
	Sinv1<- diag(sinv1)

	RR<- cov2cor(Sigma)
	Upsilon22<- TT^2 * Sigma * a22tilde(RR)
	Upsilon11<- TT * Sigma * a11tilde(RR)
	Upsilon12<- TT^(3/2) * Sigma * a12tilde(RR)

#### ADJUSTED squared sharpe ratio [adjust==1]
	sharpe.squared<- as.numeric( mu %*% Sigma.inv %*% mu )	
	ibeta<- function(x,a,b) pbeta(x,a,b) * beta(a,b) ## incomplete beta, see kanzhou2007
	sharpe.squared.adj<- ((TT-N-2)*sharpe.squared - N)/TT + 2*(sharpe.squared^(N/2)*(1+ sharpe.squared)^(-(TT-2)/2))/TT/ibeta(sharpe.squared/(1+sharpe.squared), N/2, (TT-N)/2)
## 
#### ADJUSTED psi squared [adjust==2]
	mug<- t(ones) %*% Sigma.inv %*% mu / as.numeric( t(ones) %*% Sigma.inv %*% ones )
	psi.squared<- as.numeric( sharpe.squared - mug * t(ones) %*% Sigma.inv %*% mu )
	psi.squared.adj<- ((TT-N-1)*psi.squared - (N-1))/TT + 2*(psi.squared)^((N-1)/2)*(1+psi.squared)^(-(TT-2)/2)/( TT * ibeta(psi.squared/(1+psi.squared),(N-1)/2,(TT-N+1)/2))
##################
	b<- numeric(8)
#	b[1]<- a1 * sharpe.squared.adj #t(mu) %*% Sigma.inv %*% mu
	if(adjust==1){
		b[1]<- a1 * sharpe.squared.adj #t(mu) %*% Sigma.inv %*% mu
	} else if(adjust==2) {
		b[1]<- a1 * (psi.squared.adj  + t(ones) %*% Sigma.inv %*% ones * mug^2 ) #t(mu) %*% Sigma.inv %*% mu
	} else {
		b[1]<- a1 * sharpe.squared #t(mu) %*% Sigma.inv %*% mu
	}
	b[2]<- a1 * t(ones) %*% Sigma.inv %*% mu ## the same for all elements
	b[3]<-  t(mu) %*% mu
	b[4]<-  t(ones) %*% mu
	b[5]<- a4 * t(mu) %*% Sinv2 %*% mu
	b[6]<- a4 * t(sinv2) %*% mu
	b[7]<- a6 * t(sinv1) %*% mu
	b[8]<- a1 * ones %*% mu * t(ones) %*% Sigma.inv %*% mu
	b<- b[superset]
	# A matrix:
#########
	A<- array(, dim=c(8,8))
	if(adjust==1){
		A[1,1]<- a1*a3*(sharpe.squared.adj + a2)
		A[1,2]<- a1*a3* t(mu) %*% Sigma.inv %*% ones
		A[2,2]<- a1*a3* t(ones) %*% Sigma.inv %*% ones
	} else if(adjust==2) {
		A[1,1]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones * (  (psi.squared.adj + a2)/ as.numeric(t(ones) %*% Sigma.inv %*% ones)  +  mug^2 ) #t(mu) %*% Sigma.inv %*% mu
		A[1,2]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones * mug
		A[2,2]<- a1 *a3* t(ones) %*% Sigma.inv %*% ones 
	} else {
		A[1,1]<- a1*a3*(sharpe.squared + a2)
		A[1,2]<- a1*a3* t(mu) %*% Sigma.inv %*% ones
		A[2,2]<- a1*a3* t(ones) %*% Sigma.inv %*% ones
	}
	A[1,3]<- a1*t(mu) %*% mu + a1*tr(Sigma)/TT
	A[1,4]<- a1* t(mu) %*% ones
	A[1,5]<- a5* (t(mu) %*% Sinv2 %*% mu + a2)
	A[1,6]<- a5* t(mu) %*% sinv2
	A[1,7]<- a7* t(mu) %*% sinv1
	A[1,8]<- a1*a3* ( ones %*% mu * t(mu) %*% Sigma.inv %*% ones +a2)
	A[2,3]<- a1 * t(ones) %*% mu
	A[2,4]<- a1*N
	A[2,5]<- a5 * t(mu) %*% sinv2
	A[2,6]<- a5 * t(ones) %*% sinv2
	A[2,7]<- a7* t(ones) %*% sinv1
	A[2,8]<- a1 * a3 * t(ones) %*% mu * t(ones) %*% Sigma.inv %*% ones 
	A[3,3]<- t(mu) %*% Sigma %*% mu + tr(Sigma %*% Sigma)/TT
	A[3,4]<- t(mu) %*% Sigma %*% ones
	A[3,5]<- a4*t(mu) %*% Sigma %*% Sinv2 %*% mu + a4*tr(Sigma %*% Sinv2 %*% Sigma)/TT
	A[3,6]<- a4 * t(mu) %*% Sigma %*% sinv2
	A[3,7]<- a6 * t(mu) %*% Sigma %*% sinv1
	A[3,8]<- a1 * ( t(ones) %*% mu )^2 + a1* t(ones) %*% Sigma %*% ones /TT
	A[4,4]<- t(ones) %*% Sigma %*% ones
	A[4,5]<- a4* t(ones) %*% Sigma %*% Sinv2 %*% mu
	A[4,6]<- a4* t(ones) %*% Sigma %*% sinv2 
	A[4,7]<- a6 * t(ones) %*% Sigma %*% sinv1
	A[4,8]<- a1 * t(ones) %*% mu * N 
	A[5,5]<- t(mu) %*% Sinv2 %*% Upsilon22 %*% Sinv2 %*% mu + tr( Sinv2 %*% Upsilon22 %*% Sinv2 %*% Sigma)/TT
	A[5,6]<- t(mu) %*% Sinv2 %*% Upsilon22 %*% sinv2
	A[5,7]<- t(mu) %*% Sinv2 %*% Upsilon12 %*% sinv1 
	A[5,8]<- a5 *(  t(ones) %*% mu * t(mu) %*% sinv2 + t(ones) %*% Sigma %*% sinv2 /TT )
	A[6,6]<- t(sinv2) %*% Upsilon22 %*% sinv2 
	A[6,7]<- t(sinv2) %*% Upsilon12 %*% sinv1 
	A[6,8]<- a5 *  t(ones) %*% mu * t(ones) %*% sinv2 
	A[7,7]<- t(sinv1) %*% Upsilon11 %*% sinv1 
	A[7,8]<- a7 *  t(ones) %*% mu * t(ones) %*% sinv1 
	A[8,8]<- a1*a3*(  (t(ones) %*% mu)^2 + t(ones) %*% Sigma %*% ones / TT ) * t(ones) %*% Sigma.inv %*% ones
	## upper
	A[lower.tri(A)]<-t(A)[lower.tri(A)]
	A<- A[superset, superset]
	## w.vec
	vec<- array(, dim=c(8,N) )
	vec[1,]<- Sigma.inv %*% mu
	vec[2,]<- Sigma.inv %*% ones
	vec[3,]<- mu
	vec[4,]<- ones
	vec[5,]<- Sinv2 %*% mu
	vec[6,]<- sinv2
	vec[7,]<- sinv1
	vec[8,]<- as.numeric(t(ones) %*% mu) * Sigma.inv %*% ones
	vec<- vec[superset,]
	return(list(A,b,vec))
}	#get.Abcvec




if(is.null(subset)) {
	A.index<- 1:Dcomb
	A.names<- cnames 
}	else {
	A.index<- which(apply( sapply( subset, grepl, x=cnames) , 1, all))
	A.names<- cnames[A.index] ## subset
}
A.len<- length(A.names)




adj.index<- grep("1",A.names )
A.names.adj<- gsub("1", "1'",A.names[adj.index])
adj.len<- length(adj.index)

adj2.index<- grep("1",A.names )
A.names.adj2<- gsub("1", "1''",A.names[adj2.index])
adj2.len<- length(adj2.index)
AW.len<- length(A.names)+ adj.len + adj2.len



## weights
W.hat<- array(, dim=c(AW.len, N))
#DELTA.hat<- array(, dim=c(AW.len,D))

	DELTA.hat<- array(, dim=c( A.len, D))	## on est
	DELTA.adj.hat<- array(, dim=c( adj.len, D))	## on est
	DELTA.adj2.hat<- array(, dim=c( adj2.len, D))	## on est
	##benchmark
	mu.hat<- apply(ret, 2, mean)	## exess return
	Sigma.hat<- cov(ret) * (TT-1)/TT
	Abvec.hat<- get.Abvec(mu.hat, Sigma.hat, 0)
	A.hat<- Abvec.hat[[1]]
	b.hat<- Abvec.hat[[2]]
	vec.hat<- Abvec.hat[[3]]
	scales.hat<-  rep(1, D) #1/apply(vec.hat,1, sum)
	Scales.hat<- diag(scales.hat)
	vec.hat.scaled.hat<- vec.hat*scales.hat
	## adj
	Abvec.adj.hat<- get.Abvec(mu.hat, Sigma.hat, 1)
	A.adj.hat<- Abvec.adj.hat[[1]]
	b.adj.hat<- Abvec.adj.hat[[2]]
	vec.adj.hat<- Abvec.adj.hat[[3]]
	scales.adj.hat<-  1/apply(vec.adj.hat,1, sum)
	Scales.adj.hat<- diag(scales.adj.hat)
	vec.adj.hat.scaled.hat<- vec.adj.hat*scales.adj.hat
	## adj
	Abvec.adj2.hat<- get.Abvec(mu.hat, Sigma.hat, 2)
	A.adj2.hat<- Abvec.adj2.hat[[1]]
	b.adj2.hat<- Abvec.adj2.hat[[2]]
	vec.adj2.hat<- Abvec.adj2.hat[[3]]
	scales.adj2.hat<-  1/apply(vec.adj2.hat,1, sum)
	Scales.adj2.hat<- diag(scales.adj2.hat)
	vec.adj2.hat.scaled.hat<- vec.adj2.hat*scales.adj2.hat
	for(i.D in 1:A.len){
		Ainv<- solve(A.hat[comb[[A.index[i.D]]] ,comb[[A.index[i.D]]]] )
		Cinv<- solve(Scales.hat[comb[[A.index[i.D]]] ,comb[[A.index[i.D]]]])
		bb<- b.hat[comb[[A.index[i.D]]]]
		delta.hat<- as.numeric( (Cinv %*% Ainv  - Cinv %*% Ainv %*%  HC  %*% solve( t(HC) %*% Ainv %*% HC  ) %*% t(HC) %*% Ainv    ) %*% bb /gamma + Cinv %*% Ainv %*% HC %*% solve( t(HC) %*% Ainv %*% HC  ) %*% h0 )
		DELTA.hat[ i.D,comb[[A.index[i.D]]]]<- delta.hat ## scales estimated
		W.hat[i.D,]<- delta.hat %*% vec.hat.scaled.hat[comb[[A.index[i.D]]],]
	}	
## adjusted
	for(i.D in seq_along(adj.index) ){
		Ainv<- solve(A.adj.hat[comb[[A.index[adj.index[i.D]]]] ,comb[[A.index[adj.index[i.D]]]]] )
		Cinv<- solve(Scales.adj.hat[comb[[A.index[adj.index[i.D]]]] ,comb[[A.index[adj.index[i.D]]]]])
		bb<- b.adj.hat[comb[[A.index[adj.index[i.D]]]]]
		delta.adj.hat<- as.numeric( (Cinv %*% Ainv  - Cinv %*% Ainv %*%  HC  %*% solve( t(HC) %*% Ainv %*% HC  ) %*% t(HC) %*% Ainv    ) %*% bb /gamma + Cinv %*% Ainv %*% HC %*% solve( t(HC) %*% Ainv %*% HC  ) %*% h0 )
		DELTA.adj.hat[i.D,comb[[A.index[adj.index[i.D]]]]]<- delta.adj.hat ## scales estimated
		W.hat[A.len + i.D,]<- delta.adj.hat %*% vec.adj.hat.scaled.hat[comb[[A.index[adj.index[i.D]]]],]
	}#dev.off()
	for(i.D in seq_along(adj2.index) ){
		Ainv<- solve(A.adj2.hat[comb[[A.index[adj2.index[i.D]]]] ,comb[[A.index[adj2.index[i.D]]]]] )
		Cinv<- solve(Scales.adj2.hat[comb[[A.index[adj2.index[i.D]]]] ,comb[[A.index[adj2.index[i.D]]]]])
		bb<- b.adj2.hat[comb[[A.index[adj2.index[i.D]]]]]
		delta.adj2.hat<- as.numeric( (Cinv %*% Ainv  - Cinv %*% Ainv %*%  HC  %*% solve( t(HC) %*% Ainv %*% HC  ) %*% t(HC) %*% Ainv    ) %*% bb /gamma + Cinv %*% Ainv %*% HC %*% solve( t(HC) %*% Ainv %*% HC  ) %*% h0 )
		DELTA.adj2.hat[i.D,comb[[A.index[adj2.index[i.D]]]]]<- delta.adj2.hat ## scales estimated
		W.hat[A.len +adj.len+ i.D,]<- delta.adj2.hat %*% vec.adj2.hat.scaled.hat[comb[[A.index[adj2.index[i.D]]]],]
	}#dev.off()
	## bench marks
	all.names<- paste("r:",c(A.names, A.names.adj, A.names.adj2), sep="")
	dimnames(W.hat)<- list( all.names ,dimnames(ret)[[2]]  )
	delta.hat<- rbind(DELTA.hat, DELTA.adj.hat, DELTA.adj2.hat)
	dimnames(delta.hat)<- list( all.names ,D.set)
	V<- vec.hat.scaled.hat
	dimnames(V)<- list( D.set ,dimnames(ret)[[2]]  )
	if(!detailed.output) return(W.hat[order(all.names), ]) else return( list(w=W.hat[order(all.names), ], delta=delta.hat[order(all.names), ],  V=t(V) ) )
}#8fund



