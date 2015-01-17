# return the gradient vector of C
GradiantC = function(fai, t.lst){
	d = length(fai) # d = t - 1
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai[i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	gradient = numeric(d)
	ones = rep(1,d)
	denom = rowSums(K) + ones
	for ( i in 1:d){
		gradient[i] = sum(rowSums(-1 * K * t.lst[[i]]) / denom)
	}
	gradient
}
debug(GradiantC)
# return list of t
t.gen = function(d){
	t.lst = list()
	t.lst[[d]] = matrix(rep(1:d,d),ncol = d, nrow = d,byrow=T)
	left.mask = matrix(rep(0,d^2),ncol = d, nrow = d)
	left.mask[2:d,1:(d-1)] = diag(rep(1,d-1))
	t.lst[[d]][upper.tri(left.mask)] = 0
	for ( i in 1:(d-1)){
		t.lst[[d-i]] = left.mask%*%t.lst[[d-i+1]]
		diag(t.lst[[d-i]]) = c(rep(0,i),1:(d-i))
		t.lst[[d-i]][upper.tri(left.mask)] = 0
	}
	t.lst
}

#testing 
T4 = t.gen(3)
GradiantC(c(10,3,1),T4)
# end of testing

GradiantVerify3 = function(fai){
	a = -1 * sum(fai)
	b = -1 * (fai[1] + 2*fai[2])
	gradient = numeric(2)
	gradient[1] = 1/(1+exp(a)+exp(b))*(-1*exp(a)-exp(b))
	gradient[2] = 1/(1+exp(-1*fai[2]))*(-1*exp(-1*fai[2])) - (exp(a) + 2*exp(b))/(1+exp(a)+exp(b))
	gradient
}
# testing
fai3 = c(4,2)
T3 = t.gen(2)
GradiantC(fai3,T3)
GradiantVerify3(fai3)
# end of testing


HessianVerify3 = function(fai){
	a = -1 * sum(fai)
	b = -1 * (fai[1] + 2*fai[2])
	hessian = matrix(ncol=2,nrow=2)
	# gradient = numeric(2)
	# gradient[1] = 1/(1+exp(a)+exp(b))*(-1*exp(a)-exp(b))
	# gradient[2] = 1/(1+exp(-1*fai[2]))*(-1*exp(-1*fai[2])) - (exp(a) + 2*exp(b))/(1+exp(a)+exp(b))
	# gradient
	tt = -1/(1+exp(a)+exp(b))^2 * (-1*exp(a)-exp(b)) * (-1*exp(a)-2*exp(b))
	tt = tt + 1/(1+exp(a)+exp(b)) * (exp(a)+2*exp(b))
	hessian[2,1] = tt
	hessian[1,2] = tt
	aa = -1/(1+exp(a)+exp(b))^2 * (-1*exp(a)-exp(b))^2
	aa = aa + 1/(1+exp(a)+exp(b)) * (exp(a)+exp(b))
	gg = -1*exp(-1*fai[2])^2/(1+exp(-1*fai[2]))^2 + exp(-1*fai[2])/(1 + exp(-1*fai[2]))
	gg = gg -1/(1+exp(a)+exp(b))^2 * (-1*exp(a)-2*exp(b)) * (-1*exp(a)-2*exp(b))
	gg = gg + 1/(1+exp(a)+exp(b)) * (exp(a)+4*exp(b))
	hessian[1,1] = aa
	hessian[2,2] = gg
	hessian
}
HessianVerify3(fai3)


HessianC = function(fai, t.lst){
	# same as gradient
	d = length(fai) # d = t - 1
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai[i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	# same as gradient end
	hessian = matrix(ncol=d, nrow=d)
	ones = rep(1,d)
	left.mask = matrix(rep(0,d^2),ncol = d, nrow = d)
	left.mask[1:(d-1),2:d] = diag(rep(1,d-1))
	Qq = matrix(ncol = d, nrow = d)
	W = matrix(ncol = d, nrow = d)
	for (j in 1:d){
			Qq[,j] = rowSums(-1*K*t.lst[[j]])
	}
	A = rowSums(K) + ones # A = denom !
	for (i in d:1){ # starting from the last value
		if (i == d){ 
			A.elinv = A^-1
		} else {
			A.elinv[d-i] = 0
		}
		B = rowSums(-1 * K * t.lst[[i]]) # need to change d in each step
		if(i !=d ){
			Qq[d-i,] = 0
		}
		for (j in 1:d){
			W[,j] = rowSums(t.lst[[j]]*t.lst[[i]]*K)
		}
		hessian[i,] = (-1*A.elinv^(2)*B)%*%Qq + A.elinv%*%W
	}
	#hessian[upper.tri(hessian)] = t(hessian[lower.tri(hessian)])
	hessian
}
debug(HessianC)
fai3= c(10,6)
HessianC(fai3,T3)
HessianVerify3(fai3)
fai5 = c(3,2,3,1)
T5 = t.gen(4)
GHC(fai5,T5)

############################ Calaculating Gradiant and Hessian for log(C) ############################ 
# arguments: fai in step t; t.lst returned by t.gen(d)
# result:	 a list containing gradient and hessian of log(C)
GHC = function(fai,t.lst){
	d = length(fai) # d = t - 1
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai[i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	gradient = numeric(d)
	ones = rep(1,d)
	denom = rowSums(K) + ones
	B = matrix(ncol=d,nrow=d)
	for (i in 1:d){
		B[,i] = rowSums(-1 * K * t.lst[[i]])
	}
	for ( i in 1:d){
		gradient[i] = sum(B[,i] / denom)
	}

	hessian = matrix(ncol=d, nrow=d)
	Qq = matrix(ncol = d, nrow = d)
	W = matrix(ncol = d, nrow = d)
	for (j in 1:d){
			Qq[,j] = rowSums(-1*K*t.lst[[j]])
	}
	for (i in d:1){ # starting from the last value
		if (i == d){ 
			A.elinv = denom^-1
		} else {
			A.elinv[d-i] = 0
		}
		if(i != d ){
			Qq[d-i,] = 0
		}
		for (j in 1:d){
			W[,j] = rowSums(t.lst[[j]]*t.lst[[i]]*K)
		}
		hessian[i,] = (-1*A.elinv^(2)*B[,i])%*%Qq + A.elinv%*%W
	}
	return(list(gradient=gradient,hessian=hessian))
}

# for efficiency reasons both ordering and ranking should be matrix
NewtonSolve = function(ordering=NULL,ranking=NULL,pi0,fai.init,verbose=T,epsilon=1e-8,limit=5000){
	if(verbose){
		start.time = Sys.time()
	}
	# arguments checking
	if (is.null(ordering) & is.null(ranking)){
		print("Either ranking or ordering should be provided!")
		return(NULL)
	}
	if (!is.null(ranking)){
		ordering = RankingToOrdering(ranking)
	}
	# end of checking
	n = ncol(ordering)
	nobs = nrow(ordering)
	t.lst = t.gen(n-1)
	fai.lst = list()
	fai.lst[[1]] = fai.init
	if (verbose == T){
		GHitr = list()
	}
	fai.coeff = rowSums(apply(ordering,1,WeightGivenPi,pi0))
	while(TRUE){
		itr = length(fai.lst)
#browser()
		fai.res = GHC(fai.lst[[itr]],t.lst)
		fai.gradient = -1*nobs*fai.res$gradient - fai.coeff
		fai.info = solve(nobs*fai.res$hessian)
		fai.next = fai.lst[[itr]] + fai.info%*%fai.gradient
		fai.lst[[itr+1]]=fai.next
		if (verbose==T){
			GHitr[[itr]] = list(fai.gradient=fai.gradient,fai.info=fai.info)
		}
		if ( all( abs(fai.next - fai.lst[[itr]]) < epsilon ) ){
			break
		}
		if (itr > limit){
			print(paste("does not converge in",limit,"iterations"))
			break
		}
	}
	fai.sig = solve(nobs * GHC(fai.lst[[length(fai.lst)]],t.lst)$hessian)
	
	
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai.lst[[length(fai.lst)]][i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	log_likelihood = -1*nobs*sum(log(rowSums(K)+1)) - fai.coeff %*% fai.lst[[length(fai.lst)]]
	
	
	
	
	if (verbose==T){
		end.time = Sys.time()
		return(list(fai.est=fai.lst,fai.sig=fai.sig,GHitr=GHitr,log_likelihood=log_likelihood,time.len=list(start.time=start.time,end.time=end.time)))
	}else {
		return(list(fai.est=fai.lst,fai.sig=fai.sig,log_likelihood=log_likelihood))
	}
}
debug(NewtonSolve)
perm = matrix(c(c(1,3,2,4),c(2,1,4,3),c(4,3,1,2),c(4,3,1,2)),nrow=4,byrow=T)
pi0=c(1,2,3,4)
NewtonSolve(ordering=perm,pi0=pi0,fai.init=c(1,1,1))

RankingToOrdering = function(ranking){
	obj = 1:ncol(ranking)
	ordering = matrix(obj[ranking],ncol=dim(ranking)[2],nrow=dim(ranking)[1])
	ordering
}

OrderingToRanking = function(ordering){
	ranking = t(apply(ordering,1,order))
	ranking
}
OrderingToRanking(perm)	
perm

exp(sum(log(rowSums(K)+1)))

# AllSolve = function(ordering=NULL,ranking=NULL,count,fai.init,verbose=T,epsilon=1e-8,limit=5000){
	# # try the most probable pi0
	# if(is.null(ranking)){
		# ranking = OrderingToRanking(ordering)
	# }
	# n=length(ordering[1,])
	# avg_rank = count %*% ranking;
	# pi0.est = list()
	# pi0.est[1] = list(sort(ordering[1,])[order(avg_rank)])

	# solve.result = NewtonSolve(ordering=ordering,ranking=ranking,count=count,pi0=pi0.est[[1]],fai.init=fai.init,verbose=F)
	# log_likelihood = list()
	# log_likelihood[[1]] = solve.result$log_likelihood

	# while(TRUE){
		# # verify using adjacent swap
		# target_log_likelihood = log_likelihood[[length(log_likelihood)]]
		# verify.trial = list()
		# pi0.trial.tmp = list()
		# verify.likelihood = numeric(n-1)
		# for (i in 2:n){
			# pi0.trial = pi0.est[[length(pi0.est)]]
			# pi0.trial[(i-1):i] = pi0.trial[i:(i-1)]
			# pi0.trial.tmp[[i-1]] = pi0.trial
			# verify.trial[i-1] = list(NewtonSolve(ordering=ordering,ranking=ranking,count=count,pi0=pi0.trial,fai.init=fai.init,verbose=F))
		# #browser()			
			# verify.likelihood[i-1] = verify.trial[[i-1]]$log_likelihood
		# }
		# if(max(verify.likelihood) > target_log_likelihood){
			# target_log_likelihood = max(verify.likelihood)
			# pos = which(verify.likelihood == target_log_likelihood)
			# pi0.est[[length(pi0.est)+1]] = pi0.trial.tmp[pos]
			# solve.result = verify.trial[pos]
		# } else {
			# break
		# }
	# }
		
	# return(list(solve.result = solve.result,pi0.est=pi0.est))
# }
	

# AllSolve(ordering=order.sample,fai.init=c(0.8,0.5,0.2,0.1,0.1))


AllSolve = function(ordering=NULL,ranking=NULL,count,fai.init,verbose=T,epsilon=1e-8,limit=5000){
	swap = function(arr,ind1,ind2){
		tmp = arr[ind1]
		arr[ind1] = arr[ind2]
		arr[ind2] = tmp
		arr
	}
	# try the most probable pi0
	if(is.null(ranking)){
		ranking = OrderingToRanking(ordering)
	}
	n=length(ordering[1,])
	avg_rank = count %*% ranking;
	pi0.est = list()
	pi0.est[1] = list(sort(ordering[1,])[order(avg_rank)])

	solve.result = NewtonSolve(ordering=ordering,ranking=ranking,count=count,pi0=pi0.est[[1]],fai.init=fai.init,verbose=F)
	log_likelihood = list()
	log_likelihood[[1]] = solve.result$log_likelihood

	while(TRUE){
		# verify using adjacent swap
		target_log_likelihood = log_likelihood[[length(log_likelihood)]]
		verify.trial = list()
		pi0.trial.tmp = list()
		K = choose(n,2)
		verify.likelihood = numeric(K)
		ind1=1
		ind2=2
		for (i in 1:K){
			pi0.trial = swap(pi0.est[[length(pi0.est)]],ind1,ind2)
			ind2 = ind2+1
			if(ind2 > n){
				ind1 = ind1 + 1
				ind2 = ind1 + 1
			}
			pi0.trial.tmp[[i]] = pi0.trial
		
			verify.trial[i] = list(NewtonSolve(ordering=ordering,ranking=ranking,count=count,pi0=pi0.trial,fai.init=fai.init,verbose=F))
			verify.likelihood[i] = verify.trial[[i]]$log_likelihood
		}
		if(max(verify.likelihood) > target_log_likelihood){
			target_log_likelihood = max(verify.likelihood)
			pos = which(verify.likelihood == target_log_likelihood)
			pi0.est[[length(pi0.est)+1]] = pi0.trial.tmp[pos]
			solve.result = verify.trial[pos]
		} else {
			break
		}
	}
		
	return(list(solve.result = solve.result,pi0.est=pi0.est))
}


# obtain C(fai)
faiC = function(fai, t.lst){
	d = length(fai) # d = t - 1
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai[i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	ones = rep(1,d)
	denom = exp(sum(log(rowSums(K) + ones)))
	denom
}

moments.est = function(ranking=NULL,ordering=NULL,count,size){
	# arguments checking
	if (is.null(ordering) & is.null(ranking)){
		print("Either ranking or ordering should be provided!")
		return(NULL)
	}
	if (!is.null(ordering)){
		ranking = OrderingToRanking(ordering)
	}
	# estimating pi0
	avg_rank = count %*% ranking;
	pi0.est = sort(ordering[1,])[order(avg_rank)]
	nfai = ncol(ranking)-1
	nobs = sum(count)
	# construct equation set
	hio = order(count/nobs)
	hi = sample(2:nrow(ordering),size=size,replace=T)
	probhi = sort(count/nobs,decreasing=T)[hi]
	orderhi = ordering[hi,]
	logodd = numeric(nfai)
	equmat = matrix(ncol=nfai,nrow=nfai)
	for (i in 1:nfai){
		logodd[i] = log(probhi[i]/probhi[i+1])
		kwj = WeightGivenPi(pi0.est,orderhi[i+1,])
		kwi = WeightGivenPi(pi0.est,orderhi[i,])
		equmat[i,] = kwj - kwi
	}
	fai.est = lm(logodd~equmat-1)
	fai.est
}
undebug(moments.est)
moments.est(ordering=capa.ordering,count=capa.dat$Count,size=1000)

########################## mixture model ##########################
# pi0 each row is a pi0.init (ordering)
# clu number of clusters
# fai rach row is a fai.init
MixtureSolve= function(ordering=NULL,ranking=NULL,count,fai.init,pi0.init,clu,verbose=T,epsilon=1e-8,limit=5000){
    # preparations 
    if(is.null(ranking)){
        ranking = OrderingToRanking(ordering)
    }
    tt=length(ordering[1,]) # number of objects
    n = sum(count)	# total observations
    distinctn = nrow(ranking)	# total distinct observations
    z = matrix(ncol = distinctn,nrow = clu)
    # return values
    pi0.est = matrix(ncol=tt, nrow=clu)
    pi0.est.last = matrix(ncol=tt, nrow=clu)
    p.vec = numeric(clu)
    p.vec.last = numeric(clu)
    fai = matrix(ncol=clu,nrow=tt-1)
    fai.last = matrix(ncol=clu,nrow=tt-1)
	fai.sig = list()
    # E step
	loopind=0
    while(loopind < 20){
		loopind = loopind+1
        z.tmp = matrix(ncol = distinctn,nrow = clu)
        for (i in 1:clu){
            pi.now = pi0.init[i,]
            t2 = apply(ordering,1,WeightGivenPi,pi.now) 
            z.tmp[i,] = fai.init[i,] %*% t2
        }
        z.tmp = exp(-1*z.tmp)
        sums = colSums(z.tmp)
        for (i in 1:distinctn){
            z[,i] = z.tmp[,i]/sums[i]
        }

        # M step
        p.vec = z %*% count
        p.vec = p.vec/sum(p.vec)
        for ( i in 1:clu){
			# fai.est=fai.est,fai.sig=fai.sig,pi0.est=pi0.est
            count.clu = z[i,] * count
            solve.clu = AllSolve(ordering=ordering,count=count.clu,fai.init=fai.init[i,])
            pi0.est[i,] = solve.clu$pi0.est
			# nrstep = length(solve.clu$solve.result$fai.est)
			
            fai[,i] = solve.clu$fai.est
			fai.sig[[i]] = solve.clu$fai.sig
        }
		# break?
		if (loopind == 1){
			p.vec.last = p.vec
			pi0.est.last = pi0.est
			fai.last = fai
		} else if (loopind < limit) {
			cond1 = all( p.vec.last - p.vec < epsilon)
			cond2 = all( pi0.est.last - pi0.est < epsilon)
			cond3 = all( fai.last - fai < epsilon)
			if (cond1 && cond2 && cond3){
				break
			}
		} else {
			print(paste("Algorithm did not converge in",limit,"iterations"))
			return(NULL)
		}
    } # inf loop
	return (list(p.vec=p.vec,pi0.est=pi0.est,fai=fai,fai.sig=fai.sig,iteration=loopind))
}

ordering = matrix(nrow=15,ncol=10)
for (i in 1:nrow(ordering)){
	ordering[i,] = sample(1:10,10,replace=F)
}
w = 9:1
fai = wTofai(w)
pi0 = 1:10
count=rep(1,15)
fai = rbind(fai,fai)
pi0 = rbind(pi0,rev(pi0))
debug(Mixture.solve)
Mixture.solve(ordering=ordering,count=count,fai = fai,clu=2,pi0 = pi0)


########################## simulating some data ##########################
library(permute)
sethow = how(observed = T)
perm5 = allPerms(5,sethow)
pi01 = sample(1:5)
w.true1 = rep(1,4) #exp(-0.5*c(-1,1,2,3))
dist.true1 = apply(perm5,1,KwDist,pi01,w.true1)
probs.true1 = exp(-1*dist.true1)/sum(exp(-1*dist.true1))

pi02 = sample(1:5)
w.true2 = rep(2,4) #exp(-1*c(-1,1,2,3))
dist.true2 = apply(perm5,1,KwDist,pi02,w.true2)
probs.true2 = exp(-1*dist.true2)/sum(exp(-1*dist.true2))

p.vec = c(1,0)

probs.mix = probs.true1*p.vec[1] + probs.true2*p.vec[2]

ind.samp = sample(x=nrow(perm5),size=5000,prob=probs.mix,replace=T)
ind.samp = table(ind.samp)
fai.init = moments.est(ordering =perm5.s,count = ind.samp,pi0.est = pi01,complete = F,size=500)$coef

fai.guess = rbind(fai.init,fai.init)
pi0.guess = rbind(pi01,sample(pi02))

perm5.s = as.numeric(names(ind.samp))
perm5.s = perm5[perm5.s,]

Kt = MixtureSolve(ordering=perm5.s,ranking=NULL,count=ind.samp,fai=fai.guess,clu=2,pi0=pi0.guess,verbose=F,epsilon=1e-8,limit=5000)
Kt = AllSolve(ordering=perm5.s,count=ind.samp,fai.init = fai.init,pi0.init = pi01,verbose=F,epsilon=1e-8,limit=5000)
Kt = NewtonSolve(ordering = perm5.s,count=ind.samp,fai.init = fai.init,pi0 = pi01,verbose=F)

Kt$pi0.est
Kt$p.vec
faiTow(Kt$fai[,1])
w.true1
faiTow(Kt$fai[,2])
w.true2

w.est1 = estw(Kt$fai[,1],Kt$fai.sig[[1]])
plot(main="Estimated w",x=(1:4),y=w.est1[[1]] + 2*w.est1[[2]],ylab="w",xlab="index",ylim=c(0,2),type="p")
points(x=factor(1:4),y=w.est1[[1]] - 2*w.est1[[2]])
points(x=1:4,y=w.true1,col="red",pch=3)
points(x=1:4,y=w.est1[[1]],col="blue",pch=4)
legend("topright",inset=.08, c("True w","Estimated w"), col=c("red", "blue"), pch=3:4)

w.est2 = estw(Kt$fai[,2],Kt$fai.sig[[2]])
plot(main="Estimated w",x=(1:4),y=w.est2[[1]] + 2*w.est2[[2]],ylab="w",xlab="index",ylim=c(0,5),type="p")
points(x=factor(1:4),y=w.est2[[1]] - 2*w.est2[[2]])
points(x=1:4,y=w.true2,col="red",pch=3)
points(x=1:4,y=w.est2[[1]],col="blue",pch=4)
legend("topright",inset=.08, c("True w","Estimated w"), col=c("red", "blue"), pch=3:4)

#NewtonSolve(ordering=ordering.w,count=count.clu,pi0=c(5,4,3,2,1),fai.init = fai.init.w,verbose=F)
#AllSolve(ordering=ordering.w,count=count.clu,fai.init = fai.init.w,verbose=F)
library(pmr)
dset = cbind(as.data.frame(perm5.s),as.numeric(ind.samp))
wdbm(dset,dtype="tau")


# try to use built-in optimizer

# NewtonSolve(ordering = perm5.s,count=ind.samp,fai.init = fai.init,pi0 = pi01,verbose=F)

# fai.coeff = apply(ordering,1,WeightGivenPi,pi0)%*%count

# fai.res = GHC(fai.lst[[itr]],t.lst)
# fai.gradiant = -1*nobs*fai.res$gradiant - fai.coeff
# fai.info = solve(nobs*fai.res$hessian)

# log_likelihood = -1*nobs*sum(log(rowSums(lastitr$K)+1)) - as.numeric(fai.coeff) %*% fai.lst[[length(fai.lst)]]

# ObjectiveLL = function(fai){

# nlminb(start, objective, gradient = NULL, hessian = NULL, ..., control = list(), lower = -Inf, upper = Inf)

# constraint: logical vector indicating which fai should be included ( TRUE = not set to zero)
NewtonSolveCon = function(ordering=NULL,count=NULL,pi0,fai.init,verbose=TRUE,tor=0,epsilon=1e-8,limit=5000,constraint=NULL){
	n = ncol(ordering)
	nobs = sum(count)
	if(verbose){
		start.time = Sys.time()
	}
	# arguments checking
	if(is.null(count)){
		count = rep(1,nrow(ordering))
	}
	if(!is.null(constraint)){
		contraint = rep(1,n-1)
	}
	if(class(constraint) != "logical"){
		print("constraint should be a LOGICAL vector")
		return(NULL)
	}
	# end of checking
	t.lst = t.gen(n-1)
	fai.lst = list()
	# enforce constraint on initial value
	fai.init[!constraint] = 0
	fai.lst[[1]] = fai.init
	if (verbose == T){
		GHitr = list()
	}
	fai.coeff = apply(ordering,1,WeightGivenPi,pi0)%*%count
	fai.coeff.con = fai.coeff[constraint]
	while(TRUE){
		itr = length(fai.lst)
		fai.res = GHC(fai.lst[[itr]],t.lst)
		grad = fai.res$gradiant[constraint]
		fai.gradiant = -1*nobs*grad - fai.coeff.con
		fai.info = solve(nobs*fai.res$hessian[constraint,constraint])
		fai.next.con = rep(0,n-1)
		fai.next.con[constraint] = fai.info%*%fai.gradiant
		fai.next = fai.lst[[itr]] + fai.next.con
		# handle 0 in each step
		zeroind = which(fai.next < 0)
		fai.next[zeroind] = tor
		# handle end
		fai.lst[[itr+1]]=fai.next
		if (verbose==T){
			GHitr[[itr]] = list(fai.gradiant=fai.gradiant,fai.info=fai.info)
		}

		case1 = all( abs(fai.next - fai.lst[[itr]]) < epsilon )
		if( is.na(case1)){
			print("doesn't work")
			return(NULL)
		}
		if ( case1 ){
			break
		}
		case2 = FALSE
		for(k in 1:itr){
			if(all(fai.next == fai.lst[[k]])){
				print("loop detected result might be suboptimal")
				case2=TRUE
				break
			}
		}
		if(case2 == T){
			break
		}
		if (itr > limit){
			print(paste("does not converge in",limit,"iterations"))
			break
		}
	}
	lastitr = GHC(fai.lst[[length(fai.lst)]],t.lst)
	fai.sig = solve(nobs * lastitr$hessian)
	# calculate coerced likelihood
	log_likelihood = -1*nobs*sum(log(rowSums(lastitr$K)+1)) - as.numeric(fai.coeff) %*% fai.lst[[length(fai.lst)]]
	if (verbose==T){
		end.time = Sys.time()
		return(list(fai.est=fai.lst,fai.sig=fai.sig,GHitr=GHitr,log_likelihood=log_likelihood,time.len=list(start.time=start.time,end.time=end.time)))
	}else {
		return(list(fai.est=fai.lst,fai.sig=fai.sig,log_likelihood=log_likelihood))
	}
}



# testing the algorithm
library(permute)
sethow = how(observed = T)
perm5 = allPerms(5,sethow)
pi0 = c(1,2,3,4,5)
w.true = rep(0.5,4)
dist.true = apply(perm5,1,KwDist,pi0,w.true)
probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
truesample = sample(1:gamma(6),prob=probs.true,size=10000,replace=T)
truesample = table(truesample)


testorder = matrix(nrow=nrow(perm5),ncol=ncol(perm5))
for (i in 1:nrow(perm5)){
	testorder[i,] = perm5[i,]
}
class(testorder)

testdat = new("RankData",nobj=5L,nobs=10000,ndistinct=120,ordering=testorder,ranking=RankingToOrdering(testorder),count=as.numeric(truesample))


debug(NewtonSolveCon)
test0 = NewtonSolveCon(constraint = c(rep(FALSE,3),rep(TRUE,1)),tor = 0,ordering=perm5,count=truesample,pi0=1:5,fai.init=rep(1,4),verbose=T,limit=5000)


test0$log_likelihood
test1 = NewtonSolve(tor = 0,ordering=perm5,count=truesample,pi0=1:5,fai.init=rep(1,4),verbose=FALSE,limit=5000)
test1$log_likelihood

testctrl = new("RankControl")


setwd("D:/Dropbox/study/Final Year Project/datasets/sushi")
source(file="../../simulation/functions.r")
source(file="../../simulation/OOP.r")
sushia = read.table("sushi3a.5000.10.order",skip=1)
sushia.ranking = as.matrix(sushia[,3:12])
sushia.ranking = sushia.ranking+1
debug(RankData)
sushidata = RankData(ranking=sushia.ranking)

sushisuff=ExtractSuff(sushidata,1:10)


all(sushia.ordering == sushidata@ordering)
all(sushia.ranking == sushidata@ranking)
all(sushia.count == sushidata@count)


testsuffdat = new("SuffData",nobj=5L,nobs=10000,ndistinct=120,fai.coeff=c(5646,14038,23172,30344),count=as.numeric(truesample),t.lst=t.gen(4))

testinit= new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(1:5),clu=1L)
testctrl = new("RankControl")
NewtonSolveCon(testsuffdat,testinit,c(TRUE,TRUE,FALSE,TRUE),testctrl)


NewtonStepSolve(testsuffdat,testinit,testctrl)


testinit = new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(c(2,3,4,1,5)),clu=1L)
AllSolve(testdat,testinit,testctrl)






setGeneric(name="MixtureSolve",
	def=function(dat,init,ctrl){standardGeneric("MixtureSolve")}
)

setMethod(
	f="MixtureSolve",
	signature=c(dat="RankData",init="RankInit",ctrl="RankControl"),
	definition=function(dat,init,ctrl){
		
		tt=dat@nobj # number of objects
		n = dat@nobs	# total observations
		distinctn = dat@ndistinct	# total distinct observations
		clu = init@clu
		count = dat@count
		z = matrix(ncol = distinctn,nrow = clu)
		# return values
		pi0.est = list()
		pi0.est.last = list()
		p.vec = numeric(clu)
		p.vec.last = numeric(clu)
		fai = list()
		fai.last = list()
		fai.sig = list()

		loopind=0
		while(TRUE){
			loopind = loopind+1
			# handle loopind=1
			if (loopind == 1){
				pi0.est = init@pi0.init
			}
			browser()
			# E step
			
			z.tmp = matrix(ncol = distinctn,nrow = clu)
			for (i in 1:clu){
				z.tmp[i,] = FindProb(dat,init@pi0.init[[i]],init@fai.init[[i]])
			}
			sums = colSums(z.tmp)
			for (i in 1:distinctn){
				z[,i] = z.tmp[,i]/sums[i]
			}
			
			# M step
			p.vec = z %*% count
			p.vec = p.vec/sum(p.vec)
			for ( i in 1:clu){
				#count.clu = z[i,] * count
				dat.clu = dat
				dat.clu@count = z[i,] * count
				init.clu = new("RankInit",fai.init=list(init@fai.init[[i]]),pi0.init=list(pi0.est[[i]]),clu=1L)
				solve.clu = AllSolve(dat.clu,init.clu,ctrl)
				
				pi0.est[[i]] = solve.clu$pi0.est
				# nrstep = length(solve.clu$solve.result$fai.est)
				fai[[i]] = solve.clu$fai.est
				fai.sig[[i]] = solve.clu$fai.sig
			}
			# break?
			if (loopind == 1){
				p.vec.last = p.vec
				pi0.est.last = pi0.est
				fai.last = fai
			} else if (loopind < ctrl@limit) {
				cond1 = all( p.vec.last - p.vec < ctrl@epsilon)
				cond2 = all( unlist(pi0.est.last) - unlist(pi0.est) < ctrl@epsilon)
				cond3 = all( unlist(fai.last) - unlist(fai) < ctrl@epsilon)
				if (cond1 && cond2 && cond3){
					break
				}
			} else {
				print(paste("Algorithm did not converge in",limit,"iterations"))
				return(NULL)
			}
		} # inf loop
		# fai.coeff = apply(ordering,1,WeightGivenPi,pi0)%*%count
		p.vec=as.numeric(p.vec)
		# log_likelihood = mixture_likelihood(ordering=ordering,count=count,p.vec=p.vec,pi0.est=pi0.est,fai=fai,clu=clu)
		return (list(p.vec=p.vec,pi0.est=pi0.est,fai=fai,fai.sig=fai.sig,iteration=loopind))
	}
)


MixtureSolve(testdat,testinit,testctrl)


