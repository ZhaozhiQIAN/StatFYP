############################  Express the Distance between p1 and p2 in Terms of fai ############################ 
# Important: true weights are denoted as W where w1 > w2 > ...
#            re-parametrized weights are denoted as fai where all fai > 0, and w(i) = sum(fai(i), ... , fai(t-1)
# p1,p2 are both orderings
WeightGivenPi = function(p1,p2){
	n = length(p1)
	w = numeric(length(p1)-1)
	distance = 0
	for (i in p2){
		pos1 = which(p1 == i)
		pos2 = which(p2 == i)
		relative_pos1 = (1:n - pos1)[order(p1)]
		relative_pos2 = (1:n - pos2)[order(p2)]
		Ji = which(relative_pos1 * relative_pos2 < 0)
		Ii = length(Ji)
		Li = (pos1 + pos2 + Ii)/2
		if(pos1 <= (Li-1)){
			w[pos1:(Li-1)] = w[pos1:(Li-1)] + 0.5
		}
		if (pos2<=(Li-1)){
			w[pos2:(Li-1)] = w[pos2:(Li-1)] + 0.5
		}
	}
	fai = cumsum(w)
	return(fai)
}

# find the weighted kendall distance between p1 and p2
# p1 and p2 are orderings
KwDist = function(p1, p2,w){
	n = length(p1)
	distance = 0
	for (i in p2){
		pos1 = which(p1 == i)
		pos2 = which(p2 == i)
		relative_pos1 = (1:n - pos1)[order(p1)]
		relative_pos2 = (1:n - pos2)[order(p2)]
		Ji = which(relative_pos1 * relative_pos2 < 0)
		Ii = length(Ji)
		Li = (pos1 + pos2 + Ii)/2
		c1 = ifelse(pos1<=(Li-1), sum(w[pos1:(Li-1)]),0)
		c2 = ifelse(pos2<=(Li-1), sum(w[pos2:(Li-1)]),0)
		distance = distance + (c1 + c2)/2
	}
	distance
}

############################ return list of t.lst ############################ 
# for internal use only
# to be used in GHC
# d is number of objects minus one
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

############################ transform between w and fai ############################ 
wTofai = function(w.true){
	fai.true = numeric(length(w.true))
	fai.true[1:(length(w.true)-1)] = -diff(w.true)
	fai.true[length(fai.true)] = w.true[length(w.true)]
	fai.true
}
faiTow = function(fai.true){
	w.true = rev(cumsum(rev(fai.true)))
	w.true
}

############################ transform between rankings and orderings ############################ 
# ranking: the position of objects (2 3 1 4)
# ordering: the ordered list objects (C A B D)

RankingToOrdering = function(ranking){
	if (class(ranking) %in% c("integer","numeric")){
		ranking = matrix(ranking,nrow=1)
	}
	ordering = matrix(ncol=ncol(ranking),nrow=nrow(ranking))
	for ( i in 1:nrow(ordering)){
		ordering[i,ranking[i,]] = 1:ncol(ranking)
	}
	ordering
}

OrderingToRanking = function(ordering){
	ranking = t(apply(ordering,1,order))
	ranking
}


############################ Calaculating Gradiant and Hessian for log(C) ############################ 
# arguments: fai in step t; t.lst returned by t.gen(d)
# result:	 a list containing gradient and Hessian of log(C)
GHC = function(fai,t.lst){
	d = length(fai) # d = t - 1
	K = matrix(rep(0,d^2),ncol = d, nrow = d)
	for ( i in 1:d){
		K = -1 * fai[i] * t.lst[[i]] + K
	}
	K = exp(K)
	K[upper.tri(K)] = 0
	gradiant = numeric(d)
	ones = rep(1,d)
	denom = rowSums(K) + ones
	B = matrix(ncol=d,nrow=d)
	for (i in 1:d){
		B[,i] = rowSums(-1 * K * t.lst[[i]])
	}
	for ( i in 1:d){
		gradiant[i] = sum(B[,i] / denom)
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
	return(list(gradiant=gradiant,hessian=hessian,K=K))
}

############################ use Newton's method to solve for mle of fai ############################ 
# arguments:
#	ordering/ranking: ordering or ranking matrix, each row represents a observation
#	pi0: specified mode ordering
#	fai.init: initial guess for fai
#	verbose: should intermediate Gradient vectors and inverse of information matrices be recorded
#	epsilon: if the difference between current fai and previous fai smaller than epsilon (for all elements in fai), halt
#	limit: maximum number of iterations
# returned values
#	fai.est: estimated fai
#	fai.sig: estimated covariance matrix of fai
#	log_likelihood: the log_likelihood using estimated fai
#	GHitr: intermediate Gradient vectors and inverse of information matrices (only if verbose is enabled)
#	time.len: the starting and ending time of the function (only if verbose is enabled)
NewtonSolve = function(ordering=NULL,count=NULL,pi0,fai.init,verbose=T,epsilon=1e-8,limit=5000,tor=0){
	if(verbose){
		start.time = Sys.time()
	}
	# arguments checking
	if(is.null(count)){
		count = rep(1,nrow(ordering))
	}
	# end of checking
	n = ncol(ordering)
	nobs = sum(count)
	t.lst = t.gen(n-1)
	fai.lst = list()
	fai.lst[[1]] = fai.init
	if (verbose == T){
		GHitr = list()
	}
	fai.coeff = apply(ordering,1,WeightGivenPi,pi0)%*%count
	
	while(TRUE){
		itr = length(fai.lst)
		fai.res = GHC(fai.lst[[itr]],t.lst)
		fai.gradiant = -1*nobs*fai.res$gradiant - fai.coeff
		fai.info = solve(nobs*fai.res$hessian)
		fai.next = fai.lst[[itr]] + fai.info%*%fai.gradiant
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
	# Naively coerce all negative values to 0
	tmpfai = fai.lst[[length(fai.lst)]]
	zeroind = tmpfai<0
	fai.lst[[length(fai.lst)]][zeroind] = 0
	# calculate coerced likelihood
	log_likelihood = -1*nobs*sum(log(rowSums(lastitr$K)+1)) - as.numeric(fai.coeff) %*% fai.lst[[length(fai.lst)]]
	if (verbose==T){
		end.time = Sys.time()
		return(list(fai.est=fai.lst,fai.sig=fai.sig,GHitr=GHitr,log_likelihood=log_likelihood,time.len=list(start.time=start.time,end.time=end.time)))
	}else {
		return(list(fai.est=fai.lst,fai.sig=fai.sig,log_likelihood=log_likelihood))
	}
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
	probhi = sort(count/nobs,decreasing=TRUE)[hi]
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

AllSolve_old = function(ordering=NULL,ranking=NULL,count,fai.init,verbose=T,epsilon=1e-8,limit=5000,simple=TRUE,pi0.init=NULL){
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
	if (is.null(pi0.init)){
		avg_rank = count %*% ranking;
		pi0.est = list()
		pi0.est[1] = list(sort(ordering[1,])[order(avg_rank)])
	} else {
		pi0.est = list()
		pi0.est[1] = list(pi0.init)
	}
# problem here
	solve.result = NewtonSolve(ordering=ordering,count=count,pi0=pi0.est[[1]],fai.init=fai.init,verbose=F)
	log_likelihood = list()
	log_likelihood[[1]] = solve.result$log_likelihood

	while(TRUE){
		# verify using adjacent swap
        if(length(pi0.est)%%10==0){
            print(length(pi0.est))
			print(pi0.est)
			print(log_likelihood)
        }
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
		
			verify.trial[i] = list(NewtonSolve(ordering=ordering,count=count,pi0=pi0.trial,fai.init=fai.init,verbose=F))
			verify.likelihood[i] = verify.trial[[i]]$log_likelihood
		}
		if(max(verify.likelihood) > target_log_likelihood){
			target_log_likelihood = max(verify.likelihood)
			log_likelihood[[length(log_likelihood)+1]] = target_log_likelihood
			pos = which(verify.likelihood == target_log_likelihood)
			pi0.est[[length(pi0.est)+1]] = unlist(pi0.trial.tmp[pos])
			solve.result = verify.trial[[pos]]
		} else {
			break
		}
	}
	if (simple==TRUE){
		fai.est = solve.result$fai.est[[length(solve.result$fai.est)]]
		fai.sig = solve.result$fai.sig
		pi0.est = pi0.est[[length(pi0.est)]]
		return(list(fai.est=fai.est,fai.sig=fai.sig,pi0.est=pi0.est,log_likelihood=log_likelihood[[length(log_likelihood)]]))
	}
	return(list(solve.result = solve.result,pi0.est=pi0.est,log_likelihood=log_likelihood))
}
# obtain C(fai)
faiC = function(fai,t.lst=NULL){
    if (is.null(t.lst))
	t.lst = t.gen(length(fai))
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

# count cannot be zero
moments.est = function(ordering=NULL,count,pi0.est,complete=FALSE,size){
	if (nrow(ordering) != length(count)){
		cat("ordering and count have different lenghts\n")
		return (NULL)
	}
	nfai = ncol(ordering)-1
	nobs = sum(count)
	ncount = length(count)
	# construct equation set
	prob = count/nobs
	if (complete == TRUE){
		samp = combn(length(count),2)
	} else {
		ind1.prob = rev(1:ncount)
		ind1 = sample(x=1:ncount,size=size,prob=ind1.prob,replace=T)
		ind2 = numeric(length(ind1))
		for (i in seq_along(ind1)){
			ind2[i] = sample(x=(ind1[i]+1):ncount,size=1)
		}
		samp=rbind(ind1,ind2)
		dup = duplicated(t(samp))
		samp = samp[,!dup]
	}
	nequ = ncol(samp)
	logodd = numeric(nequ)
	equmat = matrix(ncol=nfai,nrow=nequ)
	for (i in 1:nequ){
		logodd[i] = log(prob[samp[1,i]]/prob[samp[2,i]])
		kwj = WeightGivenPi(pi0.est,ordering[samp[2,i],])
		kwi = WeightGivenPi(pi0.est,ordering[samp[1,i],])
		equmat[i,] = kwj - kwi
	}
	fai.est = lm(logodd~equmat-1)
	fai.est
}

########################## mixture model ##########################
# pi0.mat each row is a pi0.init (ordering)
# clu number of clusters
# fai.mat rach row is a fai.init
MixtureSolve_old= function(ordering=NULL,ranking=NULL,count,fai.mat.init,pi0.mat.init,clu,verbose=T,epsilon=1e-8,limit=5000){
    # preparations 
	check1 = nrow(ordering)
	check2 = length(count)
	if (check1 != check2){
		cat("ordering and count have different lengths\n")
		return(NULL)
	}
    if(is.null(ranking)){
        ranking = OrderingToRanking(ordering)
    }
    tt=length(ordering[1,]) # number of objects
    n = sum(count)	# total observations
    distinctn = nrow(ranking)	# total distinct observations
    z.mat = matrix(ncol = distinctn,nrow = clu)
    # return values
    pi0.est = matrix(ncol=tt, nrow=clu)
    pi0.est.last = matrix(ncol=tt, nrow=clu)
    p.vec = numeric(clu)
    p.vec.last = numeric(clu)
    fai.mat = matrix(ncol=clu,nrow=tt-1)
    fai.mat.last = matrix(ncol=clu,nrow=tt-1)
	fai.sig.mat = list()

	loopind=0
    while(TRUE){
		loopind = loopind+1
		# handle loopind=1
		if (loopind == 1){
			pi0.est = pi0.mat.init
		}
		# E step
        z.mat.tmp = matrix(ncol = distinctn,nrow = clu)
        for (i in 1:clu){
            pi.now = pi0.mat.init[i,]
            t2 = apply(ordering,1,WeightGivenPi,pi.now) 
            z.mat.tmp[i,] = fai.mat.init[i,] %*% t2
        }
        z.mat.tmp = exp(-1*z.mat.tmp)
        sums = colSums(z.mat.tmp)
        for (i in 1:distinctn){
            z.mat[,i] = z.mat.tmp[,i]/sums[i]
        }

        # M step
        p.vec = z.mat %*% count
        p.vec = p.vec/sum(p.vec)
        for ( i in 1:clu){
			# fai.est=fai.est,fai.sig=fai.sig,pi0.est=pi0.est
            count.clu = z.mat[i,] * count
            solve.clu = AllSolve_old(ordering=ordering,count=count.clu,fai.init=fai.mat.init[i,],pi0.init=pi0.est[i,])
			  
            pi0.est[i,] = solve.clu$pi0.est
			# nrstep = length(solve.clu$solve.result$fai.est)
            fai.mat[,i] = solve.clu$fai.est
			fai.sig.mat[[i]] = solve.clu$fai.sig
        }
		# break?
		if (loopind == 1){
			p.vec.last = p.vec
			pi0.est.last = pi0.est
			fai.mat.last = fai.mat
		} else if (loopind < limit) {
			cond1 = all( p.vec.last - p.vec < epsilon)
			cond2 = all( pi0.est.last - pi0.est < epsilon)
			cond3 = all( fai.mat.last - fai.mat < epsilon)
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
	log_likelihood = mixture_likelihood(ordering=ordering,count=count,p.vec=p.vec,pi0.est=pi0.est,fai.mat=fai.mat,clu=clu)
	return (list(p.vec=p.vec,pi0.est=pi0.est,fai.mat=fai.mat,fai.sig.mat=fai.sig.mat,log_likelihood=log_likelihood,iteration=loopind))
}

mixture_likelihood = function(ordering,count,p.vec,pi0.est,fai.mat,clu){
	const = numeric(clu)
	t.lst = t.gen(ncol(ordering)-1)
	for(i in 1:clu){
		const[i] = faiC(fai.mat[,i],t.lst)
	}
	distmat = matrix(nrow=nrow(ordering),ncol=clu)
	for (i in 1:clu){
		mat1 = apply(ordering,1,WeightGivenPi,pi0.est[i,])
		distmat[,i] = t(mat1) %*% fai.mat[,i]
	}
	distmat = exp(-1*distmat) %*% (p.vec/const)
	likelihood = t(log(distmat)) %*% count
	likelihood
}