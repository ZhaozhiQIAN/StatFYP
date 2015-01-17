NewtonSolve = function(ordering=NULL,ranking=NULL,count=NULL,pi0,fai.init,verbose=T,epsilon=1e-8,limit=5000){
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
#	browser()
		itr = length(fai.lst)
		fai.res = GHC(fai.lst[[itr]],t.lst)
		fai.gradiant = -1*nobs*fai.res$gradiant - fai.coeff
		fai.info = solve(nobs*fai.res$hessian)
		fai.next = fai.lst[[itr]] + fai.info%*%fai.gradiant
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
		if (itr > limit){
			print(paste("does not converge in",limit,"iterations"))
			break
		}
	}
	lastitr = GHC(fai.lst[[length(fai.lst)]],t.lst)
	fai.sig = solve(nobs * lastitr$hessian)
	log_likelihood = -1*nobs*sum(log(rowSums(lastitr$K)+1)) - as.numeric(fai.coeff) %*% fai.lst[[length(fai.lst)]]
	if (verbose==T){
		end.time = Sys.time()
		return(list(fai.est=fai.lst,fai.sig=fai.sig,GHitr=GHitr,log_likelihood=log_likelihood,time.len=list(start.time=start.time,end.time=end.time)))
	}else {
		return(list(fai.est=fai.lst,fai.sig=fai.sig,log_likelihood=log_likelihood))
	}
}



AllSolve = function(ordering=NULL,ranking=NULL,count,fai.init,verbose=T,epsilon=1e-8,limit=5000,simple=TRUE){
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
		
			verify.trial[i] = list(NewtonSolve(ordering=ordering,ranking=ranking,count=count,pi0=pi0.trial,fai.init=fai.init,verbose=F))
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
		return(list(fai.est=fai.est,fai.sig=fai.sig,pi0.est=pi0.est))
	}
	return(list(solve.result = solve.result,pi0.est=pi0.est,log_likelihood=log_likelihood))
}


MixtureSolve= function(ordering=NULL,ranking=NULL,count,fai.mat.init,pi0.mat.init,clu,verbose=T,epsilon=1e-8,limit=5000){
    # preparations 
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
    # E step
	loopind=0
    while(loopind < 20){
		loopind = loopind+1
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
            solve.clu = AllSolve(ordering=ordering,count=count.clu,fai.init=fai.mat.init[i,])
            pi0.est[i,] = solve.clu$pi0.est
			# nrstep = length(solve.clu$solve.result$fai.est)
			# browser()
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
	return (list(p.vec=p.vec,pi0.est=pi0.est,fai.mat=fai.mat,fai.sig.mat=fai.sig.mat,iteration=loopind))
}


# count cannot be zero
moments.est = function(ordering=NULL,count,pi0.est,complete=FALSE,size){
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


RankingToOrdering = function(ranking){
	ordering = matrix(ncol=ncol(ranking),nrow=nrow(ranking))
	for ( i in 1:nrow(ordering)){
		ordering[i,ranking[i,]] = 1:ncol(ranking)
	}
	ordering
}

Pi0Est = function(ordering,count){
	ranking = OrderingToRanking(ordering)
	avg_rank = count %*% ranking;
	pi0.est = sort(ordering[1,])[order(avg_rank)]
	pi0.est
}





