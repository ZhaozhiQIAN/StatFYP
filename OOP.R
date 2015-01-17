library("Rcpp")

sourceCpp("convolve.cpp")
# handle parameter passing with objects

# class definition
setClass( "RankData",
	representation = representation(
		nobj = "integer",
		nobs = "numeric",
		ndistinct = "numeric",
		ordering = "matrix",
		ranking = "matrix",
		count = "numeric"
	)
)

setClass( "SuffData",
	representation = representation(
		nobj = "integer",
		nobs = "numeric",
		ndistinct = "numeric",
		fai.coeff = "numeric",
		count = "numeric",
		t.lst = "list"
	)
)

setClass( "RankInit",
	representation = representation(
		fai.init = "list",
		pi0.init = "list",
		clu = "integer"
	)
)

setClass( "RankControl",
	representation = representation(
		verbose = "logical",
		epsilon = "numeric", # convergence of NR and EM algorithm
		limit = "integer", # maximum number of iterations
		stepstop = "numeric",	# minimum increase to continue stepwise selection
		tor = "numeric"
	),
	prototype = prototype(
		verbose = FALSE,
		epsilon = 1e-8,
		limit = 5000L,
		stepstop = 1e-3,
		tor = 0
	)
)


# constructors
# ranking and ordering MUST be integer matrix ranging from 1 to nobj (maximum 26 objects supported)
RankData = function(ordering=NULL,ranking=NULL,...){
	if( is.null(ordering) && is.null(ranking)){
		cat("argument is missing\nno object created\n")
		return(NULL)
	} else if (is.null(ordering)){
		nobj = ncol(ranking)
		count = apply(ranking, 1, function(x){paste(letters[x], collapse="")})
		count = table(count)
		ranking = unlist(strsplit(names(count),split=""))
		ranking = as.integer(factor(ranking,levels=letters,labels=1:26))
		ranking = matrix(ranking,ncol=nobj,byrow=TRUE)
		ordering = RankingToOrdering(ranking)
		
	} else if (is.null(ranking)){
		nobj = ncol(ordering)
		count = apply(ordering, 1, function(x){paste(letters[x], collapse="")})
		count = table(count)
		ordering = unlist(strsplit(names(count),split=""))
		ordering = as.integer(factor(ordering,levels=letters,labels=1:26))
		ordering = matrix(ordering,ncol=nobj,byrow=TRUE)
		ranking = RankingToOrdering(ordering)
	}
	count = as.numeric(count)
	nobs = sum(count)
	ndistinct = nrow(ranking)
	new("RankData",nobs = nobs,nobj=nobj,ndistinct=ndistinct,ordering=ordering,ranking=ranking,count=count,...)
}


setGeneric(name="ExtractSuff",
	def=function(dat,pi0){standardGeneric("ExtractSuff")}
)

setMethod(
	f="ExtractSuff",
	signature=c("RankData","numeric"),
	definition=function(dat,pi0){
		# fai.coeff = apply(dat@ordering,1,WeightGivenPi,pi0)%*%dat@count
		fai.coeff = CWeightGivenPi(dat@ranking,pi0)
		fai.coeff = matrix(fai.coeff,ncol = dat@ndistinct,byrow = TRUE)%*%dat@count
		fai.coeff = as.numeric(fai.coeff)
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
		t.lst = t.gen(dat@nobj-1)
		new("SuffData",nobs = dat@nobs,nobj=dat@nobj,ndistinct=dat@ndistinct,fai.coeff=fai.coeff,count=dat@count,t.lst=t.lst)
	}
)

setGeneric(name="FindProb",
	def=function(dat,pi0,fai){standardGeneric("FindProb")}
)

setMethod(
	f="FindProb",
	signature=c("RankData","numeric","numeric"),
	definition=function(dat,pi0,fai){
		# distance = fai %*% apply(dat@ordering,1,WeightGivenPi,pi0)
		distance = fai %*% matrix(CWeightGivenPi(dat@ranking,pi0),ncol = dat@ndistinct,byrow = TRUE)
		C = faiC(fai)
		prob = exp(-1*distance)/C
		prob
	}
)


setGeneric(name="NewtonSolveCon",
	def=function(dat,init,constraint,ctrl){standardGeneric("NewtonSolveCon")}
)
setMethod(
	f="NewtonSolveCon",
	signature=c(dat="SuffData",init="RankInit",constraint="logical",ctrl="RankControl"),
	definition=function(dat,init,constraint,ctrl){
		n = dat@nobj
		nobs = dat@nobs
		fai.init = init@fai.init[[init@clu]]
		if(ctrl@verbose){
			start.time = Sys.time()
		}
		# end of checking
		t.lst = dat@t.lst
		fai.lst = list()
		# enforce constraint on initial value
		fai.init[!constraint] = 0
		fai.lst[[1]] = fai.init
		if (ctrl@verbose == T){
			GHitr = list()
		}
		fai.coeff.con = dat@fai.coeff[constraint]
		while(TRUE){
			itr = length(fai.lst)
			fai.res = GHC(fai.lst[[itr]],t.lst)
			grad = fai.res$gradiant[constraint]
			fai.gradiant = -1*nobs*grad - fai.coeff.con
			fai.info = try(solve(nobs*fai.res$hessian[constraint,constraint]),silent = TRUE)
			if (class(fai.info) == "try-error"){
				print ("Constrained Newton Solve: non-invertible")
				return(list(log_likelihood=-Inf))
			}
			fai.next.con = rep(0,n-1)
			fai.next.con[constraint] = fai.info%*%fai.gradiant
			fai.next = fai.lst[[itr]] + fai.next.con
			# handle 0 in each step
			zeroind = which(fai.next < 0)
			fai.next[zeroind] = ctrl@tor
			# handle end
			fai.lst[[itr+1]]=fai.next
			if (ctrl@verbose==T){
				GHitr[[itr]] = list(fai.gradiant=fai.gradiant,fai.info=fai.info)
			}

			case1 = all( abs(fai.next - fai.lst[[itr]]) < ctrl@epsilon )
			if( is.na(case1)){
				print("doesn't work")
				return(list(log_likelihood=-Inf))
			}
			if ( case1 ){
				break
			}
			if (itr > ctrl@limit){
				print(paste("does not converge in",ctrl@limit,"iterations"))
				return(list(log_likelihood=-Inf))
				break
			}
		}
		lastitr = GHC(fai.lst[[length(fai.lst)]],t.lst)
		fai.sig = solve(nobs * lastitr$hessian)
		# calculate coerced likelihood
		log_likelihood = -1*nobs*sum(log(rowSums(lastitr$K)+1)) - as.numeric(dat@fai.coeff) %*% fai.lst[[length(fai.lst)]]
		if (ctrl@verbose==T){
			end.time = Sys.time()
			return(list(fai.est=fai.lst,fai.sig=fai.sig,GHitr=GHitr,log_likelihood=log_likelihood,time.len=list(start.time=start.time,end.time=end.time)))
		}else {
			return(list(fai.est=fai.lst,fai.sig=fai.sig,log_likelihood=log_likelihood))
		}
	}
)



setGeneric(name="NewtonStepSolve",
	def=function(dat,init,ctrl){standardGeneric("NewtonStepSolve")}
)


setMethod(
	f="NewtonStepSolve",
	signature=c(dat="SuffData",init="RankInit",ctrl="RankControl"),
	definition = function(dat,init,ctrl){
		nfai = dat@nobj-1
		included = list()
		bic.prev = -Inf
		champ = list()
		for (i in 1:(nfai-1)){
			trialmat = matrix(rep(FALSE,nfai*nfai),ncol=nfai)
			diag(trialmat) = rep(TRUE,nfai)
			if (length(included) > 0){
				trialmat = trialmat[,-unlist(included),drop=FALSE]
				trialmat[unlist(included),] = TRUE
			}
			solve.result = list()
			log_likelihood = numeric(ncol(trialmat))
			for (i in 1:ncol(trialmat)){
				solve.result[[i]] = NewtonSolveCon(dat,init,trialmat[,i],ctrl)
				log_likelihood[i] = solve.result[[i]]$log_likelihood
			}
			ind = which.max(log_likelihood)
			actconst = solve.result[[ind]]$fai.est
			actconst = actconst[[length(actconst)]]
			actconst = sum(actconst!=0)
			bic = 2*max(log_likelihood) - log(dat@nobs)*actconst
			if (bic - bic.prev > ctrl@stepstop){
				bic.prev = bic
				included = list(which(trialmat[,ind]==TRUE))
				champ = solve.result[[ind]]
			} else {
				return(champ)
			}
		}
		print(included)
		return(champ)
	}
)

# modification needed
setGeneric(name="MixtureSolve",
	def=function(dat,init,ctrl){standardGeneric("MixtureSolve")}
)


setMethod(
	f = "MixtureSolve",
	signature = c(dat="RankData",init="RankInit",ctrl="RankControl"),
	definition = function(dat,init,ctrl){
		tt=dat@nobj # number of objects
		distinctn = dat@ndistinct	# total distinct observations
		clu = init@clu	# total number of clusters
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
				count.clu = z.mat[i,] * count
				solve.clu = AllSolve(ordering=ordering,count=count.clu,fai.init=fai.mat.init[i,],pi0.init=pi0.est[i,])
				
				pi0.est[i,] = solve.clu$pi0.est
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
)



setGeneric(name="AllSolve",
	def=function(dat,init,ctrl){standardGeneric("AllSolve")}
)


# init should only contain one cluster
setMethod(
	f="AllSolve",
	signature=c(dat="RankData",init="RankInit",ctrl="RankControl"),
	definition=function(dat,init,ctrl){
		swap = function(arr,ind1,ind2){
			tmp = arr[ind1]
			arr[ind1] = arr[ind2]
			arr[ind2] = tmp
			arr
		}
		n = dat@nobj
		pi0.est = init@pi0.init # init@pi0.init should be a list with one element numeric
		suffdat = ExtractSuff(dat,pi0.est[[1]])
		solve.result = NewtonStepSolve(suffdat,init,ctrl)
		log_likelihood = list()
		log_likelihood[[1]] = solve.result$log_likelihood
		while(TRUE){
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
				suffdat.branch = ExtractSuff(dat,pi0.trial)
				init.branch = init
				init.branch@pi0.init = list(pi0.trial)
				verify.trial[i] = list(NewtonStepSolve(suffdat.branch,init.branch,ctrl))
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
		fai.est = solve.result$fai.est[[length(solve.result$fai.est)]]
		fai.sig = solve.result$fai.sig
		pi0.est = pi0.est[[length(pi0.est)]]
		return(list(fai.est=fai.est,fai.sig=fai.sig,pi0.est=pi0.est,log_likelihood=log_likelihood[[length(log_likelihood)]]))
	}
)

	
