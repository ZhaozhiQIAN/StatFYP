########################## Metropolis Algorithm 1 ##########################

# construct a graph representation of all possible permutations
# try the prime number encoding
# primest <- function(n){
    # p <- 2:n
    # i <- 1
    # while (p[i] <= sqrt(n)) {
        # p <-  p[p %% p[i] != 0 | p==p[i]]
        # i <- i+1
    # }
    # p
# }

# gamma(14) = 13! will cause integer overflow
# R cannot handle the permutations of 13 objects
library(permute)
sethow = how(observed = T)
perm5 = allPerms(5,sethow)
adj.lst = matrix(nrow=gamma(6),ncol=4)
for(i in 1:nrow(adj.lst)){
	for (s in 1:ncol(adj.lst)){
		candi = perm5[i,];
		candi[s:(s+1)] = perm5[i,(s+1):s]
		for (j in 1:nrow(adj.lst)){
			if(all(perm5[j,] == candi)){
				adj.lst[i,s] = j
			}
		}
	}
}

# testing the algorithm
pi0 = c(1,2,3,4,5)
w.true = rep(0.5,4)
dist.true = apply(perm5,1,KwDist,pi0,w.true)
probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
p = exp(-1*w.true)
p = p/sum(p)
# generating 10000 samples
perm5.samp = numeric(500000)
perm5.samp[1] = 1
for (i in 1:(length(perm5.samp)-1)){
	if(i %% 10000 == 0){
		print(i)
	}
	current_pi = perm5.samp[i]
	next_pi = sample(x=adj.lst[current_pi,],prob=p,size=1)
	p.accnew = min(probs.true[next_pi]/probs.true[current_pi],1)
	perm5.samp[i+1] = sample(x=c(next_pi,current_pi),prob=c(p.accnew,1-p.accnew),size=1)
}

# generating 10000 samples by direct sampling
truesample = sample(1:gamma(6),prob=probs.true,size=10000,replace=T)
truesample = table(truesample)
perm5.table = rep(0,gamma(6))
perm5.table[1:length(table(perm5.samp))] = table(perm5.samp)
true.table =  rep(0,gamma(6))
true.table[1:length(table(truesample))] = table(truesample)
chisq.test(perm5.table,p=probs.true)
chisq.test(true.table,p=probs.true)



# adjacent transposition: acceptance rate is too high
Metropolis.simu = function(pi0, w, size){
	swap = function(arr,ind){
		tmp = arr[ind+1]
		arr[ind+1] = arr[ind]
		arr[ind] = tmp
		arr
	}
	accrate = 0
	T = length(pi0)
	samp = matrix(ncol=T,nrow=size)
	samp[1,] = pi0
	p = exp(-1*w)
	p = p/sum(p)
	pos = sample(1:(T-1),prob=p,size=(size-1),replace=TRUE)
	for(i in 2:size){
		if( i %% 1000 == 0){
			print(paste("observation",i,"is being sampled"))
		}
		current_pi = samp[i-1,]
		current.prob = exp(-1*KwDist(current_pi,pi0,w))
		next_pi = swap(current_pi,pos[i-1])
		next.prob = exp(-1*KwDist(next_pi,pi0,w))
		p.accnew = min(next.prob/current.prob,1)
		acc = sample(x=c(0,1),prob=c(p.accnew,1-p.accnew),size=1)
		if (acc==0){
			samp[i,]=next_pi
			accrate = accrate+1
		} else {
			samp[i,] = current_pi
		}
	}
	countable = table(apply(samp, 1, paste, collapse="")) 
	ordering = matrix(unlist(strsplit(names(countable),split="")),ncol=T,byrow=TRUE)
	count = as.numeric(countable)
	accrate = accrate/size
	return(list(ordering,count,accrate))
}

# version 2 use Cayley distance = 1
# favour longer distance exponentially 
Metropolis.simu.v2 = function(pi0, w, size, lambda=1){
	swap = function(arr,ind1,ind2){
		tmp = arr[ind1]
		arr[ind1] = arr[ind2]
		arr[ind2] = tmp
		arr
	}
	accrate = 0
	T = length(pi0)
	samp = matrix(ncol=T,nrow=size)
	samp[1,] = pi0
	# p = exp(-1*w)
	# p = p/sum(p)
	# pos = sample(1:(T-1),prob=p,size=(size-1),replace=TRUE)
	K = choose(T,2)
	ind1=1
	ind2=2
	pi0.propose = matrix(ncol=T,nrow=K)
	indmat = matrix(ncol=2,nrow=K)
	for ( i in 1:K){
		pi0star = swap(pi0,ind1,ind2)
		indmat[i,] = c(ind1,ind2)
		ind2 = ind2+1
		if(ind2 > T){
			ind1 = ind1 + 1
			ind2 = ind1 + 1
		}
		pi0.propose[i,] = pi0star
	}
	dist.propose = as.numeric(apply(pi0.propose,1,KwDist,pi0,w))
	p = exp(dist.propose*lambda)
	p = p/sum(p)
	pos = sample(1:K,prob=p,size=(size-1),replace=TRUE)
	for(i in 2:size){
		if( i %% 1000 == 0){
			print(paste("observation",i,"is being sampled"))
		}
		current_pi = samp[i-1,]
		current.prob = exp(-1*KwDist(current_pi,pi0,w))
		next_pi = swap(current_pi,indmat[pos[i-1],1],indmat[pos[i-1],2])
		next.prob = exp(-1*KwDist(next_pi,pi0,w))
		p.accnew = min(next.prob/current.prob,1)
		acc = sample(x=c(0,1),prob=c(p.accnew,1-p.accnew),size=1)
		if (acc==0){
			samp[i,]=next_pi
			accrate = accrate+1
		} else {
			samp[i,] = current_pi
		}
	}
	countable = table(apply(samp, 1, paste, collapse="")) 
	ordering = matrix(unlist(strsplit(names(countable),split="")),ncol=T,byrow=TRUE)
	count = as.numeric(countable)
	accrate = accrate/size
	return(list(ordering,count,accrate))
}
# 5
pi0=sample(1:5)
samp5 = Metropolis.simu.v2(pi0, w.true,10000)

wTofai(w.true)
undebug(moments.est)
perm5.mom = moments.est(ordering=perm5,count=truesample,pi0.est=pi0,complete=TRUE)
fai.ini=perm5.mom$coef

est = NewtonSolve(ordering=perm5,count=truesample,pi0=pi0,fai.init=fai.ini,verbose=FALSE)
sqrt(diag(est$fai.sig))

AllSolve(ordering=samp5[[1]],count=samp5[[2]],fai.init=fai.ini,verbose=FALSE)

# 10
pi0 = sample(1:9)
w.true = exp(-0.5*1:8)
samp10 = Metropolis.simu.v2(pi0, w.true,10000,1)

fai.iniest = moments.est(ordering=samp10[[1]],count=samp10[[2]],size=5000,pi0)
fai.ini = fai.iniest$coefficients
fai.ini


est = AllSolve(ordering=samp10[[1]],count=samp10[[2]],pi0=pi0,fai.init=fai.ini,verbose=FALSE)
fai.est = as.numeric(est$fai.est[[length(est$fai.est)]])
faiTow(fai.est)
fai.sd = sqrt(diag(est$fai.sig))
plot(x=1:4,y=fai.est + 2*fai.sd,ylim=c(-0.1,0.3))
points(x=1:4,y=fai.est - 2*fai.sd)
points(x=1:4,y=wTofai(w.true),col="red")
points(x=1:4,y=fai.est,col="blue")

#est is a object returned by NewtonSolve
estw = function(est){
	fai.est = as.numeric(est$fai.est[[length(est$fai.est)]])
	w.est = faiTow(fai.est)
	fai.cov = est$fai.sig
	n = length(w.est)
	w.var = numeric(n)
	for (i in seq_along(w.est)){
		w.var[i] = sum(fai.cov[i:n,i:n])
	}
	w.sd = sqrt(w.var)
	return(list(w.est,w.sd))
}
estw = function(fai.est,fai.cov){
	w.est = faiTow(fai.est)
	n = length(w.est)
	w.var = numeric(n)
	for (i in seq_along(w.est)){
		w.var[i] = sum(fai.cov[i:n,i:n])
	}
	w.sd = sqrt(w.var)
	return(list(w.est,w.sd))
}
w.est= estw(fai.est,est$fai.sig)

plot(main="Estimated w",x=(1:4),y=w.est[[1]] + 2*w.est[[2]],ylab="w",xlab="index",ylim=c(-0.1,0.8),type="p")
points(x=factor(1:4),y=w.est[[1]] - 2*w.est[[2]])
points(x=1:4,y=w.true,col="red",pch=3)
points(x=1:4,y=w.est[[1]],col="blue",pch=4)
legend("topright",inset=.08, c("True w","Estimated w"), col=c("red", "blue"), pch=3:4)
########################## test for goodness of fit ##########################
pi0 = 1:4
w = exp(-1*c(1:3)*2)
w = (1:3)^-1
size = 10000
test4 = Metropolis.simu(pi0, w, size)

prob.gen = function(pi0,w){
	library(permute)
	sethow = how(observed = T)
	n=length(pi0)
	perm5 = allPerms(n,sethow)
	dist.true = apply(perm5,1,KwDist,pi0,w)
	probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
	perm5 = apply(perm5, 1, paste, collapse="")
	return(list(perm5,probs.true))
}

prob4 = prob.gen(pi0,w)

all(apply(test4[[1]], 1, paste, collapse="") == prob4[[1]])

chisq.test(test4[[2]],p=prob4[[2]])

########################## clusters of probabilities? ##########################
pi0 = 1:5
w = exp(-1*c(1:5)*2)
prob5 = prob.gen(pi0,w)
plot(prob5[[2]],prob5[[2]])

pi0 = 1:6
w = pi0^(-1)
prob5 = prob.gen(pi0,w)
plot(prob5[[2]],prob5[[2]])

# Cluster due to exp decrease  of w

