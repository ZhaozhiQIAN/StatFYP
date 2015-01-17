# testing 
# four objects test
pi0 = c(1,2,3,4,5,6)
w.true = exp(-0.5*c(1:5))
fai.true = wTofai(w.true)
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

library(permute)
candidates = rbind(1:6,allPerms(1:6))
dist.true = apply(candidates,1,KwDist,pi0,w.true)
probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
T6 = t.gen(5)


obs.sample = sample(x=1:gamma(7),prob=probs.true,size=100000,replace=T)
solve.result = NewtonSolve(ordering=candidates,count=table(obs.sample),pi0=pi0,fai.init=c(0.8,0.5,0.2,0.1,0.1))

#========================= wrong pi0 =========================#
pi_nonsense = sample(pi0,size=length(pi0),replace=F)
solve.result = NewtonSolve(ordering=candidates,count=table(obs.sample),pi0=pi_nonsense,fai.init=c(0.8,0.5,0.2,0.1,0.1))
# result of fai can be negative if pi0 is not correct
# but it converges

#========================= wrong initial values of fai =========================#
fai_nonsense = runif(min=0,max=1,n=5)
solve.result = NewtonSolve(ordering=candidates,count=table(obs.sample),pi0=pi0,fai.init=fai_nonsense)
# the problem is that the algorithm is highly sensitive to the initial values of fai
# even if fai.init is relatively close to fai.true, it can
# (1) generate singular Hessian matrix
# (2) generate NA in Hessian or gradient vector

#========================= estimate fai.init =========================#
# least square method: based on log odds
fai.trial = moments.est(ordering = candidates,count=table(obs.sample),size=10000)$coefficients
solve.result = NewtonSolve(ordering=candidates,count=table(obs.sample),pi0=pi0,fai.init=fai.trial)


solve.result$fai.est[[length(solve.result$fai.est)]]
solve.result$log_likelihood


fai.true
sqrt(diag(solve.result$fai.sig))
