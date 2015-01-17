################ five object non-mixing test ###############
# true values
pi0 = c(1,2,3,4,5)
w.true = rep(4,4)
# generating samples
library(permute)
sethow = how(observed = T)
perm5 = allPerms(5,sethow)
dist.true = apply(perm5,1,KwDist,pi0,w.true)
probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
truesample = sample(1:gamma(6),prob=probs.true,size=10000,replace=T)
truesample = table(truesample)
testorder = matrix(nrow=nrow(truesample),ncol=5L)
selectedindex = as.integer(names(truesample))
for (i in 1:nrow(testorder)){
	testorder[i,] = perm5[selectedindex[i],]
}
#**** NewtonStepSolve ****#
# setting parameters
o5na_dat = new("RankData",nobj=5L,nobs=10000,ndistinct=nrow(testorder),ordering=testorder,ranking=RankingToOrdering(testorder),count=as.numeric(truesample))
o5na_suff = ExtractSuff(o5na_dat,pi0)
testctrl = new("RankControl")
testinit = new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(pi0),clu=1L)
NewtonStepSolve(o5na_suff,testinit,testctrl)
#**** AllSolve ****#
AllSolve(o5na_dat,testinit,testctrl)
testinit_wrong_pi0 = new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(sample(5,5)),clu=1L)
AllSolve(o5na_dat,testinit_wrong_pi0,testctrl)

################ five object two clusters test ###############
pi01 = c(1,2,3,4,5)
pi02 = c(2,3,1,5,4)
w1 = 


