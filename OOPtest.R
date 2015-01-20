# TESTING OBJECT-ORIENTED IMPLEMENTATION #
#         ^^^^^^^^^^^^^^^
source("functions.R")
source("OOP.R")
################ five object non-mixing test ###############
# true values
pi0 = c(1,2,3,4,5)
# w.true = c(1.2,1,0.5,0.02)
# large w: incomplete samples
# w.true = rep(4,4)
# increasing w: negative fais coerced to zero
w.true = c(2,1.2,1.8,1.5)
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
solve_result = NewtonStepSolve(o5na_suff,testinit,testctrl)
FindProb(o5na_dat,pi0,solve_result$fai.est[[6]])
#**** AllSolve ****#
AllSolve(o5na_dat,testinit,testctrl)
testinit_wrong_pi0 = new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(sample(5,5)),clu=1L)
AllSolve(o5na_dat,testinit_wrong_pi0,testctrl)

################ five object two clusters test ###############
pi01 = c(1,2,3,4,5)
pi02 = c(2,3,1,5,4)
w1 = rep(0.1,4)
w2 = c(0.8,0.5,0.2,0.01)
fai1 = wTofai(w1)
fai2 = wTofai(w2)
p_clu1 = 0.2
p_clu2 = 0.8
dist.true1 = apply(perm5,1,KwDist,pi01,w1)
probs.true1 = exp(-1*dist.true1)/sum(exp(-1*dist.true1))
dist.true2 = apply(perm5,1,KwDist,pi02,w2)
probs.true2 = exp(-1*dist.true2)/sum(exp(-1*dist.true2))
probs.true = p_clu1*probs.true1 + p_clu2*probs.true2
truesample = sample(1:gamma(6),prob=probs.true,size=10000,replace=T)
truesample = table(truesample)
testorder = matrix(nrow=nrow(truesample),ncol=5L)
selectedindex = as.integer(names(truesample))
for (i in 1:nrow(testorder)){
	testorder[i,] = perm5[selectedindex[i],]
}
o5c2_dat = new("RankData",nobj=5L,nobs=10000,ndistinct=nrow(testorder),ordering=testorder,ranking=RankingToOrdering(testorder),count=as.numeric(truesample))
# c2init = new("RankInit",fai.init=list(fai1,fai2),pi0.init = list(c(1,3,2,4,5),c(2,3,5,1,4)),clu=2L,p.init=c(0.2,0.8))
c2init = new("RankInit",fai.init=list(fai1,fai2),pi0.init = list(pi01,pi02),clu=2L,p.init=c(0.5,0.5))
c2ctrl = new("RankControl")
MixtureSolve(o5c2_dat,c2init,c2ctrl)

log(probs.true1*p_clu1 + probs.true2*p_clu2) %*% truesample

# try single cluster
c2init_single = new("RankInit",fai.init=list(rep(1,4)),pi0.init = list(pi01),clu=1L)
single_cluster_model = AllSolve(o5c2_dat,c2init_single,testctrl)
single_cluster_prob = FindProb(o5c2_dat,single_cluster_model$pi0.est,single_cluster_model$fai.est)
# EM converges to the single cluster model










