# import APA data
setwd("D://Dropbox//study//Final\ Year\ Project//datasets")
apa.dat = read.csv("APA Election Dataset.csv",header=T)
apa.dat$Ranking = as.character(apa.dat$Ranking)
attach(apa.dat)
case1 = which(nchar(Ranking) != 5)
case2 = grep("0",Ranking,value=F)
discard = unique(c(case1,case2))
complete.pos = (1:length(Ranking))[-discard]
capa.dat = apa.dat[complete.pos,]
capa.ranking = matrix(as.numeric(unlist(strsplit(capa.dat$Ranking,""))),ncol=5,byrow=T)
capa.ordering = RankingToOrdering(capa.ranking)

# convertRank(capa.ranking)
# all(t(apply(capa.ranking,1,convertRank)) == capa.ordering)

# fit model
fai.mom=moments.est(ordering = capa.ordering,count = capa.dat$Count,c(3,4,5,1,2))
fai.init = fai.mom$coef
fai.init
debug(AllSolve)
source("functions.r")
test = AllSolve(pi0.init = papersorder[5,],ordering=capa.ordering,count=capa.dat$Count,fai.init=fai.mom$coefficients,verbose=F,simple=T)
OrderingToRanking(rbind(test$pi0.est,test$pi0.est))

pistartmat = combn(5,3)
for (i in 1:ncol(pistartmat)){
    st1 = pistartmat[1,i]
    st2 = pistartmat[2,i]
    st3 = pistartmat[3,i]
    pi0mat.guess = rbind(papersorder[st1,],papersorder[st2,],papersorder[st3,])
    faimat.guess = rbind(fai.init,fai.init,fai.init)
    test2 = MixtureSolve(ordering=capa.ordering,count=capa.dat$Count,clu=3,fai.mat.init = faimat.guess,pi0.mat.init = pi0mat.guess)
    test2$pi0.est
    test2$fai.mat
    test2$p.vec
    log2=mixture_likelihood(ordering=capa.ordering,count=capa.dat$Count,clu=3,pi0.est=test2$pi0.est,fai.mat=test2$fai.mat,p.vec=test2$p.vec)
    bic2 = 2*log2-log(sum(capa.dat$Count))*(sum(test2$fai.mat!=0)+length(test2$p.vec)*2)
    if(bic2>bic3){break}
}

pistartmat = combn(5,2)
for (i in 5:7){
    st1 = pistartmat[1,i]
    st2 = pistartmat[2,i]
    pi0mat.guess2 = rbind(papersorder[st1,],papersorder[st2,])
    faimat.guess2 = rbind(fai.init,fai.init)
    test3 = MixtureSolve(ordering=capa.ordering,count=capa.dat$Count,clu=2,fai.mat.init = faimat.guess2,pi0.mat.init = pi0mat.guess2)
    test3$pi0.est
    test3$fai.mat
    test3$p.vec
    log3=mixture_likelihood(ordering=capa.ordering,count=capa.dat$Count,clu=2,pi0.est=test3$pi0.est,fai.mat=test3$fai.mat,p.vec=test3$p.vec)
    bic3 = 2*log3-log(sum(capa.dat$Count))*(sum(test3$fai.mat!=0)+length(test3$p.vec)*2)
    if(bic2<bic3){break}
}
# i=5,6,7
pistartmat = combn(5,4)
for (i in 1:ncol(pistartmat)){
    st1 = pistartmat[1,i]
    st2 = pistartmat[2,i]
    st3 = pistartmat[3,i]
    st4 = pistartmat[4,i]
    pi0mat.guess4 = rbind(papersorder[st1,],papersorder[st2,],papersorder[st3,],papersorder[st4,])
    faimat.guess4 = rbind(fai.init,fai.init,fai.init,fai.init)
    test4 = MixtureSolve(ordering=capa.ordering,count=capa.dat$Count,clu=4,fai.mat.init = faimat.guess4,pi0.mat.init = pi0mat.guess4)
    test4$pi0.est
    test4$fai.mat
    test4$p.vec
    log4=mixture_likelihood(ordering=capa.ordering,count=capa.dat$Count,clu=4,pi0.est=test4$pi0.est,fai.mat=test4$fai.mat,p.vec=test4$p.vec)
    bic4 = 2*log4-log(sum(capa.dat$Count))*(sum(test4$fai.mat!=0)+length(test4$p.vec)*2)
    if(bic4>bic3){break}
}
# does not work    


for(i in 1:5){
    test1 = AllSolve(pi0.init = papersorder[i,],ordering=capa.ordering,count=capa.dat$Count,fai.init=fai.init,verbose=F,simple=T)
    log1 = test1$log_likelihood
    bic1 = 2*log1-log(sum(capa.dat$Count))*(sum(test1$fai.est!=0)+1)
    if(bic1>bic3){break}
}

papersorder = as.matrix(read.table("orders.txt",header=F))
papersorder = RankingToOrdering(papersorder)
papersrank = OrderingToRanking(papersorder)
# test3 is the best

w.est=estw(test3$fai.mat[,1],test3$fai.sig.mat[[1]])
plot(main="Estimated w",x=(1:4),y=w.est[[1]] + 2*w.est[[2]],ylab="w",xlab="index",ylim=c(-0.5,0.6),type="p")
points(x=factor(1:4),y=w.est[[1]] - 2*w.est[[2]])
points(x=1:4,y=w.est[[1]],col="blue",pch=4)

w.est=estw(test3$fai.mat[,2],test3$fai.sig.mat[[2]])
plot(main="Estimated w",x=(1:4),y=w.est[[1]] + 2*w.est[[2]],ylab="w",xlab="index",ylim=c(-0.5,0.6),type="p")
points(x=factor(1:4),y=w.est[[1]] - 2*w.est[[2]])
points(x=1:4,y=w.est[[1]],col="blue",pch=4)


# evaluation of the goodness of fit

FitProb = function(p.vec,pi0.est,fai.mat,ordering=NULL){
	clu = length(p.vec)
	nobj = nrow(fai.mat)+1
	if (is.null(ordering)){
		library(permute)
		sethow = how(observed = T)
		perm = allPerms(nobj,sethow)
	} else {
		perm=ordering
	}
	if(clu == 1){
		pi0.est = matrix(pi0.est,nrow=1)
		fai.mat = matrix(fai.mat,ncol=1)
	}
	clu.prob = matrix(ncol=clu,nrow = nrow(perm))
	for( i in seq_along(p.vec)){
		pi0 = pi0.est[i,]
		fai = fai.mat[,i]
		w.true = faiTow(fai)
		dist.true = apply(perm,1,KwDist,pi0,w.true)
		probs.true = exp(-1*dist.true)/sum(exp(-1*dist.true))
		clu.prob[,i] = probs.true
	}
	prob = clu.prob %*% p.vec
	prob
}
# debug(FitProb)
mixprob = FitProb(test3$p.vec,test3$pi0.est,test3$fai.mat,ordering = capa.ordering)
mixcount = mixprob*sum(capa.dat$Count)
mixres = sum((mixcount-capa.dat$Count)^2/mixcount)

library(StatMethRank)
library(Rankcluster)
paperdist = matrix(nrow=nrow(capa.ranking),ncol=nrow(papersorder))
for (i in 1:nrow(capa.ranking)){
	for (j in 1:nrow(papersorder)){
		rank1 = capa.ranking[i,]
		rank2 = papersorder[j,]
		paperdist[i,j] = distCayley(rank1, rank2)
	}
}
paperlambda = c(0.16,0.79,1.52,1.81,1.72)
paperprob.vec = c(0.42,0.31,0.12,0.08,0.07)
paperprob.est = matrix(ncol=ncol(paperdist),nrow=nrow(paperdist))
papernorm = list()
for (i in 1:ncol(paperprob.est)){
	paperprob.est[,i] = exp(-1*paperdist[,i]*paperlambda[i])
	papernorm[[i]] = sum(paperprob.est[,i])
	paperprob.est[,i] = paperprob.est[,i]/sum(paperprob.est[,i])
}
papernorm = unlist(papernorm)
paperprob= paperprob.est %*% paperprob.vec
papercount = paperprob * sum(capa.dat$Count)
paperres = sum((papercount-capa.dat$Count)^2/papercount)

# paper's ordering should be treated as ranking and vice versa
# find BIC
paperlog = log(as.numeric(paperprob)) %*% capa.dat$Count
paperBIC = 2*paperlog - log(sum(capa.dat$Count))*10
bic3
-27373.08*2-log(sum(capa.dat$Count))*14
hivote = which(capa.dat$Count>100)
hirank = capa.ranking[hivote,]
hidist = matrix(nrow = nrow(hirank),ncol=nrow(papersrank))
for (i in 1:nrow(papersrank)){
	hidist[,i] = apply(hirank,1,distCayley,papersorder[i,])
}
hiprob.est = matrix(ncol=ncol(hidist),nrow=nrow(hidist))
for (i in 1:ncol(hiprob.est)){
	hiprob.est[,i] = exp(-1*hidist[,i]*paperlambda[i])
	hiprob.est[,i] = hiprob.est[,i]/papernorm[i]
}
hiprob = hiprob.est %*% paperprob.vec
hiprob %*% sum(capa.dat$Count)



library(pmr)
# create dset
capa.dset = data.frame(capa.ordering)
capa.dset = cbind(capa.dset,capa.dat$Count)
# unweighted model
tau.fit = dbm(capa.dset,dtype="tau")
rho.fit = dbm(capa.dset,dtype="rho")
# weighted model
wtau.fit = wdbm(capa.dset,dtype="tau")
wtau.fit@min*-1
wrho.fit = wdbm(capa.dset,dtype="rho2")
luce.fit = pl(capa.dset)

test$pi0.est
test$log_likelihood
test[[1]]
OrderingToRanking(ordering = matrix(unlist(test$pi0.est),nrow=3,byrow=T))
w.est=estw(test[[1]][[1]])
plot(main="Estimated w",x=(1:4),y=w.est[[1]] + 2*w.est[[2]],ylab="w",xlab="index",ylim=c(-0.5,0.6),type="p")
points(x=factor(1:4),y=w.est[[1]] - 2*w.est[[2]])
points(x=1:4,y=w.est[[1]],col="blue",pch=4)

library(StatMethRank)
mwtau.fit = mwdbm(capa.dset,G=3,dset.agg = TRUE, dtype = "Kendall")
mwtau.fit2 = mwdbm(capa.dset,G=2,dset.agg = TRUE, dtype = "Kendall",iter=20)


# top three rankings
# our model
top3.clu1 = sum(capa.dat$Count[grepl("241..",capa.dat$Ranking)])
top3.clu2 = sum(capa.dat$Count[grepl("452..",capa.dat$Ranking)])
top3.obs = colSums(matrix(capa.dat$Count,nrow=2))
top3.mix = colSums(matrix(mixcount,nrow=2))
sum((top3.mix - top3.obs)^2/top3.obs)
# paper model
top3.paper = colSums(matrix(papercount,nrow=2))
top3.paper = colSums(matrix(papercount,nrow=2))
sum((top3.paper - top3.obs)^2/top3.obs)