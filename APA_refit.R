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
apa_obj = new("RankData",nobj=5L,nobs=sum(capa.dat$Count),ndistinct=120,ranking=capa.ranking,ordering=RankingToOrdering(capa.ranking),count=capa.dat$Count)

testctrl = new("RankControl")
testinit = new("RankInit",fai.init=list(rep(1,4)),pi0.init=list(c(2,3,4,1,5)),clu=1L)
apa_c1 = AllSolve(apa_obj,testinit,testctrl)

c2init = new("RankInit",fai.init=list(rep(0.1,4),rep(0.1,4)),pi0.init = list(c(2,3,4,1,5),c(2,5,1,4,3)),clu=2L,p.init=c(0.5,0.5))
c2ctrl = new("RankControl")
two_cluster_model = MixtureSolve(apa_obj,c2init,c2ctrl)
lapply(two_cluster_model$pi0.est,RankingToOrdering)
prob.est1 = FindProb(apa_obj,two_cluster_model$pi0.est[[1]],two_cluster_model$fai[[1]])
prob.est2 = FindProb(apa_obj,two_cluster_model$pi0.est[[2]],two_cluster_model$fai[[2]])
prob.est = two_cluster_model$p[1]*prob.est1+two_cluster_model$p[2]*prob.est2
two_cluster_expectation = prob.est*apa_obj@nobs
plot(1:120,apa_obj@count,main="cluster_exp2",ylab="expectation/observation")
points(1:120,two_cluster_expectation,col="blue",pch=3)
legend(0,150,c("Observation","Estimated 2c"),pch=c(1,3),col=c("black","blue"))
PlotExp(apa_obj,two_cluster_model)

c3init = new("RankInit",fai.init=list(rep(0.1,4),rep(0.1,4),rep(0.1,4)),pi0.init = list(c(3,4,5,1,2),c(2,3,1,5,4),c(4,2,5,3,1)),clu=3L,p.init=rep(1,3)/3)
c3ctrl = new("RankControl")
three_cluster_model = MixtureSolve(apa_obj,c3init,c3ctrl)
PlotExp(apa_obj,three_cluster_model)