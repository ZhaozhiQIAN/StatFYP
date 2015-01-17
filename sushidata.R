# load sushi data complete part
setwd("D:/Dropbox/study/Final Year Project/datasets/sushi")
source(file="../../simulation/functions.r")
sushia = read.table("sushi3a.5000.10.order",skip=1)
sushia.ranking = as.matrix(sushia[,3:12])
head(sushia.ranking)
dupp = duplicated(sushia.ranking)
sum(dupp)

sushia.count = table(apply(sushia.ranking, 1, paste, collapse=""))
sushia.ranking = matrix(as.numeric(unlist(strsplit(names(sushia.count),split=""))),ncol=10,byrow=TRUE)
sushia.ranking = sushia.ranking+1
sushia.count = as.numeric(sushia.count)
sushia.ordering = RankingToOrdering(sushia.ranking)
meanrank = colMeans(sushia.ranking)
pi0.mean = order(meanrank)
fai.mom=moments.est(ordering = sushidata@ordering,count = sushidata@count,pi0.mean)
fai.mom=moments.est(ordering = sushia.ordering,count = sushia.count,pi0.mean)
fai.init = fai.mom$coef
fai.init

pi0mat.guess = rbind(pi0.mean,rev(pi0.mean))
faimat.guess = rbind(fai.init,fai.init)
test2 = MixtureSolve(ordering=sushia.ordering,count=sushia.count,clu=2,fai.mat.init = faimat.guess,pi0.mat.init = pi0mat.guess)


test1 = AllSolve(ordering=sushia.ordering,count=sushia.count,fai.init=fai.init,verbose=F,pi0.init=pi0.mean)

test0 = NewtonSolve(tor = 0,ordering=sushia.ordering,count=sushia.count,pi0=pi0.mean,fai.init,verbose=T,limit=5000)

for(i in 1:9){
fai1 = sapply(test0[[1]],function(x){a=x[i]
a})
len=length(fai1)
plot(fai1[(len-200):len],ylab=paste("fai",i),main=paste("last 200 itr of fai",i))
}


dup = duplicated( matrix(unlist(test0[[1]]),nrow=9,byrow=TRUE))
sum(dup)