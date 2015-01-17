library("Rcpp")

sourceCpp("convolve.cpp")
x=c(1,2,3)
y=c(2,3,4)
y = sum(y)
convolveCpp(x, y)

# project matrix

mat = diag(1:4)
ExtDiag(mat)
diag(mat)

p1=c(1,3,2,5,4)
p2=c(4,5,3,2,1)
r1=c(1,3,2,5,4)
r2=c(5,4,3,1,2)

sgn = sign(outer(r1,r1,'-')*outer(r2,r2,'-'))
apply(sgn,1,function(x) sum(x<0))
r1m=rbind(r1,r1)
ret = CWeightGivenPi(r1m,r2)
matrix(ret,nrow=2)

ret = CWeightGivenPiV2(r1m,r2)
matrix(ret,nrow=2)


rankmat = matrix(ncol=10, nrow=10000)
for (i in 1:nrow(rankmat)){
	tmp=sample(1:10,10,replace=FALSE)
	rankmat[i,] = tmp;
}

ordermat = RankingToOrdering(rankmat)

pi0=sample(1:10,10,replace=FALSE)
pi0rank=RankingToOrdering(pi0)

func1=function(){
	m1=apply(ordermat,1,WeightGivenPi,pi0)
}
t1=system.time(func1())
m1=func1()

func2=function(){
	m2=CWeightGivenPi(rankmat,pi0rank)
	m2=matrix(m2,ncol=10000,byrow=TRUE)
	m2
}
t2=system.time(func2())
m2=func2()

func3=function(){
	m3=CWeightGivenPiV2(rankmat,pi0rank)
	m3=matrix(m3,ncol=10000,byrow=TRUE)
	m3
}
t3=system.time(func3())
m3=func3()

t1[[3]]/t2[[3]]
t2[[3]]/t3[[3]]

all(m1==m2)
all(m2==m3)


