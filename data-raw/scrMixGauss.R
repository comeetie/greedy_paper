
library(greed)
library(future)

library(MASS)
N=600
X = rbind(mvrnorm(N/3,c(-5,0),diag(2)),mvrnorm(N/3,c(0,5),diag(2)),mvrnorm(N/3,c(5,0),diag(2)))
sol=greed(X)
plot(X,col=sol@cl)
No = 20
X = rbind(mvrnorm(N-No,c(0,0),diag(2)),mvrnorm(No,c(0,0),diag(2)*10))
sol=greed(X)
plot(X,col=sol@cl)
x1 = rnorm(N/2)
x2 = rnorm(N/2)
noise = 0.2
X = rbind(cbind(x1,x1+rnorm(N/2)*noise),cbind(x2,-x2+rnorm(N/2)*noise))
sol=greed(X)
plot(X,col=sol@cl)

crux = function(N,m,noise){
  x1 = rnorm(N/2)
  x2 = rnorm(N/2)
  X = rbind(cbind(x1,x1+rnorm(N/2)*noise),cbind(x2,-x2+rnorm(N/2)*noise))
  X+rep(m,each=N)
}

noise=0.2
X = rbind(crux(1000,c(-10,-10),noise),crux(1000,c(10,-10),noise),crux(1000,c(10,10),noise),crux(1000,c(-10,10),noise))
plan(multisession)
sol=greed(X,verbose = TRUE)
sol@icl
plot(X,col=sol@cl)
sol=greed:::fit_greed(new("mvmreg",beta=0.001,N0=5),list(X=matrix(1,nrow = 4000),Y=X),rep(1:8,each=500))
sol@icl
plot(X,col=sol@cl)
tail(sol@train_hist)

K=10
X=c()
cl=c()
for (k in 1:K){
  mu=rep(0,20)
  mu[k]=5
  X=rbind(X,mvrnorm(300,mu,diag(20)))
  cl=c(cl,rep(k,300))
}

library(MASS)


set.seed(100)
x <- matrix(rpois(70, 6), 10, 7) 
x.new <- makemultdata(x, cuts = 5)
multmixmodel.sel(x.new$y, comps = c(1,2), epsilon = 1e-03)



