
library(greed)
library(future)
plan(multisession)
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


mo=new("gmm",N0=10)
Xc=cbind(X[,1],rep(1,4*Nt))
Y=matrix(X[,2],nrow=4*Nt)
K=20
data=list(X=Xc,Y=Y,N=4*Nt,move_mat=greed:::as.sparse(matrix(1,ncol=K,nrow=K)))
so=greed:::fit_greed(new("mvmreg",N0=2,beta=0.01),data,sample(1:K,4*Nt,replace = TRUE),verbose = TRUE,type='swap')
plot(X,col=so@cl)
sol=greed_cond(Xc,Y)
plot(X,col=sol@cl)


K=8
data=list(Y=greed:::zscore(X),X=matrix(1,nrow = Nt*4,ncol=1),N=Nt*4,move_mat=greed:::as.sparse(matrix(1,ncol=K,nrow=K)))
so=greed:::fit_greed(model,data,rep(1:8,each=Nt/2),verbose = TRUE,type='none')
sum(unlist(purrr::map(so@obs_stats$regs,function(r){r$log_evidence})))
so@icl
K=4
data=list(Y=greed:::zscore(X),X=matrix(1,nrow = Nt*4,ncol=1),N=Nt*4,move_mat=greed:::as.sparse(matrix(1,ncol=K,nrow=K)))
so=greed:::fit_greed(new("gmm",N0=2,beta=0.001),data,sol@cl,verbose = TRUE,type='none')
sum(unlist(purrr::map(so@obs_stats$regs,function(r){r$log_evidence})))
so@icl




plan(multisession)
sol=greed(X,verbose = TRUE,alg=new("hybrid",pop_size=40,prob_mutation=0.5),K=15)
sol@icl
plot(X,col=sol@cl)


sol=greed(X)
K=8
model=new("gmm",N0=2,tau=0.0001,mu=c(0,0),epsilon=100*diag(2))
data=greed:::preprocess(model,X,K)
cl=rep(1:8,each=200)
sol=greed:::fit_greed(model,data,cl,type='none')
sol@icl

tail(sol@train_hist)


cl=sol@cl
cl[cl==3]=sample(c(3,5),sum(cl==3),replace = TRUE)
plot(X,col=cl)
so=greed:::fit_greed(new("gmm",N0=2,beta=0.00001),data,cl,verbose = TRUE,type = 'swap')


so=greed:::fit_greed(new("gmm",N0=50),data,rep(1:8,each=500),verbose = TRUE,type='none')
sum(unlist(purrr::map(so@obs_stats$regs,function(r){r$log_evidence})))
so@icl
so=greed:::fit_greed(new("gmm",N0=50),data,rep(1:4,each=1000),verbose = TRUE,type='none')
sum(unlist(purrr::map(so@obs_stats$regs,function(r){r$log_evidence})))
so@icl


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



new("gmm",N0=4)

set.seed(100)
x <- matrix(rpois(70, 6), 10, 7) 
x.new <- makemultdata(x, cuts = 5)
multmixmodel.sel(x.new$y, comps = c(1,2), epsilon = 1e-03)



