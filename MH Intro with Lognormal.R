#
# lognormal scaling example for ASTIN 2017 Workshop- Glenn Meyers
#
rm(list = ls())      # clear workspace"
#
setwd("~/Dropbox/ASTIN 2017")
#
log.proposal.den=function(x,m,s){
  d=dnorm(x,m,s,log=T)
  return(d)
}
log.prior.den=function(mu){
  d=dnorm(mu,8,1,log=T)
  return(d)
}
#
# generate data
#
set.seed(12345)
nobs=25
sigma=1
y=round(rlnorm(nobs,8,sigma))
print(y)
#
#
thinned=NULL
nburn=1000
sigma.prop=0.01 #examples with 100,0.01 and 0.4
thin=1
n=nburn+thin*10000
for (i in 1:nburn/thin){
  thinn=c(T,rep(F,thin-1))
  thinned=c(thinned,thinn)
}
thinn=c(T,rep(F,thin-1))
for (i in 1:10000){
  thinned=c(thinned,thinn)
}
indext=subset(1:length(thinned),thinned==T)
#
# MH run1
#
set.seed(123)
mu1=rep(0,n)
mu1[1]=7 #starting value
yp1=rep(0,n)
for (i in 2:n){
  mu1[i]=mu1[i-1]
  mu1.prop=rnorm(1,mu1[i-1],sigma.prop)
  r=sum(dlnorm(mu1[i],sigma,log=T))+nobs*log.prior.den(mu1.prop)-
    sum(dlnorm(mu1[i-1],sigma,log=T))-nobs*log.prior.den(mu1[i-1])+
    log.proposal.den(mu1[i-1],mu1.prop,sigma.prop)-
    log.proposal.den(mu1.prop,mu1[i-1],sigma.prop)
  u=log(runif(1))
  if (u<r){
    mu1[i]=mu1.prop
  }
  yp1[i]=rlnorm(1,mu1[i],sigma)
}
#
# MH run2
#
set.seed(321)
mu2=rep(0,n)
mu2[1]=9 #starting value
yp2=rep(0,n)
for (i in 2:n){
  mu2[i]=mu2[i-1]
  mu2.prop=rnorm(1,mu2[i-1],sigma.prop)
  r=sum(dlnorm(mu2[i],sigma,log=T))+nobs*log.prior.den(mu2.prop)-
    sum(dlnorm(mu2[i-1],sigma,log=T))-nobs*log.prior.den(mu2[i-1])+
    log.proposal.den(mu2[i-1],mu2.prop,sigma.prop)-
    log.proposal.den(mu2.prop,mu2[i-1],sigma.prop)
  u=log(runif(1))
  if (u<r){
    mu2[i]=mu2.prop
  }
  yp2[i]=rlnorm(1,mu2[i],sigma)
}

#
# layer average severity statistics
#
library(actuar)
LEV=levlnorm(25000,mu1[thinned[1001:10000]],sigma)-
  levlnorm(10000,mu1[thinned[1001:10000]],sigma)
X=c(yp1[thinned[1001:11000]],yp2[thinned[1001:11000]])
X=pmax(X-10000,0)
X=pmin(X,15000)
Xpos=subset(X,X>0)
Xzero=sum(X==0)
par(mfrow=c(2,1))
hist(LEV,main="Parameter Risk",cex.main=0.9,breaks=25,cex.lab=0.9,cex.sub=0.9,
     xlab=expression(paste("Predictive Distribution of E[",italic(X),"]")),
     sub=paste("Predictive Mean =",round(mean(LEV))))
hist(Xpos,main="Total (Process + Parameter) Risk",cex.main=0.9,
     breaks=25,cex.lab=0.9,cex.sub=0.9,
     xlab=expression(paste("Predictive Distribution of ",italic(X))),
     sub=paste("Predictive Mean =",round(mean(X)),".There are",Xzero,"zero outcomes"))
#
# trace plots
#
y_range=range(mu1,mu2)
par(mfrow=c(1,1))
plot(indext,mu1[thinned],type="l",ylim=y_range,col="blue",lwd=2,
     main=paste("sigma.prop =",sigma.prop),
     xlab="Iteration", ylab=expression(mu))
par(new=T)
plot(indext,mu2[thinned],type="l",ylim=y_range,col="red",lwd=2,
     main="",xlab="",ylab="")
abline(v=nburn,lwd=10,col="gray")
