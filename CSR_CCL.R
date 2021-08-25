#
# Script to run the CSR_CCL Model on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
####### Most of the script is in functions - They are called at the end.
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
setwd("~/Dropbox/ASTIN 2017")
#
# get packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(parallel)
library(doParallel)
#
# stan script
#
scodeU = "
data{
int<lower=1> len_data;
int<lower=0,upper=1> wne1[len_data];
real logprem[len_data];
real logcpd[len_data];
real loginc[len_data];
int<lower=1,upper=10> w[len_data];
int<lower=1,upper=10> d[len_data];
}
parameters{
real r_alpha[9];
real r_beta_p[9];
real r_beta_i[9];
real logelr_p;
real logelr_i;
real <lower=0,upper=100000> a_p[10];
real <lower=0,upper=100000> a_i[10];
real <lower=0,upper=1> r_rho;
real gamma;
}
transformed parameters{
real alpha[10];
real beta_p[10];
real beta_i[10];
real speedup[10];
real rho;
real sig2_p[10];
real sig2_i[10];
real sig_p[10];
real sig_i[10];
real mu_p[len_data];
real mu_i[len_data];
alpha[1] = 0;
for (i in 2:10) alpha[i] = r_alpha[i-1];
for (i in 1:9) beta_p[i] = r_beta_p[i];
for (i in 1:9) beta_i[i] = r_beta_i[i];
beta_p[10] = 0;
beta_i[10] = 0;
speedup[1] = 1;
for (i in 2:10) speedup[i] = speedup[i-1]*(1-gamma);
sig2_p[10] = gamma_cdf(1/a_p[10],1,1);
for (i in 1:9) sig2_p[10-i] = sig2_p[11-i]+gamma_cdf(1/a_p[i],1,1);
for (i in 1:10) sig_p[i] = sqrt(sig2_p[i]);
sig2_i[10] = gamma_cdf(1/a_i[10],1,1);
for (i in 1:9) sig2_i[10-i] = sig2_i[11-i]+gamma_cdf(1/a_i[i],1,1);
for (i in 1:10) sig_i[i] = sqrt(sig2_i[i]);
mu_p[1] = logprem[1]+logelr_p+beta_p[d[1]]*speedup[w[1]];
mu_i[1] = logprem[1]+logelr_i+beta_i[d[1]];
rho = -2*r_rho+1;
for (i in 2:len_data){
mu_p[i] = logprem[i]+logelr_p+alpha[w[i]]+beta_p[d[i]]*speedup[w[i]];
mu_i[i] = logprem[i]+logelr_i+alpha[w[i]]+beta_i[d[i]]+
rho*(loginc[i-1]-mu_i[i-1])*wne1[i];
}
}
model {
r_alpha ~ normal(0,3.162);
r_beta_p ~ normal(0,3.162);
r_beta_i ~ normal(0,3.162);
a_p ~ inv_gamma(1,1);
a_i ~ inv_gamma(1,1);
logelr_p ~ normal(-.4,3.162);
logelr_i ~ normal(logelr_p,0.05);
gamma ~ normal(0,0.05);
r_rho ~ beta(2,2);
for (i in 1:len_data) {
logcpd[i] ~ normal(mu_p[i],sig_p[d[i]]);
loginc[i] ~ normal(mu_i[i],sig_i[d[i]]);
}
}
generated quantities{
vector[len_data] log_lik;
for (i in 1:len_data){
log_lik[i] = normal_lpdf(logcpd[i]|mu_p[i],sig_p[d[i]])+
normal_lpdf(loginc[i]|mu_i[i],sig_i[d[i]]);
}
}
"
#
# function to get Schedule P triangle data given ins group
#
ins.line.data=function(g.code){
  b=subset(a,a$GRCODE==g.code)
  name=b$GRNAME
  grpcode=b$GRCODE
  w=b$AccidentYear
  d=b$DevelopmentLag
  cum_incloss=b[,6]
  cum_pdloss=b[,7]
  bulk_loss=b[,8]
  dir_premium=b[,9]
  ced_premium=b[,10]
  net_premium=b[,11]
  single=b[,12]
  posted_reserve97=b[,13]
  # get incremental paid losses - assume data is sorted by ay and lag
  inc_pdloss=numeric(0)
  for (i in unique(w)){
    s=(w==i)
    pl=c(0,cum_pdloss[s])
    ndev=length(pl)-1
    il=rep(0,ndev)
    for (j in 1:ndev){
      il[j]=pl[j+1]-pl[j]
    }
    inc_pdloss=c(inc_pdloss,il)
  }
  data.out=data.frame(grpcode,w,d,net_premium,dir_premium,ced_premium,
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,posted_reserve97)
  return(data.out)
}

#
# initialization function for scodeU
#
initU=function(chain_id){
  set.seed(123+chain_id)
  list(r_alpha=rnorm(9,0,0.2),r_beta_p=runif(9),
       r_beta_i=runif(9),a_i=runif(10),a_p=runif(10),
       logelr_i=runif(1,-0.75,-0.25),logelr_p=runif(1,-0.75,-0.25),
       gamma=rnorm(1,0,0.1),r_rho=runif(1))
}
#
#  data for the model
#
pars.list=c("alpha","beta_p","beta_i","gamma","logelr_p",
            "logelr_i","rho","sig_p","sig_i","log_lik")
#
# dummy data for compiling
#
data.dummy=list(len_data = 55,
                logprem  = rep(8,55),
                logcpd   = rep(8,55),
                loginc   = rep(8,55),
                w        = c(1:10,1:9,1:8,1:7,1:6,1:5,1:4,1:3,1,2,1),
                d        = c(rep(1,10),rep(2,9),rep(3,8),rep(4,7),rep(5,6),
                             rep(6,5),rep(7,4),rep(8,3),rep(9,2),10),
                wne1     =rep(1,55))

#
# compile the univariate model
#
fitU = stan(model_code=scodeU,data=data.dummy,seed=123,init=initU,chains=0)
#
# set up function to run stan model and create output
#
model_function=function(grpcode){
  #
  # read and aggregate the insurer data and
  # set up training and test data frames
  # for line 1
  #
  a=read.csv(insurer.data)
  cdata=ins.line.data(grpcode)
  w=cdata$w-1987
  d=cdata$d
  #
  # sort the data in order of d, then w within d
  #
  o1=100*d+w
  o=order(o1)
  w=w[o]
  d=d[o]
  premium=cdata$net_premium[o]
  cpdloss=cdata$cum_pdloss[o]
  cpdloss=pmax(cpdloss,1)
  incloss=cdata$cum_incloss[o]
  incloss=pmax(incloss,1)
  adata=data.frame(grpcode,w,d,premium,cpdloss,incloss)
  rdata=subset(adata,(adata$w+adata$d)<12)
  hdata=subset(adata,(adata$w+adata$d)>11)
  wne1=ifelse(rdata$w==1,0,1)
  premium=rdata$premium
  #
  #  data for the model
  #
  data.u=list(len_data = length(rdata$premium),
              logprem  = log(rdata$premium),
              logcpd   = log(rdata$cpdloss),
              loginc   = log(rdata$incloss),
              w        = rdata$w,
              d        = rdata$d,
              wne1     = wne1)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<65)){
    fitU1=stan(fit = fitU, data = data.u,init=initU,
               seed = 123,iter=stan_iter,thin=stan_thin,
               chains = 4,pars=pars.list,
               control=list(adapt_delta=.9999),
               refresh=0)
    fitU1_summary=as.matrix(summary(fitU1)$summary)[1:54,c(1,3,10)]
    mrh=subset(fitU1_summary,is.na(fitU1_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin=stan_thin/2
  #
  # goodness of fit statistics for comparing models
  #
  #
  # extract information from stan output to process in R
  #
  b=extract(fitU1)
  alpha=b$alpha
  beta_p=b$beta_p
  beta_i=b$beta_i
  gamma=b$gamma
  logelr_p=b$logelr_p
  logelr_i=b$logelr_i
  rho=b$rho
  sig_p=b$sig_p
  sig_i=b$sig_i
  #
  # simulate outcomes for d=10 using parallel processing
  #
  set.seed(12345)
  at_p.wd10=matrix(0,length(logelr_p),10)
  for (w in 1:10){
    at_p.wd10[,w]=rlnorm(length(logelr_p),log(premium[w])+logelr_p+alpha[,w],
                         sig_p[,10])
  }
  at_i.wd10=matrix(0,length(logelr_i),10)
  mu_i=matrix(0,length(logelr_i),10)
  mu_i[,1]=log(premium[1])+logelr_i+beta_i[,10]
  at_i.wd10[,1]=rep(rdata$incloss[55],length(logelr_i))
  for (w in 2:10){
    mu_i[,w]=log(premium[w])+logelr_i+alpha[,w]+rho*(log(at_i.wd10[,w-1])-mu_i[,w-1])
    at_i.wd10[,w]=rlnorm(length(logelr_i),mu_i[,w],sig_i[,10])
  }
  Premium=subset(rdata,rdata$d==1)$premium
  ss_p.wd10=rep(0,10)
  ms_p.wd10=rep(0,10)
  ss_i.wd10=rep(0,10)
  ms_i.wd10=rep(0,10)
  #
  ms_p.wd10[1]=mean(at_p.wd10[,1])
  ms_i.wd10[1]=mean(at_p.wd10[,1])
  for (w in 2:10){
    ms_p.wd10[w]=mean(at_p.wd10[,w])
    ss_p.wd10[w]=sd(at_p.wd10[,w])
    ms_i.wd10[w]=mean(at_i.wd10[,w])
    ss_i.wd10[w]=sd(at_i.wd10[,w])
  }
  Pred.IP_CSR=rowSums(at_p.wd10)
  ms_p.td10=mean(Pred.IP_CSR)
  ss_p.td10=sd(Pred.IP_CSR)
  IP_CSR.Estimate=round(ms_p.wd10)
  IP_CSR.SE=round(ss_p.wd10)
  IP_CSR.CV=round(IP_CSR.SE/IP_CSR.Estimate,4)
  act=sum(subset(adata$cpdloss,adata$d==10)[1:10])
  pct.IP_CSR=sum(Pred.IP_CSR<=act)/length(Pred.IP_CSR)*100
  #
  # put IP_CSR accident year statistics into a data frame
  #
  Pred.IP_CSR=rowSums(at_p.wd10[,1:10])
  ms_p.td10=mean(Pred.IP_CSR)
  ss_p.td10=sd(Pred.IP_CSR)
  IP_CSR.Estimate=round(ms_p.wd10)
  IP_CSR.SE=round(ss_p.wd10)
  IP_CSR.CV=round(IP_CSR.SE/IP_CSR.Estimate,4)
  act_p=sum(subset(adata$cpdloss,adata$d==10)[1:10])
  pct.IP_CSR=sum(Pred.IP_CSR<=act_p)/length(Pred.IP_CSR)*100
  IP_CSR.Estimate=c(IP_CSR.Estimate,round(ms_p.td10))
  IP_CSR.SE=c(IP_CSR.SE,round(ss_p.td10))
  IP_CSR.CV=c(IP_CSR.CV,round(ss_p.td10/ms_p.td10,4))
  Premium=c(Premium,sum(Premium))
  Outcome_P=subset(adata$cpdloss,adata$d==10)
  Outcome_P=c(Outcome_P,sum(Outcome_P))
  Group=rep(grpcode,11)
  IP_CSR.Pct=c(rep(NA,10),pct.IP_CSR)
  W=c(1:10,"Total")
  risk_P=data.frame(Group,W,Premium,IP_CSR.Estimate,IP_CSR.SE,
                    IP_CSR.CV,Outcome_P,IP_CSR.Pct)
  #write.csv(risk,file=outfilename,row.names=F)
  #
  # put IP_CCL accident year statistics into a data frame
  #
  Pred.IP_CCL=rowSums(at_i.wd10)
  ms_i.td10=mean(Pred.IP_CCL)
  ss_i.td10=sd(Pred.IP_CCL)
  IP_CCL.Estimate=round(ms_i.wd10)
  IP_CCL.SE=round(ss_i.wd10)
  IP_CCL.CV=round(IP_CCL.SE/IP_CCL.Estimate,4)
  act_I=sum(subset(adata$incloss,adata$d==10)[1:10])
  pct.IP_CCL=sum(Pred.IP_CSR<=act_I)/length(Pred.IP_CCL)*100
  IP_CCL.Estimate=c(IP_CCL.Estimate,round(ms_i.td10))
  IP_CCL.SE=c(IP_CCL.SE,round(ss_i.td10))
  IP_CCL.CV=c(IP_CCL.CV,round(ss_i.td10/ms_i.td10,4))
  #Premium=c(Premium,sum(Premium))
  Outcome_I=subset(adata$incloss,adata$d==10)
  Outcome_I=c(Outcome_I,sum(Outcome_I))
  Group=rep(grpcode,11)
  IP_CCL.Pct=c(rep(NA,10),pct.IP_CCL)
  W=c(1:10,"Total")
  risk_P=data.frame(Group,W,Premium,IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,
                    Outcome_P,IP_CSR.Pct)
  risk_I=data.frame(Group,W,Premium,IP_CCL.Estimate,IP_CCL.SE,IP_CCL.CV,
                    Outcome_I,IP_CCL.Pct)
  risk_IP=data.frame(Group,W,Premium,
                     IP_CCL.Estimate,IP_CCL.SE,IP_CCL.CV,Outcome_I,
                     IP_CCL.Pct,
                     IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,Outcome_P,
                     IP_CSR.Pct)
  #
  Group=grpcode
  Premium=Premium[11]
  IP_CSR.Estimate=IP_CSR.Estimate[11]
  IP_CSR.SE=IP_CSR.SE[11]
  IP_CSR.CV=IP_CSR.CV[11]
  Outcome_I=Outcome_I[11]
  IP_CSR.Pct=IP_CSR.Pct[11]
  #
  IP_CCL.Estimate=IP_CCL.Estimate[11]
  IP_CCL.SE=IP_CCL.SE[11]
  IP_CCL.CV=IP_CCL.CV[11]
  Outcome_P=Outcome_P[11]
  IP_CCL.Pct=IP_CCL.Pct[11]
  mean_rho=round(mean(rho),4)
  sd_rho=round(sd(rho),4)
  mean_gamma=round(mean(gamma),4)
  sd_gamma=round(sd(gamma),4)
  loglik1=extract_log_lik(fitU1)
  loo1 <- loo(loglik1)
  LOOIC=round(loo1$looic,3)
  LOOIC_SE=round(loo1$se_looic,3)
  SumStats=data.frame(Group,Premium,
                      IP_CCL.Estimate,IP_CCL.SE,IP_CCL.CV,Outcome_I,
                      IP_CCL.Pct,mean_rho,sd_rho,
                      IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,Outcome_P,
                      IP_CSR.Pct,mean_gamma,sd_gamma,
                      LOOIC,LOOIC_SE,max_Rhat,stan_thin)
  output=list(ParmSummary=fitU1_summary,
              risk_P     =risk_P,
              risk_I     =risk_I,
              risk_IP    =risk_IP,
              Pred_IP_CCL=Pred.IP_CCL,
              Pred_IP_CSR=Pred.IP_CSR,
              SumStats =SumStats,
              b        =b)
  return(output)
}
#
# Single triangle
#
insurer.data="~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv"
g=353
# outfilename=paste("CSR_CCL CA",g,".csv")
a=read.csv(insurer.data)
co_model=model_function(g)
print(co_model$ParmSummary)
print(co_model$risk_P)
print(co_model$risk_I)
par(mfrow=c(2,1))
rng=range(co_model$Pred_IP_CCL,co_model$Pred_IP_CSR)
hist(co_model$Pred_IP_CCL,xlab="Simulated Outcomes",xlim=rng,
     main="Predictive Distribution of Outcomes for Incurred Losses",
     sub=paste("Mean =",co_model$SumStats$IP_CCL.Estimate,
               " SE =",co_model$SumStats$IP_CCL.SE))
hist(co_model$Pred_IP_CSR,xlab="Simulated Outcomes",xlim=rng,
     main="Predictive Distribution of Outcomes for Paid Losses",
     sub=paste("Mean =",co_model$SumStats$IP_CSR.Estimate,
               " SE =",co_model$SumStats$IP_CSR.SE))
# write.csv(co_model$risk_IP,file=outfilename)
#
# loop through the set of insurere
#
# stats=NULL
# outfile="CSR_CCL Summary Statistics.csv"
# #
# # CA LOB
# #
# insurer.data="~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv"
# insurer.list="~/Dropbox/CAS Loss Reserve Database/Selected Insurers CA.csv"
# a=read.csv(insurer.data)
# l=read.csv(insurer.list)
# lob=substring(insurer.list,55,56)
# for (g in l$grpcode){
#   print(paste("Processing Insurer",g,"of",lob))
#   co_model=model_function(g)
#   stats=rbind(stats,data.frame(lob,co_model$SumStats))
#   write.csv(stats,file=outfile)
#   t1=Sys.time()
#   print(t1-t0)
# }
# #
# # PA LOB
# #
# insurer.data="~/Dropbox/CAS Loss Reserve Database/ppauto_pos.csv"
# insurer.list="~/Dropbox/CAS Loss Reserve Database/Selected Insurers PA.csv"
# a=read.csv(insurer.data)
# l=read.csv(insurer.list)
# lob=substring(insurer.list,55,56)
# for (g in l$grpcode){
#   print(paste("Processing Insurer",g,"of",lob))
#   co_model=model_function(g)
#   stats=rbind(stats,data.frame(lob,co_model$SumStats))
#   write.csv(stats,file=outfile)
#   t1=Sys.time()
#   print(t1-t0)
# }
# #
# # WC LOB
# #
# insurer.data="~/Dropbox/CAS Loss Reserve Database/wkcomp_pos.csv"
# insurer.list="~/Dropbox/CAS Loss Reserve Database/Selected Insurers WC.csv"
# a=read.csv(insurer.data)
# l=read.csv(insurer.list)
# lob=substring(insurer.list,55,56)
# for (g in l$grpcode){
#   print(paste("Processing Insurer",g,"of",lob))
#   co_model=model_function(g)
#   stats=rbind(stats,data.frame(lob,co_model$SumStats))
#   write.csv(stats,file=outfile)
#   t1=Sys.time()
#   print(t1-t0)
# }
# #
# # OL LOB
# #
# insurer.data="~/Dropbox/CAS Loss Reserve Database/othliab_pos.csv"
# insurer.list="~/Dropbox/CAS Loss Reserve Database/Selected Insurers OL.csv"
# a=read.csv(insurer.data)
# l=read.csv(insurer.list)
# lob=substring(insurer.list,55,56)
# for (g in l$grpcode){
#   print(paste("Processing Insurer",g,"of",lob))
#   co_model=model_function(g)
#   stats=rbind(stats,data.frame(lob,co_model$SumStats))
#   write.csv(stats,file=outfile)
#   t1=Sys.time()
#   print(t1-t0)
# }
#
