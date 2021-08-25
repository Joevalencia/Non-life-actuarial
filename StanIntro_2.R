#
# Gamma example with covariates for ASTIN 2017 Workshop- Glenn Meyers
#
# preliminary housekeeping
#
rm(list = ls())      # clear workspace
t0=Sys.time()
setwd("~/Dropbox/ASTIN 2017")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
#
# inputs
#
stan_thin=1
set.seed(123456)
#
# get simulated data
#
num.obs=500
x=runif(num.obs)
beta=1+x
y=rgamma(num.obs,1,beta)
data.model1=list(num_obs  = num.obs,
                 y        = y)
#
# Stan script
#
scode = "
data{
int<lower=1> num_obs;
real<lower=0> y[num_obs];
}
parameters{
real<lower=0> beta;
real<lower=0> alpha;
}
model {
alpha ~ gamma(1,1);
beta~ gamma(1,1);
y ~ gamma(alpha,beta);
}
generated quantities{
vector[num_obs] log_lik;
for (i in 1:num_obs) log_lik[i] = gamma_lpdf(y[i]|alpha,beta);
}
"
#
# initialization function for scode
#
init.parms=function(chain_id){
  set.seed(123+chain_id)
  list(alpha =rgamma(1,shape=1,rate=1),
       beta=rgamma(1,shape=1,rate=1))
}
#
# list of parameters in output
#
pars.list=c("alpha","beta","log_lik")
fit.model1 = stan(model_code=scode,data=data.model1,
                  seed = 123,iter=5000*stan_thin,thin=stan_thin,
                  #              control=list(adapt_delta=.999,max_treedepth=50),
                  chains = 4,pars=pars.list,verbose=FALSE)
fit.model1_summary=as.matrix(summary(fit.model1)$summary)[1:2,c(1,3,10)]
max_Rhat=round(max(fit.model1_summary[,3]),4)
print(fit.model1_summary)
print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
#
# extract posterior distribution
#
b=extract(fit.model1,pars=pars.list)
alpha=b$alpha
beta=b$beta
pairs(data.frame(alpha,beta))
#
# goodness of fit statistics for comparing models
#
loglik1=extract_log_lik(fit.model1)
loo1 <- loo(loglik1)
print(loo1)
#
t1=Sys.time()
print(t1-t0)
loom=matrix(0,3,2)
loom[1,1]=loo1$elpd_loo
loom[1,2]=loo1$se_elpd_loo
loom[2,1]=loo1$p_loo
loom[2,2]=loo1$se_p_loo
loom[3,1]=loo1$looic
loom[3,2]=loo1$se_looic