###Simulate dynamics for Fig. 6.###

#load some libraries
library(rjags)
library(mvtnorm)

#load the posterior distributions from the Bayesian fie
post_param<-readRDS("Output/parameter_posterior_samples.RData")
#take a peek
summary(post_param)

#combine chains and get the mean and variance-covariance matrix of parameters
post_param_all<-do.call("rbind",post_param)
#set parameters of multivariate normal distribution as first five columns of parameters
MU<-apply(post_param_all[,-6],2,mean)
SIGMA<-cov(post_param_all[,-6])

#create 1000 sets of parameter drawing from these distributions
rand_param<-rmvnorm(1000,MU,SIGMA)

#create parameter set with reduced silver carp growth rate
rand_param_low<-rand_param
rand_param_low[,"b[1]"]<-rand_param_low[,"b[1]"]-0.1*rand_param_low[,"b[1]"]

#create parameter set with increased silver carp growth rate
rand_param_high<-rand_param
rand_param_high[,"b[1]"]<-rand_param_high[,"b[1]"]+0.1*rand_param_high[,"b[1]"]

#create a closure function to build a model that can be iterated given parameter values
f<-function(par)function(x)c(b=par[1]*x[1]*exp(par[2]*x[1]+par[3]*x[2]),s=par[4]*x[2]*exp(par[5]*x[2]))
#as a test, create a function d() from closure f() with the first row of simulated parameters
d<-f(rand_param[1,])
#give d() some initial conditions and get the output 
d(c(.1,.1))

#create a function to take a model, initial conditions and number of times to be iterated
#and returns all the steps
NestList <- function(f,x,n) {
  stopifnot(n>0)
  res <- matrix(x, nrow=n+1,ncol = 2)
  if (n == 1L) return(res)
  for (i in seq_len(n)) res[i+1,] <- f(res[i,])
  res
}


#####begin simulations#####

#create a list to hold answers and simulate 1000 time step for each parameter combination
# at estimated silver carp growth rate
ans<-list()
#run a 100 step simulation of joint dynamics for each parameter combination
for(i in 1:dim(rand_param)[1]){
  d<-f(rand_param[i,])
  ans[[i]]<-NestList(d,c(.1,.1),100)
}

#collect all the bighead carp sims
bh_ans<-sapply(ans,function(x)x[,1])
bh_dyn<-t(apply(bh_ans,1,quantile,c(.025,.5,.975)))
#collect all the silver carp sims
sv_ans<-sapply(ans,function(x)x[,2])
sv_dyn<-t(apply(sv_ans,1,quantile,c(.025,.5,.975)))


#do the same for each parameter combination at reduced silver carp growth rate
ans_low<-list()
for(i in 1:dim(rand_param_low)[1]){
  d<-f(rand_param_low[i,])
  ans_low[[i]]<-NestList(d,c(.1,.1),100)
}

#collect all the bighead carp sims
bhlw_ans<-sapply(ans_low,function(x)x[,1])
bhlw_dyn<-t(apply(bhlw_ans,1,quantile,c(.025,.5,.975)))
#collect all the silver carp sims
svlw_ans<-sapply(ans_low,function(x)x[,2])
svlw_dyn<-t(apply(svlw_ans,1,quantile,c(.025,.5,.975))) 

#finally, do the same for each parameter combination at increased silver carp growth rate
ans_high<-list()
for(i in 1:dim(rand_param_high)[1]){
  d<-f(rand_param_high[i,])
  ans_high[[i]]<-NestList(d,c(.1,.1),100)
}

#collect all the bighead carp sims
bhhg_ans<-sapply(ans_high,function(x)x[,1])
bhhg_dyn<-t(apply(bhhg_ans,1,quantile,c(.025,.5,.975)))

#collect all the silver carp sims
svhg_ans<-sapply(ans_high,function(x)x[,2])
svhg_dyn<-t(apply(svhg_ans,1,quantile,c(.025,.5,.975))) 

#save the output for later plotting
write.csv(sv_dyn,"Data/sv_dyn.csv",row.names = F)
write.csv(bh_dyn,"Data/bh_dyn.csv",row.names = F)
write.csv(svlw_dyn,"Data/svlw_dyn.csv",row.names = F)
write.csv(bhlw_dyn,"Data/bhlw_dyn.csv",row.names = F)
write.csv(svhg_dyn,"Data/svhg_dyn.csv",row.names = F)
write.csv(bhhg_dyn,"Data/bhhg_dyn.csv",row.names = F)

