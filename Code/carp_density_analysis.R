
#load data
den_dat<-read.csv("Data/density_data_2012_2022.csv")
#factor pool
den_dat$Pool<-factor(den_dat$Pool,levels = c("Alton","LaGrange","Peoria","Starved Rock","Marseilles","Dresden"))
#create variable for upper pool 
upper<-rep(0,dim(den_dat)[1])
upper[den_dat$Pool%in%levels(den_dat$Pool)[4:6]]<-1
#add upper to data.frame
df<-data.frame(den_dat,upper)
#reorder and 
df<-df[order(df$Pool,df$Year),c("Year","Pool","SVCP","BHCP","upper")]
min_date<-which(df$Year==2012)
#find max data for each pool
tapply(den_dat$Year,INDEX = den_dat$Pool,FUN = max,na.rm=T)
max_date<-which(df$Year==2022)

#create a data.frame with adjacent years in same row
dat_skeleton<-list()

for(i in levels(df$Pool)) {
  dat_skeleton[[i]] <- data.frame(
    pool = df[which(df$Pool == i& df$Year != "2022"), "Pool"],
    year = df[which(df$Pool == i& df$Year != "2022"), "Year"],
    sv_t0 = df[which(df$Pool ==i & df$Year != "2022"), "SVCP"],
    sv_t1 = df[which(df$Pool ==i & df$Year != "2012"), "SVCP"],
    bh_t0 = df[which(df$Pool ==i & df$Year != "2022"), "BHCP"],
    bh_t1 = df[which(df$Pool ==i & df$Year != "2012"), "BHCP"],
    upper = df[which(df$Pool ==i & df$Year != "2012"), "upper"]
  )
}
#combine
df_adj<-as.data.frame(do.call(rbind,dat_skeleton))
df_adj$upper<-factor(df_adj$upper)

#remove all NA
df_adj<-df_adj[-which(apply(df_adj,1,function(x)any(is.na(x)))),]
names(df_adj)
#set zeta as the constant for adding to log(Bighead Carp+zeta) 
zeta=0.01

####Silver Carp####
#fit a log-linear Ricker model to Silver Carp as directed by M_1 and M_3
summary(f_sv_m1_m3<-lm(log(sv_t1)~offset(log(sv_t0))+sv_t0,data = df_adj))
#fit a log-linear Ricker model to Silver Carp as directed by M_2 and M_4
summary(f_sv_m2_m4<-lm(log(sv_t1)~offset(log(sv_t0))+sv_t0+bh_t0,data = df_adj))

####Bighead Carp####      
#fit a log-linear ricker model to Bighead Carp as directed by M_1 and M_2
summary(f_bh_m1_m2<-lm(log(bh_t1+zeta)~offset(log(bh_t0+zeta))+bh_t0,data = df_adj))
#fit a log-linear ricker model to Bighead Carp as directed by M_3 and M_4
summary(f_bh_m3_m4<-lm(log(bh_t1+zeta)~offset(log(bh_t0+zeta))+bh_t0+sv_t0,data = df_adj))

####Model Selection####
#best of the 4 models can be tested using two tests
#(a) P(sv|M_2)>P(sv|M_1) AND (b) P(bh|M_3)>P(bh|M_1)
#if if(!a&!b)->M_1; if(a~!b)->M_2; if(~a&b)->M_3; if(a&b)->M_4 
#test of (a): this gives a p-value of the null hypothesis that the simpler model is no different
p_a<-1-pchisq(chi_a<--2*(logLik(f_sv_m1_m3)[1]-logLik(f_sv_m2_m4)[1]),1)
#test of (b): this gives a p-value of the null hypothesis that the simpler model is no different
p_b<-1-pchisq(chi_b<--2*(logLik(f_bh_m1_m2)[1]-logLik(f_bh_m3_m4)[1]),1)
#results give Wilkes Theorem and alpha=0.05
p_a<=0.05
p_b<=0.05
###(F,T) supports model M_3 ###
chi_a
chi_b
M_1_D<--2*(logLik(f_bh_m1_m2)[1]+logLik(f_sv_m1_m3)[1])
M_2_D<--2*(logLik(f_bh_m1_m2)[1]+logLik(f_sv_m2_m4)[1])
M_3_D<--2*(logLik(f_bh_m3_m4)[1]+logLik(f_sv_m1_m3)[1])
M_4_D<--2*(logLik(f_bh_m3_m4)[1]+logLik(f_sv_m2_m4)[1])
1-pchisq(M_1_D-M_2_D,1)
1-pchisq(M_1_D-M_3_D,1)
1-pchisq(M_1_D-M_4_D,2)
anova(f_sv_m1_m3,f_sv_m2_m4,test="Chisq")

#calculate rho (i.e., residual correlation)
cor.test(exp(resid(f_bh_m3_m4)),exp(resid(f_sv_m1_m3)))
cor.test((resid(f_bh_m3_m4)),(resid(f_sv_m1_m3)))
#monte carlo resampling estimate of uncertainty around rho
library(mvtnorm)
pred_sv<-cbind(rep(1,5000),rmvnorm(5000,coef(f_sv_m1_m3),vcov(f_sv_m1_m3)))%*%t(cbind(f_sv_m1_m3$model[,2],rep(1,54),f_sv_m1_m3$model[,3]))
pred_bh<-cbind(rep(1,5000),rmvnorm(5000,coef(f_bh_m3_m4),vcov(f_bh_m3_m4)))%*%t(cbind(f_bh_m3_m4$model[,2],rep(1,54),f_bh_m3_m4$model[,3:4]))
res_sv<-sweep(exp(t(pred_sv)),1,exp(f_sv_m1_m3$model[,1]),FUN = "-")
res_bh<-sweep(exp(t(pred_bh)),1,exp(f_bh_m3_m4$model[,1]),FUN = "-")
hist(diag(cor(res_sv,res_bh)))
quantile(diag(cor(res_sv,res_bh)),c(.025,.5,.975))

#same analysis with SEM
library(lavaan)
mod<-"
      lsv_t1~1*lsv_t0+sv_t0
      lbh_t1~1*lbh_t0+bh_t0
      lsv_t1~~lbh_t1
      sv_t0~~bh_t0
      "
dat<-data.frame(lsv_t1=log(df_adj$sv_t1),lsv_t0=log(df_adj$sv_t0),sv_t0=df_adj$sv_t0,sv_t1=df_adj$sv_t1,
                bh_t1=df_adj$bh_t1,
                lbh_t1=log(df_adj$bh_t1+zeta),lbh_t0=log(df_adj$bh_t0+zeta),bh_t0=df_adj$bh_t0)
fit<-sem(mod,data = dat,meanstructure = TRUE)
sem_coef<-parameterestimates(fit)
sem_coef[which(sem_coef$lhs=="lsv_t1"&sem_coef$rhs=="sv_t0"),c("ci.lower","est","ci.upper")]
mod1<-"
      lsv_t1~1*lsv_t0+sv_t0
      lbh_t1~sv_t0+1*lbh_t0+bh_t0
      lsv_t1~~lbh_t1
      sv_t0~~bh_t0
      "
1-.30/cor(df_adj$sv_t1,df_adj$bh_t1)

fit1<-sem(mod1,data = dat,meanstructure = TRUE)
summary(fit1)
sem_coef<-parameterestimates(fit1)
sqrt(diag(vcov(fit1)))
sem_coef[c(16,5,3,15,2),"est"]
vcov(fit1)[c(11,3,2,10,1),c(11,3,2,10,1)]

pred_sv_sem<-cbind(rep(1,5000),rmvnorm(5000,sem_coef[c(15,2),"est"],vcov(fit1)[c(10,1),c(10,1)]))%*%t(cbind(f_sv_m1_m3$model[,2],rep(1,54),f_sv_m1_m3$model[,3]))
pred_bh_sem<-cbind(rep(1,5000),rmvnorm(5000,sem_coef[c(16,5,3),"est"],vcov(fit1)[c(11,3,2),c(11,3,2)]))%*%t(cbind(f_bh_m3_m4$model[,2],rep(1,54),f_bh_m3_m4$model[,3:4]))
res_sv_sem<-sweep(exp(t(pred_sv_sem)),1,exp(f_sv_m1_m3$model[,1]),FUN = "-")
res_bh_sem<-sweep(exp(t(pred_bh_sem)),1,exp(f_bh_m3_m4$model[,1]),FUN = "-")
hist(diag(cor(res_sv_sem,res_bh_sem)))
quantile(diag(cor(res_sv_sem,res_bh_sem)),c(.025,.5,.975))

#corr
0.170/(sqrt(0.503)*sqrt(0.535))
0.007/(sqrt(0.054)*sqrt(0.004))

####Bayesian non-linear fit###
library(rjags)
mod<-'model{
for(i in 1:n){
y[i,]~dmnorm.vcov(mu[i,],T[,])
#sv_t1[i]~dnorm(mu_sv[i],tau_sv)
mu[i,1]<-b[1]*sv_t0[i]*exp(b[2]*sv_t0[i])
mu[i,2]<-a[1]*bh_t0[i]*exp(a[2]*bh_t0[i]+a[3]*sv_t0[i])
}

b[1] ~ dnorm(0,.001)
b[2] ~ dnorm(0,.001)
for(j in 1:3){
a[j] ~ dnorm(0,.001)
}
lambda[1]~dgamma(.001,.001)
lambda[2]~dgamma(.001,.001)
r_xy~dunif(-1,1)
sigma[1]<-1/sqrt(lambda[1])
sigma[2]<-1/sqrt(lambda[2])
T[1,1]<-1/lambda[1]
T[1,2]<-r_xy*sigma[1]*sigma[2]
T[2,1]<-r_xy*sigma[1]*sigma[2]
T[2,2]<-1/lambda[2]
}'
#jdat$bh_ob[which(jdat$bh_ob==0)]<-NA
jdat<-list('n'=dim(df_adj)[1],"y"=cbind(df_adj$sv_t1,df_adj$bh_t1),
           "bh_t0"=df_adj$bh_t0,"sv_t0"=df_adj$sv_t0)
tot.mod<-jags.model(textConnection(mod),data = jdat,n.chains = 5,n.adapt = 1000)
#burn-in
update(tot.mod,1000)
#get samples from it for what we are interested V, S and the obs errors
fit.mod<-coda.samples(tot.mod,c('b','a','r_xy'),n.iter = 5000)
#save these for future use
saveRDS(fit.mod,"Output/parameter_posterior_samples.RData")
#get predictions
fit.mu<-coda.samples(tot.mod,c('mu'),n.iter = 5000)
saveRDS(fit.mu,"Output/prediction_posterior_samples.RData")
mu_samp<-summary(fit.mu)

####this joins results from the different analyses####
#int sv
int_sv<-rbind(exp(sem_coef[which(sem_coef$lhs=="lsv_t1"&sem_coef$op=="~1"),c("ci.lower","est","ci.upper")]), #sem
c(exp(coef(f_sv_m1_m3)[1]-1.96*sqrt(vcov(f_sv_m1_m3)[1,1])),exp(coef(f_sv_m1_m3)[1]),exp(coef(f_sv_m1_m3)[1]+1.96*sqrt(vcov(f_sv_m1_m3)[1,1]))),
summary(fit.mod)$quantiles["b[1]",c("2.5%","50%","97.5%")])

#int bh
int_bh<-rbind(exp(sem_coef[which(sem_coef$lhs=="lbh_t1"&sem_coef$op=="~1"),c("ci.lower","est","ci.upper")]), #sem
c(exp(coef(f_bh_m3_m4)[1]-1.96*sqrt(vcov(f_bh_m3_m4)[1,1])),exp(coef(f_bh_m3_m4)[1]),exp(coef(f_bh_m3_m4)[1]+1.96*sqrt(vcov(f_bh_m3_m4)[1,1]))),
summary(fit.mod)$quantiles["a[1]",c("2.5%","50%","97.5%")])

#intra sv
intra_sv<-rbind((sem_coef[which(sem_coef$lhs=="lsv_t1"&sem_coef$rhs=="sv_t0"),c("ci.lower","est","ci.upper")]),
c(coef(f_sv_m1_m3)[2]-1.96*sqrt(vcov(f_sv_m1_m3)[2,2]),coef(f_sv_m1_m3)[2],coef(f_sv_m1_m3)[2]+1.96*sqrt(vcov(f_sv_m1_m3)[2,2])),
summary(fit.mod)$quantiles["b[2]",c("2.5%","50%","97.5%")])

#intra bh
intra_bh<-rbind((sem_coef[which(sem_coef$lhs=="lbh_t1"&sem_coef$rhs=="bh_t0"),c("ci.lower","est","ci.upper")]),
c((coef(f_bh_m3_m4)[2]-1.96*sqrt(vcov(f_bh_m3_m4)[2,2])),(coef(f_bh_m3_m4)[2]),(coef(f_bh_m3_m4)[2]+1.96*sqrt(vcov(f_bh_m3_m4)[2,2]))),
summary(fit.mod)$quantiles["a[2]",c("2.5%","50%","97.5%")])

#inter bh<-sv
inter_bh<-rbind((sem_coef[which(sem_coef$lhs=="lbh_t1"&sem_coef$rhs=="sv_t0"),c("ci.lower","est","ci.upper")]),
c((coef(f_bh_m3_m4)[3]-1.96*sqrt(vcov(f_bh_m3_m4)[3,3])),(coef(f_bh_m3_m4)[3]),(coef(f_bh_m3_m4)[3]+1.96*sqrt(vcov(f_bh_m3_m4)[3,3]))),
summary(fit.mod)$quantiles["a[3]",c("2.5%","50%","97.5%")])

#corr
corr_et1<-rbind(quantile(diag(cor(res_sv_sem,res_bh_sem)),c(.025,.5,.975)),quantile(diag(cor(res_sv,res_bh)),c(.025,.5,.975)),summary(fit.mod)$quantiles["r_xy",c("2.5%","50%","97.5%")])

1-pchisq(-2*(logLik(fit)[1]-logLik(fit1)[1]),1)

####look at goodness of fit####
cor(df_adj$sv_t1,mu_samp$quantiles[55:108,"50%"])
cor(df_adj$sv_t1,mu_samp$quantiles[1:54,"50%"])
plot(df_adj$sv_t1,mu_samp$quantiles[1:54,"50%"])
abline(0,1)
points(df_adj$bh_t1,mu_samp$quantiles[55:108,"50%"],pch=21,bg="lightblue")
abline(0,1)

#look at correlation among estimates
round(cov(do.call(rbind,fit.mod)[,1:5]),9)
###create Figure 1: Observed dynamics###

pool<-unique(df_adj$pool)
pretty_pool<-c("Alton","La Grange", "Peoria", "Starved Rock","Marseilles" , "Dresden Island")
names(pretty_pool)<-pool
cols<-c("#41afaa",
        "#466eb4",
        "#00a0e1",
        "#e6a532",
        "#d7642c",
        "#af4b91")
names(cols)<-pool
ylims<-c(4,4,4,1,1,1)
names(ylims)<-pool

pdf("Output/Fig1.pdf")
  par(mfrow=c(2,1))
  #plot (a) Silver carp. initiate plot with first pool 
  plot(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$SVCP[den_dat$Pool%in%pool[1]],type='b', col=cols[1],ylim=c(0,4),
     ylab="",xlab="Year",bty="l",main="Silver Carp")
    arrows(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$SVCP[den_dat$Pool%in%pool[1]],den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$SVCP_Low_95CI[den_dat$Pool%in%pool[1]],
       length = 0)
    arrows(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$SVCP[den_dat$Pool%in%pool[1]],den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$SVCP_Up_95CI[den_dat$Pool%in%pool[1]],
       length = 0)
    mtext(side = 2,text = bquote("Density (No./1000"~m^3~")"),padj = -2.5)
    #add the other pools
    for(i in pool[-1]){
      lines(den_dat$Year[den_dat$Pool%in%i],den_dat$SVCP[den_dat$Pool%in%i],type='b', col=cols[i])
      arrows(den_dat$Year[den_dat$Pool%in%i],den_dat$SVCP[den_dat$Pool%in%i],den_dat$Year[den_dat$Pool%in%i],den_dat$SVCP_Low_95CI[den_dat$Pool%in%i],
         length = 0)
      arrows(den_dat$Year[den_dat$Pool%in%i],den_dat$SVCP[den_dat$Pool%in%i],den_dat$Year[den_dat$Pool%in%i],den_dat$SVCP_Up_95CI[den_dat$Pool%in%i],
         length = 0)
    }
    legend("topright",legend = pretty_pool[pool],col = cols,bty='n',lty=1)
    mtext("(a)",side=3,adj=-.1)
    #plot (b) Bighead carp. initiate plot with first pool 
    plot(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$BHCP[den_dat$Pool%in%pool[1]],type='b', col=cols[1],ylim=c(0,1),
     ylab="",xlab="Year",bty="l",main="Bighead Carp")
      arrows(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$BHCP[den_dat$Pool%in%pool[1]],den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$BHCP_Low_95CI[den_dat$Pool%in%pool[1]],
       length = 0)
      arrows(den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$BHCP[den_dat$Pool%in%pool[1]],den_dat$Year[den_dat$Pool%in%pool[1]],den_dat$BHCP_Up_95CI[den_dat$Pool%in%pool[1]],
       length = 0)
      mtext(side = 2,text = bquote("Density (No./1000"~m^3~")"),padj = -2.5)
      #add the other pools
      for(i in pool[-1]){
        lines(den_dat$Year[den_dat$Pool%in%i],den_dat$BHCP[den_dat$Pool%in%i],type='b', col=cols[i])
        arrows(den_dat$Year[den_dat$Pool%in%i],den_dat$BHCP[den_dat$Pool%in%i],den_dat$Year[den_dat$Pool%in%i],den_dat$BHCP_Low_95CI[den_dat$Pool%in%i],
         length = 0)
        arrows(den_dat$Year[den_dat$Pool%in%i],den_dat$BHCP[den_dat$Pool%in%i],den_dat$Year[den_dat$Pool%in%i],den_dat$BHCP_Up_95CI[den_dat$Pool%in%i],
         length = 0)
      }
      mtext("(b)",side=3,adj=-.1)
dev.off()

## create Figure 4: predicted vs observed plots
#get correlations
sv_cor<-round(cor(df_adj$sv_t1,exp(predict(f_sv_m1_m3))),3)
bh_cor<-round(cor(df_adj$bh_t1,exp(predict(f_bh_m3_m4))),3)

pdf("Output/Fig4.pdf",height = 5,width = 8)
  par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(5,4,1,1))
#plot (a) Silver carp
    plot(df_adj$sv_t1~exp(predict(f_sv_m1_m3)-zeta),pch=21,bg="grey",bty='l',
     main="Silver Carp",ylab="",xlab="Predicted")
    mtext(side = 2,text = bquote("Density (No./1000"~m^3~")"),padj = -2.5)
    legend(x = 0,y=1.5,legend = bquote(rho == .(sv_cor)),bty='n')
    abline(0,1)
    mtext("(a)",side=3,adj=-0.2)
#plot (b) Bighead carp   
  plot(df_adj$bh_t1~exp(predict(f_bh_m3_m4)-zeta),pch=21,bg="grey",bty='l',
     main="Bighead Carp",xlab="Predicted",ylab="")
    abline(0,1)
    legend(x=0,y=0.45,legend = bquote(rho == .(bh_cor)),bty='n')
    mtext("(b)",side=3,adj=-0.2)
dev.off()

####create Figure 5: parameter values plot####
labs=c("A",bquote(alpha[0]),bquote(alpha[1]),"B",bquote(beta[0]),bquote(rho[E[t1]]))

pdf("Output/Fig5.pdf")
#regression
  plot(c(int_bh[2,2],intra_bh[2,2],inter_bh[2,2],int_sv[2,2],intra_sv[2,2],corr_et1[2,2]),seq(.8,5.8,by=1),ylab="Parameter",
     main="",ylim=c(0,6.5),xlim=c(-7,2),pch=21,bty="l",xaxt='n',yaxt='n',xlab="Value",cex.lab=1.25)
  arrows(c(int_bh[2,2],intra_bh[2,2],inter_bh[2,2],int_sv[2,2],intra_sv[2,2],corr_et1[2,2]),seq(.8,5.8,by=1),
       c(int_bh[2,1],intra_bh[2,1],inter_bh[2,1],int_sv[2,1],intra_sv[2,1],corr_et1[2,1]),seq(.8,5.8,by=1),length=0,lwd=2)
  arrows(c(int_bh[2,2],intra_bh[2,2],inter_bh[2,2],int_sv[2,2],intra_sv[2,2],corr_et1[2,2]),seq(.8,5.8,by=1),
       c(int_bh[2,3],intra_bh[2,3],inter_bh[2,3],int_sv[2,3],intra_sv[2,3],corr_et1[2,3]),seq(.8,5.8,by=1),length=0,lwd=2)
  abline(v=0,lty=2)
  points(c(int_bh[2,2],intra_bh[2,2],inter_bh[2,2],int_sv[2,2],intra_sv[2,2],corr_et1[2,2]),
       seq(.8,5.8,by=1),pch=21,bg=c(rep("darkgrey",3),rep("darkgrey",2)))
#add sem
  points(c(int_bh[1,2],intra_bh[1,2],inter_bh[1,2],int_sv[1,2],intra_sv[1,2],corr_et1[1,2]),1:6,pch=22)
  arrows(c(int_bh[1,2],intra_bh[1,2],inter_bh[1,2],int_sv[1,2],intra_sv[1,2],corr_et1[1,2]),1:6,
       c(int_bh[1,1],intra_bh[1,1],inter_bh[1,1],int_sv[1,1],intra_sv[1,1],corr_et1[1,1]),1:6,length=0,lwd=2)
  arrows(c(int_bh[1,2],intra_bh[1,2],inter_bh[1,2],int_sv[1,2],intra_sv[1,2],corr_et1[1,2]),1:6,
       c(int_bh[1,3],intra_bh[1,3],inter_bh[1,3],int_sv[1,3],intra_sv[1,3],corr_et1[1,3]),1:6,length=0,lwd=2)
  points(c(int_bh[1,2],intra_bh[1,2],inter_bh[1,2],int_sv[1,2],intra_sv[1,2],corr_et1[1,2]),
       1:6,pch=22,bg=c(rep("darkgrey",3),rep("darkgrey",2)))
#add nl-bayes
  points(c(int_bh[3,2],intra_bh[3,2],inter_bh[3,2],int_sv[3,2],intra_sv[3,2],corr_et1[3,2]),seq(1.2,6.2,by=1),pch=23)
  arrows(c(int_bh[3,2],intra_bh[3,2],inter_bh[3,2],int_sv[3,2],intra_sv[3,2],corr_et1[3,2]),seq(1.2,6.2,by=1),
       c(int_bh[3,1],intra_bh[3,1],inter_bh[3,1],int_sv[3,1],intra_sv[3,1],corr_et1[3,1]),seq(1.2,6.2,by=1),length=0,lwd=2)
  arrows(c(int_bh[3,2],intra_bh[3,2],inter_bh[3,2],int_sv[3,2],intra_sv[3,2],corr_et1[3,2]),seq(1.2,6.2,by=1),
       c(int_bh[3,3],intra_bh[3,3],inter_bh[3,3],int_sv[3,3],intra_sv[3,3],corr_et1[3,3]),seq(1.2,6.2,by=1),length=0,lwd=2)
  points(c(int_bh[3,2],intra_bh[3,2],inter_bh[3,2],int_sv[3,2],intra_sv[3,2],corr_et1[3,2]),
       seq(1.2,6.2,by=1),pch=23,bg=c(rep("darkgrey",3),rep("darkgrey",2)))
  axis(side = 2,at = 1:6,labels =do.call(expression,labs),las=1,cex.axis=1.25)
  axis(side = 1, at=seq(-6,2,1),seq(-6,2,1),cex.axis=1.25)
  legend("topleft",legend = c("Log-Linear Regression","SEM","Non-Linear Bayes"),pch=c(21,22,23),bty="n",pt.bg = "darkgrey")
dev.off()


####create figure 6: dynamics plot####
#load simulation results
sv_dyn<-read.csv("Data/sv_dyn.csv",header = T)
bh_dyn<-read.csv("Data/bh_dyn.csv",header = T)

svlw_dyn<-read.csv("Data/svlw_dyn.csv",header = T)
bhlw_dyn<-read.csv("Data/bhlw_dyn.csv",header = T)

svhg_dyn<-read.csv("Data/svhg_dyn.csv",header = T)
bhhg_dyn<-read.csv("Data/bhhg_dyn.csv",header = T)

pdf("Output/Fig6.pdf",height = 5,width = 10)
  par(mfrow=c(1,3),oma=c(3,4,0,0))
  #add plot (a) observed Silver carp growth rate
  plot(1:101,sv_dyn[,2],lwd=2,ylim=c(0,.25),type="l",ylab="",xlab='Time',bty="l",main="Observed SC Growth Rate (B)",cex.lab=1.25,cex.axis=1.25)
    polygon(c(1:101,rev(1:101)),c(sv_dyn[,1],rev(sv_dyn[,3])),col=adjustcolor("grey",alpha.f=0.5),lwd=.5,lty=1 )
    polygon(c(1:101,rev(1:101)),c(bh_dyn[,1],rev(bh_dyn[,3])),col=adjustcolor("lightblue",alpha.f=0.5),lwd=.5,lty=2)
    lines(1:101,bh_dyn[,2],lwd=2,lty=2,ylim=c(0,.5),type="l")
    lines(1:101,sv_dyn[,2],lwd=2,lty=1,ylim=c(0,.5),type="l")
    legend("topleft",legend = c("Silver Carp","Bighead Carp"),lwd=2,lty=c(1,2),bty='n',cex = 1.25)
    mtext(side = 2,text = bquote("Density (No./1000"~m^3~")"),padj = -2.5)
    mtext("(a)",side = 3,adj=-0.3)
  #add plot (b) reduced Silver carp growth rate
  plot(1:101,svlw_dyn[,2],lwd=2,ylim=c(0,.25),type="l",ylab="",xlab='Time',bty="l",main="Reduced SC Growth Rate (B)",cex.lab=1.25,cex.axis=1.25)
    polygon(c(1:101,rev(1:101)),c(svlw_dyn[,1],rev(svlw_dyn[,3])),col=adjustcolor("grey",alpha.f=0.5),lwd=.5,lty=1 )
    polygon(c(1:101,rev(1:101)),c(bhlw_dyn[,1],rev(bhlw_dyn[,3])),col=adjustcolor("lightblue",alpha.f=0.5),lwd=.5,lty=2)
    lines(1:101,bhlw_dyn[,2],lwd=2,lty=2,ylim=c(0,.5),type="l")
    lines(1:101,svlw_dyn[,2],lwd=2,lty=1,ylim=c(0,.5),type="l")
    mtext("(b)",side = 3,adj=-0.3)
    
  #add plot (c) increased Silver carp growth rate
  plot(1:101,svhg_dyn[,2],lwd=2,ylim=c(0,.25),type="l",ylab="",xlab='Time',bty="l",main="Increased SC Growth Rate (B)",cex.lab=1.25,cex.axis=1.25)
    polygon(c(1:101,rev(1:101)),c(svhg_dyn[,1],rev(svhg_dyn[,3])),col=adjustcolor("grey",alpha.f=0.5),lwd=.5,lty=1 )
    polygon(c(1:101,rev(1:101)),c(bhhg_dyn[,1],rev(bhhg_dyn[,3])),col=adjustcolor("lightblue",alpha.f=0.5),lwd=.5,lty=2)
    lines(1:101,bhhg_dyn[,2],lwd=2,lty=2,ylim=c(0,.5),type="l")
    lines(1:101,svhg_dyn[,2],lwd=2,lty=1,ylim=c(0,.5),type="l")
    mtext("(c)",side = 3,adj=-0.3)
dev.off()

