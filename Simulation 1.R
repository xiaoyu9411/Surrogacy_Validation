#1 and 2 years

#generate data

library(MASS)
library(survival)
library(survRM2)
library(nlme)
library(lme4)
library(mvmeta)


R2t1<-0.90
R2t2<-0.90

vart1<-1
vart2<-1
covt1<-0.75
covt2<-0.75
ntrial<-30 #number of trials in total
tau<-c(1,2)
varbij<-0.2
varbeij<-0.3

vare<-0.1
trialtau2<-30 #number of trials avaliable at second time point


D1<-0.1*matrix(c(vart1,covt1,0,0,covt1,vart1,0,0,0,0,vart1,sqrt(R2t1),0,0,sqrt(R2t1),vart1),
               nrow=4, ncol=4)
D2<-0.4*matrix(c(vart2,covt2,0,0,covt2,vart2,0,0,0,0,vart2,sqrt(R2t2),0,0,sqrt(R2t2),vart2),
               nrow=4, ncol=4)

R<-matrix(c(0.05,0.04,0,0,0.04,0.05,0,0,0,0,0.07,0.06,0,0,0.06,0.07),nrow = 4,ncol=4)
meaneff<-c(10,8,2,3,20,18,3,5)  #12, 24 months

D<-rbind(cbind(D1,R),cbind(t(R),D2))
mu<-rep(0,8)
#is.positive.definite(D)



nsim<-5000
R2all<-matrix(nrow=nsim,ncol=length(tau)*3)
colnames(R2all)<-c("R21","R21U","R21L","R22","R22U","R22L")
numtau<-length(tau)
numend<-2


it<-0
for (l in 1:nsim){
  it<-it+1
  set.seed(l*11)
  
  betaall<-matrix(NA,ncol=numtau*numend*2,nrow=ntrial)
  covall<-list()
  for (i in 1:ntrial){
    covall[[i]]<-matrix(0,ncol=numtau*numend*2,nrow=numtau*numend*2)
  }
  
  for ( j in 1:ntrial){
    theta<-matrix(,nrow=1,ncol=6+length(tau))
    n <- 500
    randomeff<-mvrnorm(n=1,mu=mu, Sigma=D)#mt1,ms1,at1,bt1,mt2,ms2,at2,bt2
    
    beta<-meaneff+randomeff
    trt <- rbinom(n,1,0.5) # trt group
    timeT<-c()
    timeS<-c()
    
    glmm<-matrix(,,ncol=6+length(tau))
    colnames(glmm)<-c("theta","id","trt","It","Is","Itau1","Itau2","endpoint")
    for (i in 1:n){
      meanYTt1<-beta[1]+beta[3]*trt[i]
      meanYTt2<-beta[5]+beta[7]*trt[i]
      meanYSt1<-beta[2]+beta[4]*trt[i]
      meanYSt2<-beta[6]+beta[8]*trt[i]
      ds<-matrix(,nrow=4,ncol=6+length(tau))
      bij<-rnorm(1,mean = 0,sd=sqrt(varbij))
      bsij<-rnorm(1,mean = 0,sd=sqrt(varbeij))
      btij<-rnorm(1,mean = 0,sd=sqrt(varbeij))
      ds[1,1]<-meanYTt1+bij+btij+rnorm(1,mean = 0,sd=sqrt(vare))
      ds[2,1]<-meanYTt2+bij+btij+rnorm(1,mean = 0,sd=sqrt(vare))
      ds[3,1]<-meanYSt1+bij+bsij+rnorm(1,mean = 0,sd=sqrt(vare))
      ds[4,1]<-meanYSt2+bij+bsij+rnorm(1,mean = 0,sd=sqrt(vare))
      ds[,2]<-rep(i,nrow(ds))
      ds[,3]<-rep(trt[i],nrow(ds))
      ds[,4]<-c(1,1,0,0)
      ds[,5]<-c(0,0,1,1)
      ds[,6]<-c(1,0,1,0)
      ds[,7]<-c(0,1,0,1)
      ds[,8]<-c(1,1,2,2)
      glmm<-rbind(glmm,ds)
    }
    glmm<-glmm[-1,]
    glmm<-as.data.frame(glmm)
    numtau<-length(tau)
    numend<-2
    
    
    fit<-lmer(theta~Itau1*It+Itau1*Is+
                trt*Itau1*It+trt*Itau1*Is+
                Itau2*It+Itau2*Is+
                trt*Itau2*It+trt*Itau2*Is+
                -1-trt-Itau1-Itau2-It-Is
              -trt*Itau1-trt*Itau2-trt*It-trt*Is+
                (1|id)+(1|id:endpoint),
              data=glmm,REML=FALSE)
    
    s<-summary(fit)
    
    betaall[j,1:(2*2*length(tau))]<-s$coefficients[c(1,2,5,6,3,4,7,8),1]
    covall[[j]][1:nrow(s$vcov),1:ncol(s$vcov)]<-as.matrix(s$vcov)
  }
  
  #when less trials are avaliable the the second time point
  #betaall[(trialtau2+1):ntrial,5:8]<-NA
  #for (m in (trialtau2+1):ntrial){
   # covall[[m]][5:8,5:8]<-NA
  #}
  

  
  
  REmod <- mvmeta(betaall, covall, method="mm")
  cov<-REmod$Psi
  R2t<-c()
  R2tu<-c()
  
  R2tl<-c()
  for (i in 1:length(tau)){
    cov1<-cov[((i-1)*4+1):((i-1)*4+4),((i-1)*4+1):((i-1)*4+4)]
    R2<-t(c(cov1[2,3],cov1[4,3]))%*%solve(matrix(c(cov1[2,2],cov1[2,4],cov1[4,2],cov1[4,4]),nrow=2,ncol=2))%*% c(cov1[2,3],cov1[4,3])
    R2t[i]<-R2/cov1[3,3]
    N<-nrow(betaall)
    R2tu[i]<-R2t[i]+1.96*sqrt((4*R2t[i]*(1-R2t[i])^2)/(N-3))
    R2tl[i]<-R2t[i]-1.96*sqrt((4*R2t[i]*(1-R2t[i])^2)/(N-3))
  }
  
  R2all[it,c(1,4)]<-R2t
  R2all[it,c(2,5)]<-R2tu
  R2all[it,c(3,6)]<-R2tl
  write.table(R2all[it,],file="result1.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}





