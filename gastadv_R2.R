library(surrosurv)
library(survival)
library(survRM2)
library(mvmeta)
data("gastadv") #load the data
gastadv$trt<-ifelse(gastadv$trt==0.5,1,0)
gastadv$timeS<-gastadv$timeS
gastadv$timeT<-gastadv$timeT


tau<-c(6*30,12*30,18*30, 24*30) #time points of interest
trial<-unique(gastadv$trialref)


numtau<-length(tau) #number of time points
numend<-2 #two endpoints

beta<-matrix(NA,ncol=numtau*numend*2,nrow=length(trial))#store the trial-specific coefficients
covall<-list() #store the within-study covariance matrix

for (i in 1:length(trial)){
  covall[[i]]<-matrix(0,ncol=numtau*numend*2,nrow=numtau*numend*2)
}

#first stage
for ( i in 1:length(trial)){
  #overall RMST:
  dat<-gastadv[gastadv$trialref==trial[i],]
  fut<-min(max(dat[which(dat$trt==0),]$timeT), max(dat[which(dat$trt==1),]$timeT))
  fus<-min(max(dat[which(dat$trt==0),]$timeS), max(dat[which(dat$trt==1),]$timeS))
  fu<-min(fut,fus)
  tautrial<-tau[tau<=fu]
  #pseudo observations:
  pseudot<-rep(NA,length(dat$id)*length(tautrial))
  
  objt<-survfit(Surv(timeT, statusT)~1, data=dat)
  overallt<-rep(NA,length(tautrial))
  for (k in 1:length(tautrial)){
    overallt[k]<-summary(objt, rmean=tautrial[k], print.rmean=T)$table[5]
  }
  for (j in 1:length(dat$id)){
    index<-dat$id[j]
    sampt<-dat[dat$id!=index,]
    objpartt<-survfit(Surv(timeT, statusT)~1, data=sampt)
    for (k in 1:length(tautrial)){
      partt<-summary(objpartt, rmean=tautrial[k], print.rmean=T)$table[5]
      pseudot[(j-1)*length(tautrial)+k]<-length(dat$id)*overallt[k]-(length(dat$id)-1)*partt
    }
  }
  
  
  pseudos<-rep(NA,length(dat$id)*length(tautrial))
  objs<-survfit(Surv(timeS, statusS)~1, data=dat)
  overalls<-rep(NA,length(tautrial))
  for (k in 1:length(tautrial)){
    overalls[k]<-summary(objs, rmean=tautrial[k], print.rmean=T)$table[5]
  }
  for (j in 1:length(dat$id)){
    index<-dat$id[j]
    samps<-dat[dat$id!=index,]
    objparts<-survfit(Surv(timeS, statusS)~1, data=samps)
    for (k in 1:length(tautrial)){
      parts<-summary(objparts, rmean=tautrial[k], print.rmean=T)$table[5]
      pseudos[(j-1)*length(tautrial)+k]<-length(dat$id)*overalls[k]-(length(dat$id)-1)*parts
    }
  }
  
  glmmdat<-as.data.frame(matrix(,ncol=9,nrow=2*nrow(dat)*length(tautrial)))
  colnames(glmmdat)<-c("theta","id","trt","It","Is","Itau1","Itau2","Itau3","Itau4")

  
  for ( k in 1:nrow(dat)){
    for (m in 1:length(tautrial)){
      glmmdat$theta[length(tautrial)*2*(k-1)+m]<-pseudot[(k-1)*length(tautrial)+m]
      glmmdat$id[length(tautrial)*2*(k-1)+m]<-k
      glmmdat$trt[length(tautrial)*2*(k-1)+m]<-dat$trt[k]
      glmmdat$It[length(tautrial)*2*(k-1)+m]<-1
      glmmdat$Is[length(tautrial)*2*(k-1)+m]<-0
      glmmdat[length(tautrial)*2*(k-1)+m,6]<-ifelse(m==1,1,0)
      glmmdat[length(tautrial)*2*(k-1)+m,6]<-ifelse(length(tautrial)>=1,glmmdat[length(tautrial)*2*(k-1)+m,6],NA)
      glmmdat[length(tautrial)*2*(k-1)+m,7]<-ifelse(m==2,1,0)
      glmmdat[length(tautrial)*2*(k-1)+m,7]<-ifelse(length(tautrial)>=2,glmmdat[length(tautrial)*2*(k-1)+m,7],NA)
      glmmdat[length(tautrial)*2*(k-1)+m,8]<-ifelse(m==3,1,0)
      glmmdat[length(tautrial)*2*(k-1)+m,8]<-ifelse(length(tautrial)>=3,glmmdat[length(tautrial)*2*(k-1)+m,8],NA)
      glmmdat[length(tautrial)*2*(k-1)+m,9]<-ifelse(m==4,1,0)
      glmmdat[length(tautrial)*2*(k-1)+m,9]<-ifelse(length(tautrial)>=4,glmmdat[length(tautrial)*2*(k-1)+m,9],NA)
      
      
      glmmdat$theta[length(tautrial)*2*(k-1)+length(tautrial)+m]<-pseudos[(k-1)*length(tautrial)+m]
      glmmdat$id[length(tautrial)*2*(k-1)+length(tautrial)+m]<-k
      glmmdat$trt[length(tautrial)*2*(k-1)+length(tautrial)+m]<-dat$trt[k]
      glmmdat$It[length(tautrial)*2*(k-1)+length(tautrial)+m]<-0
      glmmdat$Is[length(tautrial)*2*(k-1)+length(tautrial)+m]<-1
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,6]<-ifelse(m==1,1,0)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,6]<-ifelse(length(tautrial)>=1,glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,6],NA)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,7]<-ifelse(m==2,1,0)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,7]<-ifelse(length(tautrial)>=2,glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,7],NA)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,8]<-ifelse(m==3,1,0)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,8]<-ifelse(length(tautrial)>=3,glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,8],NA)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,9]<-ifelse(m==4,1,0)
      glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,9]<-ifelse(length(tautrial)>=4,glmmdat[length(tautrial)*2*(k-1)+length(tautrial)+m,9],NA)
      
      
      
      
       }

    
  }
  

  
  glmmdat$endpoint<-ifelse(glmmdat$It==1,1,2)
  for (l in 1:nrow(glmmdat)){
    if(!is.na(glmmdat$Itau1[l])){if(glmmdat$Itau1[l]==1){glmmdat$time[l]<-1}}
    if(!is.na(glmmdat$Itau2[l])){if(glmmdat$Itau2[l]==1){glmmdat$time[l]<-2}}
    if(!is.na(glmmdat$Itau3[l])){if(glmmdat$Itau3[l]==1){glmmdat$time[l]<-3}}
    if(!is.na(glmmdat$Itau4[l])){if(glmmdat$Itau4[l]==1){glmmdat$time[l]<-4}}
  }
  
  
  
 
library(nlme)
library(lme4)

  
  if(length(tautrial)==1){ 
    fit<-lmer(theta~Itau1*It+
               Itau1*Is+
               trt*Itau1*It+
               trt*Itau1*Is+
               -1-trt-Itau1-It-Is
             -trt*Itau1-trt*It-trt*Is+
               (1|id),
             data=glmmdat,REML = FALSE)   
    s<-summary(fit)
    
    beta[i,1:(2*2*length(tautrial))]<-s$coefficients[,1] 
    covall[[i]][1:nrow(s$vcov),1:ncol(s$vcov)]<-as.matrix(s$vcov)
  }  
  
  
  

  

  
  if(length(tautrial)==2){ 
    fit<-lmer(theta~Itau1*It+Itau1*Is+
                trt*Itau1*It+trt*Itau1*Is+
                Itau2*It+Itau2*Is+
               trt*Itau2*It+trt*Itau2*Is+
               -1-trt-Itau1-Itau2-It-Is
             -trt*Itau1-trt*Itau2-trt*It-trt*Is+
             (1|id)+(1|id:endpoint),
             data=glmmdat,REML=FALSE
             )   
    s<-summary(fit)
    
    beta[i,1:(2*2*length(tautrial))]<-s$coefficients[c(1,2,5,6,3,4,7,8),1]
    covall[[i]][1:nrow(s$vcov),1:ncol(s$vcov)]<-as.matrix(s$vcov)
  }
  
  
  
if(length(tautrial)==3){ 
    fit<-lmer(theta~Itau1*It+Itau1*Is+
                trt*Itau1*It+trt*Itau1*Is+
                Itau2*It+Itau2*Is+
                trt*Itau2*It+ trt*Itau2*Is+
               Itau3*It+Itau3*Is+
              trt*Itau3*It+trt*Itau3*Is+
             -1-trt-Itau1-Itau2-Itau3-It-Is
             -trt*Itau1-trt*Itau2-trt*Itau3-trt*It-trt*Is+
               (1|id)+(1|id:endpoint),
             data=glmmdat,REML=FALSE)   
    s<-summary(fit)
    
    beta[i,1:(2*2*length(tautrial))]<-s$coefficients[c(1,2,7,8,3,4,9,10,5,6,11,12),1]  
    covall[[i]][1:nrow(s$vcov),1:ncol(s$vcov)]<-as.matrix(s$vcov)
    
  }   
  
  
  
if(length(tautrial)==4){ 
fit<-lmer(theta~Itau1*It+Itau1*Is+
            trt*Itau1*It+trt*Itau1*Is+
            Itau2*It+Itau2*Is+
            trt*Itau2*It+trt*Itau2*Is+
           Itau3*It+Itau3*Is+
            trt*Itau3*It+trt*Itau3*Is+
            Itau4*It+Itau4*Is+
           trt*Itau4*It+trt*Itau4*Is
           -1-trt-Itau1-Itau2-Itau3-Itau4-It-Is
           -trt*Itau1-trt*Itau2-trt*Itau3-trt*Itau4-trt*It-trt*Is+
           (1|id)+(1|id:endpoint),,data=glmmdat,REML=FALSE)   
s<-summary(fit)

beta[i,1:(2*2*length(tautrial))]<-s$coefficients[c(1,2,9,10,3,4,11,12,5,6,13,14,7,8,15,16),1] 
covall[[i]][1:nrow(s$vcov),1:ncol(s$vcov)]<-as.matrix(s$vcov)

}  
  


  

}


#second stage
REmod <- mvmeta(beta, covall, method="mm",control = list(maxiter=1000))

cov<-REmod$Psi
R2t<-c()
R2tu<-c()


R2tl<-c()
for (i in 1:length(tau)){
cov1<-cov[((i-1)*4+1):((i-1)*4+4),((i-1)*4+1):((i-1)*4+4)]
R2<-t(c(cov1[2,3],cov1[4,3]))%*%solve(matrix(c(cov1[2,2],cov1[2,4],cov1[4,2],cov1[4,4]),nrow=2,ncol=2))%*% c(cov1[2,3],cov1[4,3])
R2t[i]<-R2/cov1[3,3]
R2tu[i]<-R2t[i]+1.96*sqrt((4*R2t[i]*(1-R2t[i])^2)/(N[i]-3))

R2tl[i]<-R2t[i]-1.96*sqrt((4*R2t[i]*(1-R2t[i])^2)/(N[i]-3))
}

R2t  #estimated R^2 at each time point
R2tu  #upper bound of analytical 95% CI
R2tl  #lower bound of analytical 95% CI

#Using Clayton copula model:

data("gastadv")
allSurroRes <- surrosurv(gastadv, c('clayton'), verbose = TRUE)





