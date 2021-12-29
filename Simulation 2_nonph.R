library(MASS)
library(survival)
library(survRM2)
library(nlme)
library(lme4)
library(mvmeta)
library(surrosurv)


R2 = 0.9
N1=10
N2=5
ni=250
baseVars = c(.2,.2)
baseCorr = 0.5
alphas = log(0.7)
betat = log(0.8)
alphaVar = 0.1
betaVar = 0.1
mstS = 3*365.25
mstT = 4*365.25 
muS <- log(log(2) / mstS)*(-1)
muT <- log(log(2) / mstT)*(-1)
d_ab <- sqrt(R2 * alphaVar * betaVar)
d_ST <- sqrt(baseCorr) * sqrt(prod(baseVars))
gammase1<-0.8
gammate1<-0.8
gammasc1<-0.8
gammatc1<-0.8

gammase2<-1.2
gammate2<-1.2
gammasc2<-0.8
gammatc2<-0.8

Sigma_trial <- matrix(c(baseVars[1], d_ST, 0, 0,
                        d_ST, baseVars[2], 0, 0,
                        0, 0, alphaVar, d_ab,
                        0, 0, d_ab, betaVar), 4)




nsim<-5000
R2c<-matrix(nrow=nsim,ncol=1)
tau<-c(5*365,10*365,15*365)
R2all<-matrix(nrow=nsim,ncol=length(tau))
followup<-matrix(nrow=nsim,ncol=3)
  
it<-0
ind<-0



for (l in 1:nsim){
  it<-it+1
  set.seed(l*2021)
  
  trialref <- unlist(mapply(rep, 1:N1, each = 250))
  data <- data.frame(trialref = factor(trialref),
                     trt = rbinom(n = length(trialref), size = 1, prob = 0.5))
  
  data$id <- factor(mapply(paste, data$trialref,
                           unlist(lapply(table(data$trialref), function(x)
                             1:x)),
                           #rep(1:ni, N),
                           sep = '.'))

  pars <- mvrnorm(n = N1,
                  mu = c(muS, muT, alphas, betat),
                  Sigma = Sigma_trial)
  rownames(pars) <- levels(data$trialref)
  
  pars <- pars[data$trialref,]
  library('msm')
  Y <- rtnorm(n = nrow(data), lower = 0)
  Y <- cbind(Y, Y)
  
  lambdaS <- rexp(n = nrow(data), rate = 1)
  lambdaT <- rexp(n = nrow(data), rate = 1)
  
  deltaS <- exp(pars[, 1] + pars[, 3] * data$trt)
  deltaT <- exp(pars[, 2] + pars[, 4] * data$trt)
  
  
  for (j in 1:nrow(data)){
    if(data$trt[j]==0){
      data$S[j] <-deltaS[j] * (Y[j, 1] * sqrt(2 * lambdaS[j])) ^ (1 / gammasc1)
      data$T[j] <-deltaT[j] * (Y[j, 1] * sqrt(2 * lambdaT[j])) ^ (1 / gammatc1)
    }
    if(data$trt[j]==1){
      data$S[j] <-deltaS[j] * (Y[j, 1] * sqrt(2 * lambdaS[j])) ^ (1 / gammase1)
      data$T[j] <-deltaT[j] * (Y[j, 1] * sqrt(2 * lambdaT[j])) ^ (1 / gammate1)
    }
  }
  
  
  
  trialref <- unlist(mapply(rep, (N1+1):(N1+N2), each = 250))
  data2 <- data.frame(trialref = factor(trialref),
                     trt = rbinom(n = length(trialref), size = 1, prob = 0.5))
  
  data2$id <- factor(mapply(paste, data2$trialref,
                           unlist(lapply(table(data2$trialref), function(x)
                             1:x)),
                           #rep(1:ni, N),
                           sep = '.'))
  
  pars <- mvrnorm(n = N2,
                  mu = c(muS, muT, alphas, betat),
                  Sigma = Sigma_trial)
  rownames(pars) <- levels(data2$trialref)
  
  pars <- pars[data2$trialref,]
  library('msm')
  Y <- rtnorm(n = nrow(data2), lower = 0)
  Y <- cbind(Y, Y)
  
  lambdaS <- rexp(n = nrow(data2), rate = 1)
  lambdaT <- rexp(n = nrow(data2), rate = 1)
  
  deltaS <- exp(pars[, 1] + pars[, 3] * data2$trt)
  deltaT <- exp(pars[, 2] + pars[, 4] * data2$trt)
  
  
  for (j in 1:nrow(data2)){
    if(data2$trt[j]==0){
      data2$S[j] <-deltaS[j] * (Y[j, 1] * sqrt(2 * lambdaS[j])) ^ (1 / gammasc2)
      data2$T[j] <-deltaT[j] * (Y[j, 1] * sqrt(2 * lambdaT[j])) ^ (1 / gammatc2)
    }
    if(data2$trt[j]==1){
      data2$S[j] <-deltaS[j] * (Y[j, 1] * sqrt(2 * lambdaS[j])) ^ (1 / gammase2)
      data2$T[j] <-deltaT[j] * (Y[j, 1] * sqrt(2 * lambdaT[j])) ^ (1 / gammate2)
    }
  }
  
  
  
 data<-rbind(data,data2) 
  
  
  
  
  #data$C <- rep(Inf, nrow(data))
  
  #random censoring
  findK <- Vectorize(function(logK) {
    (mean(runif(10 * length(data$T), 0, exp(logK)) < data$T) - 0.2) ^ 2
  })
  suppressWarnings({
    k <- exp(optim(log(max(data$T)), findK)$par)
  })
  data$C <- runif(length(data$T), 0, k)
  
  
  
  #censor at a timepoint
  data$C <- pmin(data$C, 15*365)
  
  data$timeT <- pmin(data$C, data$T)
  data$statusT <- mapply('>=', data$C, data$T) * 1
  data$timeS <- pmin(data$C, data$S)
  data$statusS <- mapply('>=', data$C, data$S) * 1
  
  for (i in 1:nrow(data)){
    if (data$timeS[i]>data$timeT[i])
      {data$statusS[i]<-0}
  }
  
  data$timeS<-pmin(data$timeS,data$timeT)
  data <- data[, !(names(data) %in% c('C', 'T', 'S'))]

  FUT<-c()
  FUS<-c()

  FU<-c()
  
  trial<-unique(data$trialref)
  
  for ( i in 1:length(unique(data$trialref))){
    #overall RMST:
    dat<-data[data$trialref==trial[i],]
    FUT[i]<- min(max(dat[which(dat$trt==0),]$timeT), 
                 max(dat[which(dat$trt==1),]$timeT))
    FUS[i]<-min(max(dat[which(dat$trt==0),]$timeS), 
                max(dat[which(dat$trt==1),]$timeS))
    FU[i]<-min(FUT[i],FUS[i])
  }  
  
  
if(min(FU)>=5*365){

  ind<-ind+1
  followup[ind,1]<-sum(FU>5*365)
  followup[ind,2]<-sum(FU>10*365)
  followup[ind,3]<-sum(FU==15*365)
  
  write.table(followup[ind,],file="resultf1n.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
  
allSurroRes <- surrosurv(data, c('clayton'), verbose = TRUE) 
 
  R2c[ind,1]<-allSurroRes$Clayton$adj$R2
  write.table(R2c[ind,],file="resultc1n.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)

 


 trial<-unique(data$trialref)
  numtau<-length(tau)
  numend<-2
  beta<-matrix(NA,ncol=numtau*numend*2,nrow=length(trial))
  covall<-list()
  
  for (i in 1:length(trial)){
    covall[[i]]<-matrix(0,ncol=numtau*numend*2,nrow=numtau*numend*2)
  }
  
  for ( i in 1:length(trial)){
    #overall RMST:
    dat<-data[data$trialref==trial[i],]
    fut<-min(max(dat[which(dat$trt==0),]$timeT), max(dat[which(dat$trt==1),]$timeT))
    fus<-min(max(dat[which(dat$trt==0),]$timeS), max(dat[which(dat$trt==1),]$timeS))
    fu<-min(fut,fus)
    tautrial<-tau[tau<=fu]
    
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
  
  
  #devtools::install_github( "cran/mvnmle")
  #library(mvnmle)
  
  #result<-mlest(beta,gradtol = 0.00001,iterlim=10000)
  
  #tau 1
  #cov<-result$sigmahat
  er<-tryCatch(mvmeta(beta, covall, method="mm"),error=function(e) 1)
  if (typeof(er)=="list"){
  REmod <- mvmeta(beta, covall, method="mm")
  cov<-REmod$Psi
  R2t<-c()
  
  R2tu<-c()
  
  
  R2tl<-c()
  for (i in 1:length(tau)){
    cov1<-cov[((i-1)*4+1):((i-1)*4+4),((i-1)*4+1):((i-1)*4+4)]
    R2<-t(c(cov1[2,3],cov1[4,3]))%*%solve(matrix(c(cov1[2,2],cov1[2,4],cov1[4,2],cov1[4,4]),nrow=2,ncol=2))%*% c(cov1[2,3],cov1[4,3])
    R2t[i]<-R2/cov1[3,3]
    N<-nrow(beta)

  }

  R2all[ind,]<-R2t
  write.table(R2all[ind,],file="resultall1n.txt",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
}
}
}


