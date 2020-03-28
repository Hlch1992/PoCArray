########################################################################################################
#############################          posterior positive rate            ##############################
########################################################################################################
post_p=function(n1Gobs,n1,mu,var,beta=0.4,alpha=0.05) {
  #### assuming a beta distribution for probability of sucess 
  #### calculate the parameter in beta distribution given the prior mean and variance
  a=mu^2*((1-mu)/var-1/mu)
  b=mu^2*((1-mu)/var-1/mu)*(1/mu-1)
  integrand <- function(p) {dbeta(p,a,b)*
      dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)}
  dino=integrate(integrand, lower = 0, upper = 1)$value
  #### calculate the posterior p value given the observed trials
  expectation=function(p) {p*dbeta(p,a,b)*
      dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)/dino}
  integrate(expectation, lower = 0, upper = 1)$value
}

########################################################################################################
#########          Strategy multiurn : Simulation with frequentist and bayesian rule           #########
########################################################################################################
multiurnsim=function(p=c(.1,.3,.3,.5),i=4,c0=0,n1.range=1:5,C3_C2=8,n=5,N=15,alpha=.05,beta=.4,I=0.05,baye=F,pest=NA,pvar=NA) {
  if (length(p)!=i) {print("number of urns not match");break}
  if (max(n1.range)>n) {print("more than n trials in stage 1");break}
  if (n>N) {print("more than N trials in POC array");break}
  if (baye & (is.na(pest) | is.na(pvar))) {print("prior information missing");break}
  if (baye & length(pest)!=length(p)) {print("number of prior not match");break}
  
  #### simulate all the trials
  obsout=t(sapply(p,function(x) sample(rbinom(N,1,(x*(1-beta)+(1-x)*alpha)))))
  res=matrix(NA,5,length(n1.range),dimnames = list(c("Eff","benefit","COST","COST in Phase2","Type 3 error"),paste("n1 =",n1.range)))
  
  #### investigate the number of trials in stage 1
  for(n1 in n1.range){
    
    #### the sucess trial number in stage 1
    if(n1==1) {sucessn1=sucessn1copy=obsout[,1]
    }else sucessn1=sucessn1copy=rowSums(obsout[,1:n1])
    
    #### calculate the posterior p under a bayesian rule
    if(baye)     post.p=sapply(1:i,function(x) post_p(sucessn1[x],n1,pest[x],pvar[x]))
    
    #### benefit in stage 1
    E_s1=sum(sucessn1*p*(1-beta)/(p*(1-beta)+(1-p)*alpha))
    
    #### cost in stage 1
    COST_s1=i*n1
    
    #### cost in stage 2
    COST_s2=min(i*(n-n1),sum(sucessn1>0)*(N-n1))
    
    #### cost in phase 3 from stage 1
    COST_31=sum(sucessn1)*C3_C2
    
    #### output 
    if (n1==n)      {
      E_s2=0
      COST_32=0
      COST=sum((1/(1+I))^c(0:4))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3
    }else  {
      
      #### trials for stage 2 given the observations in stage 1
      forstage2=c()
      pstage2=c()
      
      if(baye) {
        ####  a bayesian rule
        while(sum(sucessn1copy*post.p)>0) {
          forstage2=c(forstage2,as.vector(obsout[which(sucessn1copy*post.p==max(sucessn1copy*post.p)),-c(1:n1)]))
          pstage2=c(pstage2,rep(p[which(sucessn1copy*post.p==max(sucessn1copy*post.p))],(N-n1)))
          sucessn1copy[which(sucessn1copy*post.p==max(sucessn1copy*post.p))]=0  }
      } else {
        
        #### a frequentist rule
        while(sum(sucessn1copy)>0){
          forstage2=c(forstage2,as.vector(obsout[which(sucessn1copy==max(sucessn1copy)),-c(1:n1)]))
          pstage2=c(pstage2,rep(p[which(sucessn1copy==max(sucessn1copy))],(N-n1)))
          sucessn1copy[which(sucessn1copy==max(sucessn1copy))]=0  }
      }
      forstage2=forstage2[1:min(i*(n-n1),sum(sucessn1>0)*(N-n1))]
      pstage2=pstage2[1:min(i*(n-n1),sum(sucessn1>0)*(N-n1))]
      
      #### benefit in stage 2
      E_s2=sum(forstage2*pstage2*(1-beta)/(pstage2*(1-beta)+(1-pstage2)*alpha))
      
      #### cost in phase 3 from stage 2
      COST_32=sum(forstage2)*C3_C2
      
      COST=sum((1/(1+I))^c(0:6))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3+
        sum((1/(1+I))^c(2,3))*COST_s2/2+sum((1/(1+I))^c(4,5,6))*COST_32/3
    }
    res[,n1]=c((((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2)/COST),((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2),COST,
                    (sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3))*COST_s2/2),((N*sum(p)*(1-beta)-E_s1-E_s2)/(N*sum(p)*(1-beta))))
  }
  t(res)
}

eff_multiurn=function(p=c(.1,.3,.3,.5),i=4,c0=0,n1.range=1:5,C3_C2=8,n=5,N=15,alpha=.05,beta=.4,I=0.05,baye=F,pest=NA,pvar=NA,nsim=10000) {
  set.seed(1)
  res=multiurnsim(p=p,i=i,c0=c0,C3_C2=C3_C2,n1.range=n1.range,n=n,N=N,alpha=alpha,beta=beta,I=I,baye=baye,pest=pest,pvar=pvar)
  for (t in 2:nsim) {
    set.seed(t)
    res=res+multiurnsim(p=p,i=i,c0=c0,C3_C2=C3_C2,n1.range=n1.range,n=n,N=N,alpha=alpha,beta=beta,I=I,baye=baye,pest=pest,pvar=pvar)
  }
  res=cbind(res,res[,1]/nsim)
  res[,1]=res[,2]/res[,3]
  res[,2:5]=res[,2:5]/nsim
  res
}

#########   heatmaps (frequentist rule)
hmplot=function(simplify=F,var=0.01,baye=F,nsim=1000,pfix=c(.3,.4,.5),prior=T) {
  max.m=0
  p.range=seq(0.1,0.5,0.01)
  c0.range=seq(0,0.4,0.01)
  colors=rainbow(5)
  weightf=ifelse(simplify,1.05,1)
  for (p in p.range) {
    for (c0 in c0.range) {
      if (prior)      {eff=eff_multiurn(p=c(pfix,p),pvar = rep(var,4),pest = c(pfix,p),c0=c0,baye = baye,nsim = nsim)
      } else eff=eff_multiurn(p=c(pfix,p),pvar = rep(var,4),pest = (.6- c(pfix,p)),c0=c0,baye = baye,nsim = nsim)
      res=c(eff[1,1]*weightf,eff[2:4,1],eff[5,1]*weightf)
      max.m=c(max.m,which.max(res))
    }
  }
  figure=matrix(max.m[-1],41,41,dimnames = list(value_p1=p.range,value_p2=p.range))
  colormatrix=colors[figure]
  hmp=plot(rep(c0.range,41),rep(p.range,each=41),col=colormatrix,xlab = "C0", ylab = "value_p",
           xlim = c(0,.45),main = paste("p = (",pfix[1],pfix[2],pfix[3],'); simplify =',simplify))
  hmp=legend("topright",c("con","2","3","4","agg"),col = colors,pch = 15)
}
