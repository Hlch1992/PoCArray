########################################################################################################
#################################          Parameters Setup            #################################
########################################################################################################
p.range=c(0.1,0.3,0.5)
C3_C2.range=c(6,8,10)
alpha=0.05
beta=0.4

########################################################################################################
#############################          posterior positive rate            ##############################
########################################################################################################
post_p=function(n1Gobs,n1,mu,var,beta=0.4,alpha=0.05) {
  a=mu^2*((1-mu)/var-1/mu)
  b=mu^2*((1-mu)/var-1/mu)*(1/mu-1)
  integrand <- function(p) {dbeta(p,a,b)*
      dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)}
  dino=integrate(integrand, lower = 0, upper = 1)$value
  expectation=function(p) {p*dbeta(p,a,b)*
      dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)/dino}
  integrate(expectation, lower = 0, upper = 1)$value
}

########################################################################################################
################################          Strategy 1 : 1 Urn            ################################
########################################################################################################
eff_stra1=function(p,c0,C3_C2=8,n1.range=c(1:n),alpha=0.05,beta=0.4,r=(1/1.05)) {
  res=matrix(NA,5,length(n1.range),dimnames = list(c("Eff","benefit","COST","COST in Phase2","Type 3 error"),paste("n1 =",n1.range)))
  for (n1 in n1.range) {
    E_s1=n1*p*(1-beta)
    COST_s1=n1
    COST_s2=(n-n1)*(1-dbinom(0,n1,(p*(1-beta)+alpha*(1-p))))
    COST_31=n1*(p*(1-beta)+alpha*(1-p))*C3_C2
    E_s2=COST_s2*p*(1-beta)
    COST_32=COST_s2*(p*(1-beta)+alpha*(1-p))*C3_C2
    COST=ifelse (n1==n,   sum(r^c(0:4))*c0+sum(r^c(0,1))*COST_s1/2+sum(r^c(2,3,4))*COST_31/3 ,
                 sum(r^c(0:6))*c0+sum(r^c(0,1))*COST_s1/2+sum(r^c(2,3,4))*COST_31/3+sum(r^c(2,3))*COST_s2/2+sum(r^c(4,5,6))*COST_32/3 )
    res[,n1]=c(((r^2*E_s1+r^4*E_s2)/COST),(r^2*E_s1+r^4*E_s2),COST,
                   (sum(r^c(0,1))*COST_s1/2+sum(r^c(2,3))*COST_s2/2),((N*p-n1*p-COST_s2*p)/(N*p)))
  }
  t(res)
}


########################################################################################################
########################          Strategy 3_1 : 2 Urns with unequal p          ########################
########################################################################################################

eff_stra3=function(p1,p2,c0,C3_C2=8,n1.range=c(1:n),alpha=0.05,beta=0.4,I=0.05) {
  res_sce1=matrix(NA,5,length(n1.range),dimnames = list(c("Eff","benefit","COST","COST in Phase2","Type 3 error"),paste("n1 =",n1.range)))
  for (n1 in n1.range) {
    E_s1=n1*p1*(1-beta)+n1*p2*(1-beta)
    COST_s1=2*n1
    COST_s2=2*(n-n1)*(1-dbinom(0,n1,(p1*(1-beta)+alpha*(1-p1)))*dbinom(0,n1,(p2*(1-beta)+alpha*(1-p2))))
    COST_31=(n1*(p1*(1-beta)+alpha*(1-p1))+n1*(p2*(1-beta)+alpha*(1-p2)))*C3_C2
    E_s2_sce1=(1-beta)*(n-n1)*(p1+p2)*sum(dbinom(1:n1,n1,(p1*(1-beta)+alpha*(1-p1)))*dbinom(1:n1,n1,(p2*(1-beta)+alpha*(1-p2))))
    COST_3_sce1=(n-n1)*sum(dbinom(1:n1,n1,(p1*(1-beta)+alpha*(1-p1)))*dbinom(1:n1,n1,(p2*(1-beta)+alpha*(1-p2))))*((1-beta)*p1+alpha*(1-p1)+(1-beta)*p2+alpha*(1-p2))

    for (n1Gobs1 in 1:n1) {
      E_s2_sce1=E_s2_sce1+(1-beta)*(2*n-2*n1)*p1*dbinom(n1Gobs1,n1,(p1*(1-beta)+alpha*(1-p1)))*sum(dbinom(0:(n1Gobs1-1),n1,(p2*(1-beta)+alpha*(1-p2))))
      COST_3_sce1=COST_3_sce1+((1-beta)*p1+alpha*(1-p1))*(2*n-2*n1)*dbinom(n1Gobs1,n1,(p1*(1-beta)+alpha*(1-p1)))*sum(dbinom(0:(n1Gobs1-1),n1,(p2*(1-beta)+alpha*(1-p2))))
    }
    for (n1Gobs2 in 1:n1) {
      E_s2_sce1=E_s2_sce1+(1-beta)*(2*n-2*n1)*p2*dbinom(n1Gobs2,n1,(p2*(1-beta)+alpha*(1-p2)))*sum(dbinom(0:(n1Gobs2-1),n1,(p1*(1-beta)+alpha*(1-p1))))
      COST_3_sce1=COST_3_sce1+((1-beta)*p2+alpha*(1-p2))*(2*n-2*n1)*dbinom(n1Gobs2,n1,(p2*(1-beta)+alpha*(1-p2)))*sum(dbinom(0:(n1Gobs2-1),n1,(p1*(1-beta)+alpha*(1-p1))))
    }
    COST_3_sce1=COST_3_sce1*C3_C2
    if (n1==n)      COST=sum((1/(1+I))^c(0:4))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3
    else  COST=sum((1/(1+I))^c(0:6))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3+
      sum((1/(1+I))^c(2,3))*COST_s2/2+sum((1/(1+I))^c(4,5,6))*COST_3_sce1/3
    res_sce1[,n1]=c((((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2_sce1)/COST),((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2_sce1),COST,
                    (sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3))*COST_s2/2),((N*(p1+p2)-n1*(p1+p2)-E_s2_sce1/(1-beta))/(N*(p1+p2))))
  }
  t(res_sce1)
}



########################################################################################################
#######################          Strategy 3_2 : 2 Urns on bayesian rule          #######################
########################################################################################################
eff_baye_sim=function(p1,p2,p1est,var1,p2est,var2,c0=0,C3_C2=8,n1.range=c(1:n),alpha=0.05,beta=0.4,I=0.05,baye=T) {
  res_sce2=matrix(NA,5,length(n1.range),dimnames = list(c("Eff","benefit","COST","COST in Phase2","Type 3 error"),paste("n1 =",n1.range)))
  obstrail1=sample(rbinom(N,1,(p1*(1-beta)+alpha*(1-p1))))
  obstrail2=sample(rbinom(N,1,(p2*(1-beta)+alpha*(1-p2))))
  for(n1 in n1.range){
    E_s1=sum(obstrail1[1:n1]*p1*(1-beta)/(p1*(1-beta)+(1-p1)*alpha),obstrail2[1:n1]*p2*(1-beta)/(p2*(1-beta)+(1-p2)*alpha))
    COST_s1=2*n1
    if (n==n1) {
      E_s2_sce2=0;COST_s2=0;COST_31=C3_C2*sum(obstrail1[1:n1],obstrail2[1:n1]);COST_3_sce2=0
    ####### no positive trail in stage 1
    }else if (sum(obstrail1[1:n1],obstrail2[1:n1])==0){
      E_s2_sce2=0;COST_s2=0;COST_31=0;COST_3_sce2=0
    ####### drug 2 has positive trial(s) in stage 1
    }else if (sum(obstrail1[1:n1])==0) {
      E_s2_sce2=sum(obstrail2[(n1+1):(2*n-n1)])*p2*(1-beta)/(p2*(1-beta)+(1-p2)*alpha);COST_s2=2*(n-n1);
      COST_31=C3_C2*sum(obstrail2[1:n1]);COST_3_sce2=C3_C2*sum(obstrail2[(n1+1):(2*n-n1)])
    ####### drug 1 has positive trial(s) in stage 1
    }else if (sum(obstrail2[1:n1])==0) {
      E_s2_sce2=sum(obstrail1[(n1+1):(2*n-n1)])*p2*(1-beta)/(p2*(1-beta)+(1-p2)*alpha);COST_s2=2*(n-n1);
      COST_31=C3_C2*sum(obstrail1[1:n1]);COST_3_sce2=C3_C2*sum(obstrail2[(n1+1):(2*n-n1)])
    ####### drug 1 and 2 have positive trials in stage 1
    }else {
      if (baye) {
      post.p1=post_p(sum(obstrail1[1:n1]),n1,p1est,var1)
      post.p2=post_p(sum(obstrail2[1:n1]),n1,p2est,var2)
      } else { post.p1=sum(obstrail1[1:n1]);    post.p2=sum(obstrail2[1:n1]) }

      ####### drug 1 and 2 have equally positive posterier p in stage 1
      if (post.p1==post.p2) {
        E_s2_sce2=sum(obstrail1[(n1+1):n]*p1*(1-beta)/(p1*(1-beta)+(1-p1)*alpha),obstrail2[(n1+1):n]*p2*(1-beta)/(p2*(1-beta)+(1-p2)*alpha));COST_s2=2*(n-n1);
        COST_31=C3_C2*sum(obstrail1[1:n1],obstrail2[1:n1]);COST_3_sce2=C3_C2*sum(obstrail1[(n1+1):n]+obstrail2[(n1+1):n])
      ####### drug 2 has more positive posterier p after stage 1
      }else if (post.p1<post.p2) {
        E_s2_sce2=sum(obstrail2[(n1+1):(2*n-n1)])*p2*(1-beta)/(p2*(1-beta)+(1-p2)*alpha);COST_s2=2*(n-n1);
        COST_31=C3_C2*sum(obstrail1[1:n1],obstrail2[1:n1]);COST_3_sce2=C3_C2*sum(obstrail2[(n1+1):(2*n-n1)])
      ####### drug 1 has more positive posterier p sfter stage 1
      }else {
        E_s2_sce2=sum(obstrail1[(n1+1):(2*n-n1)])*p1*(1-beta)/(p1*(1-beta)+(1-p1)*alpha);COST_s2=2*(n-n1);
        COST_31=C3_C2*sum(obstrail1[1:n1],obstrail2[1:n1]);COST_3_sce2=C3_C2*sum(obstrail1[(n1+1):(2*n-n1)])
      }
    }
    if (n1==n)      COST=sum((1/(1+I))^c(0:4))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3
    else  COST=sum((1/(1+I))^c(0:6))*c0+sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3,4))*COST_31/3+
      sum((1/(1+I))^c(2,3))*COST_s2/2+sum((1/(1+I))^c(4,5,6))*COST_3_sce2/3
    res_sce2[,n1]=c((((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2_sce2)/COST),((1/(1+I))^2*E_s1+(1/(1+I))^4*E_s2_sce2),COST,
                    (sum((1/(1+I))^c(0,1))*COST_s1/2+sum((1/(1+I))^c(2,3))*COST_s2/2),((N*(p1+p2)-E_s1/(1-beta)-E_s2_sce2/(1-beta))/(N*(p1+p2))))

  }
  t(res_sce2)
}


eff_baye=function(p1,p2,p1est,var1,p2est,var2,c0=0,C3_C2=8,n1.range=c(1:n),alpha=0.05,beta=0.4,I=0.05,baye=T,nsim=10000) {
  set.seed(1)
  res=eff_baye_sim(p1,p2,p1est,var1,p2est,var2,c0=c0,C3_C2=C3_C2,n1.range=n1.range,alpha=alpha,beta=beta,I=I,baye=baye)
  for (t in 2:nsim) {
    set.seed(t)
    res=res+eff_baye_sim(p1,p2,p1est,var1,p2est,var2,c0=c0,C3_C2=C3_C2,n1.range=n1.range,alpha=alpha,beta=beta,I=I,baye=baye)
  }
  res=cbind(res,res[,1]/nsim)
  res[,1]=res[,2]/res[,3]
  res[,2:5]=res[,2:5]/nsim
  res
}

########################################################################################################
################################          functions for plots            ###############################
########################################################################################################
### Sensitivities of priors
proplot=function(mu,var,beta=0.4,alpha=0.05) {
  n1=10
  col=rainbow(11)
  a=mu^2*((1-mu)/var-1/mu)
  b=a*(1/mu-1)
  prom=matrix(NA,11,101)
  for (n1Gobs in 0:n1) {
    integrand <- function(p) {dbeta(p,a,b)*dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)}
    expectation=function(p) {dbeta(p,a,b)*dbinom(n1Gobs,n1,(1-beta)*p+(1-p)*alpha)/integrate(integrand, lower = 0, upper = 1)$value}
    prom[(n1Gobs+1),]=sapply(seq(0,1,0.01), expectation)
  }
  proplot=plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),a,b),ylim = range(prom[,-c(1:10)]),type = 'l',main = paste("Probability of positive: mu =",mu,"var =",var),xlab = 'p',ylab='probability')
  for (n1Gobs in 0:n1) {
    proplot=lines(seq(0,1,0.01),prom[(n1Gobs+1),],col=col[n1Gobs+1])
  }
}

#########   heatmaps (frequentist rule)
hmplot=function(p1,simplify,var=0.01,baye=F,rev=F) {
  max.m=0
  p2.range=seq(0.1,0.5,0.01)
  c0.range=seq(0,0.4,0.01)
  colors=rainbow(11)
  weightf=ifelse(simplify,1.05,1)
  for (p2 in p2.range) {
    for (c0 in c0.range) {
      eff1=eff_stra1(p1,c0)
      eff2=eff_stra1(p2,c0)
      eff3=ifelse(matrix(baye,10,5),ifelse(matrix(rev,10,5),eff_baye(p1,p2,p2,var,p1,var,c0,nsim=nsim),eff_baye(p1,p2,p1,var,p2,var,c0,nsim=nsim)),eff_stra3(p1,p2,c0))
      res=c((((eff1[1,1]*eff1[1,3]+eff2[1,1]*eff2[1,3])/(eff1[1,3]+eff2[1,3]))*weightf),eff3[1,1]*weightf,eff3[2:9,1],eff3[10,1]*weightf)
      max.m=c(max.m,which.max(res))
    }
  }
  figure=matrix(max.m[-1],41,41,dimnames = list(value_p1=p2.range,value_p2=p2.range))
  colormatrix=colors[figure]
  hmp=plot(rep(c0.range,41),rep(p2.range,each=41),col=colormatrix,xlab = "C0", ylab = "value_p2",main = paste("p1 =",p1,'var =', var,'simplify =',simplify))
  hmp=legend("topright",c("con1","con2","2","3","4","5","6","7","8","9","agg"),col = colors,pch = 15)
}

#########   heatmaps (bayesian rule)
bayehmplot=function(p1,var,simplify) {
  max.m=0
  p2.range=seq(0.1,0.5,0.01)
  c0.range=seq(0,0.4,0.01)
  colors=rainbow(6)
  weightf=ifelse(simplify,1.05,1)
  for (p2 in p2.range) {
    for (c0 in c0.range) {
      eff1=eff_stra1(p1,c0)
      eff2=eff_stra1(p2,c0)
      eff3=eff_baye(p1,p2,p1,var,p2,var,c0,nsim = nsim)
      res=c((((eff1[1,1]*eff1[1,3]+eff2[1,1]*eff2[1,3])/(eff1[1,3]+eff2[1,3]))*weightf),eff3[1,1]*weightf,eff3[2:4,1],eff3[5,1]*weightf)
      max.m=c(max.m,which.max(res))
    }
  }
  figure=matrix(max.m[-1],41,41,dimnames = list(value_p1=p2.range,value_p2=p2.range))
  colormatrix=colors[figure]
  hmp=plot(rep(c0.range,41),rep(p2.range,each=41),col=colormatrix,xlab = "C0", ylab = "value_p2",main = paste("p1 =",p1,'simplify =',simplify))
  hmp=legend("topright",c("con1","con2","2","3","4","agg"),col = colors,pch = 15)
  hmp
}

#########   plot for comparson on bayesian rule
complot=function(p1,p2,c0,explore) {
  colors=rainbow(10)
  switch (explore,efficiency = {column=1},benefit={column=2},cost={column=3},type3error={column=5})
  upper=eff_baye(p1,p2,p1,0.004,p2,0.004,c0,nsim=nsim)[,column]
  lower=eff_baye(p1,p2,p2,0.004,p1,0.004,c0,nsim=nsim)[,column]
  plot(seq(1,n),upper,col=colors[1],ylab = explore,xlab = "n1",ylim = range(c(upper,lower)),main = paste("mu1 =",p1,"mu2 =",p2,'c0 =',c0),type = "l")
  lines(seq(1,n),eff_baye(p1,p2,p1,0.01,p1,0.01,c0,nsim=nsim)[,column],type = 'l',col=colors[4])
  lines(seq(1,n),lower,type = 'l',col=colors[7])
  lines(seq(1,n),eff_baye(p1,p2,p1,0.04,p2,0.04,c0,nsim=nsim)[,column],type = 'l',col=colors[2])
  lines(seq(1,n),eff_baye(p1,p2,p1,0.04,p1,0.04,c0,nsim=nsim)[,column],type = 'l',col=colors[5])
  lines(seq(1,n),eff_baye(p1,p2,p2,0.04,p1,0.04,c0,nsim=nsim)[,column],type = 'l',col=colors[8])
  #  lines(seq(1,n),eff_baye(p1,0.004,p2,0.004,p1,p2,c0)[,column],type = 'l',col=colors[3])
  #  lines(seq(1,n),eff_baye(p1,0.004,p2,0.004,p1,p1,c0)[,column],type = 'l',col=colors[6])
  #  lines(seq(1,n),eff_baye(p1,0.004,p2,0.004,p2,p1,c0)[,column],type = 'l',col=colors[9])
  lines(seq(1,n),eff_stra3(p1,p2,c0)[,column],type = 'l',col='black')
  #  if (explore=='efficiency')  legend(ifelse(c0==0 | (p2-p1)==0.4,'topright','bottomright'),c('true prior','same prior','wrong prior','flat prior'),col = c(colors[c(2,5,8)],'black'),lty = 1)
  #  else if (explore=='type3error') legend('topright',c('true prior','same prior','wrong prior','flat prior'),col = c(colors[c(2,5,8)],'black'),lty = 1)
  #  else  legend('bottomright',c('true prior','same prior','wrong prior','flat prior'),col = c(colors[c(2,5,8)],'black'),lty = 1)
}

