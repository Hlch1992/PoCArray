source('functions.R')

### few indications
n=5
N=15
### many indications
n=10
N=30




### frequntist rule without resource reallocation
scenario=matrix(c(rep(c(0.1,0.3,0.5),each=5),rep(seq(0,0.4,0.1),3)),ncol=2)
x=apply(scenario,1,function(x) eff_stra1(x[1],c0=x[2]))
### scenario 2
#y=apply(scenario,1,function(x) eff_stra2(x[1],c0=x[2]))


### frequntist rule with resource reallocation
scenarioz=matrix(c(rep(0.1,15),rep(0.3,15),rep(0.5,15),rep(c(rep(0.1,5),rep(0.3,5),rep(0.5,5)),3),rep(seq(0,0.4,0.1),9)),ncol=3)
z=apply(scenarioz,1,function(x) eff_stra3(x[1],x[2],c0=x[3]))


### frequent rule summary
scenario3=matrix(c(rep(0.1,9),rep(0.3,9),rep(0.5,9),rep(c(rep(0.1,3),rep(0.3,3),rep(0.5,3)),3),rep(seq(0,0.4,0.2),9)),ncol=3)
scenario3list= split(scenario3, rep(1:nrow(scenario3), ncol(scenario3)))
scenario3list= split(scenario3, rep(1:nrow(scenario3), ncol(scenario3)))
w=lapply(scenario3list, function(x) eff_stra3(x[1],x[2],c0=x[3]))
### ratio outcome over aggresive strategy
toagg=lapply(w, function(x) t(t(x)/x[n,]))
efficiency=matrix(unlist(lapply(toagg, function(x) x[,1])),ncol=n,byrow = T)
benefit=matrix(unlist(lapply(toagg, function(x) x[,2])),ncol=n,byrow = T)
cost=matrix(unlist(lapply(toagg, function(x) x[,3])),ncol=n,byrow = T)
type3error=matrix(unlist(lapply(toagg, function(x) x[,5])),ncol=n,byrow = T)
write.csv(round(efficiency,4),'efficiency.csv')
write.csv(round(benefit,4),'benefit.csv')
write.csv(round(cost,4),'cost.csv')
write.csv(round(type3error,4),'type3error.csv')


### comparson of bayesian rule to frequent rule
### testing: small number of simulation
nsim=10#00

prior=lapply(scenario3list, function(x)
  array(c(eff_baye(x[1],x[2],x[1],0.04,x[2],0.04,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[1],0.04,x[1],0.04,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[2],0.04,x[1],0.04,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[1],0.01,x[2],0.01,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[1],0.01,x[1],0.01,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[2],0.01,x[1],0.01,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[1],0.004,x[2],0.004,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[1],0.004,x[1],0.004,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3]),
          eff_baye(x[1],x[2],x[2],0.004,x[1],0.004,c0=x[3],nsim=nsim)[,1:5]/eff_stra3(x[1],x[2],c0=x[3])),
        dim = c(n,5,9)))
priorform=lapply(prior, function(x) cbind(rbind(x[,,1],x[,,4],x[,,7]),rbind(x[,,2],x[,,5],x[,,8]),rbind(x[,,3],x[,,6],x[,,9])))
priorformreduce=cbind(rbind(priorform[[4]],priorform[[5]],priorform[[6]]),rbind(priorform[[7]],priorform[[8]],priorform[[9]]),
                            rbind(priorform[[16]],priorform[[17]],priorform[[18]]))
b=c(2,7,12,17,22,27,32,37,42)
write.csv(round(priorformreduce[,c(b-1)],4),'efficiencybaye.csv')
write.csv(round(priorformreduce[,b],4),'benefitbaye.csv')
write.csv(round(priorformreduce[,c(b+1)],4),'costbaye.csv')
write.csv(round(priorformreduce[,c(b+3)],4),'type3errorbaye.csv')



### Sensitivities of priors
proplot(0.5,0.01)
proplot(0.5,0.04)
proplot(0.5,0.004)
proplot(0.3,0.01)
proplot(0.3,0.04)
proplot(0.3,0.004)
proplot(0.1,0.01)
proplot(0.1,0.04)
proplot(0.1,0.004)
legend('topright',c('prior','0 outof 10','1 outof 10','2 outof 10','3 outof 5','4 outof 5','5 outof 5'),col=c('black',col),lty=1)



plot(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.5,0.01,c0=0)[,1],ylab = "efficiency",xlab = "n1",ylim = c(0.064,0.079),main = "prior mu1 = prior mu2, var=0.01",type = "l")
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.5,0.01,c0=0.1)[,1],type = 'l',col='blue')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.5,0.01,c0=0.2)[,1],type = 'l',col='red')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.5,0.01,c0=0.3)[,1],type = 'l',col='yellow')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.5,0.01,c0=0.4)[,1],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.1','c0=0.2','c0=0.3','c0=0.4'),col=c('black','blue','red','yellow','green'),lty =1)

plot(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.3,0.01,c0=0)[,1],ylab = "efficiency",xlab = "n1",ylim = c(0.064,0.079),main = "prior mu1 = prior mu2, var=0.01",type = "l")
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.3,0.01,c0=0.1)[,1],type = 'l',col='blue')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.3,0.01,c0=0.2)[,1],type = 'l',col='red')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.3,0.01,c0=0.3)[,1],type = 'l',col='yellow')
lines(seq(1,n),eff_baye(0.3,0.5,0.3,0.01,0.3,0.01,c0=0.4)[,1],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.1','c0=0.2','c0=0.3','c0=0.4'),col=c('black','blue','red','yellow','green'),lty =1)

plot(seq(1,n),eff_baye(0.3,0.5,0.5,0.01,0.3,0.01,c0=0)[,1],ylab = "efficiency",xlab = "n1",ylim = c(0.064,0.079),main = "reverse truth prior, var=0.01",type = "l")
lines(seq(1,n),eff_baye(0.3,0.5,0.5,0.01,0.3,0.01,c0=0.1)[,1],type = 'l',col='blue')
lines(seq(1,n),eff_baye(0.3,0.5,0.5,0.01,0.3,0.01,c0=0.2)[,1],type = 'l',col='red')
lines(seq(1,n),eff_baye(0.3,0.5,0.5,0.01,0.3,0.01,c0=0.3)[,1],type = 'l',col='yellow')
lines(seq(1,n),eff_baye(0.3,0.5,0.5,0.01,0.3,0.01,c0=0.4)[,1],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.1','c0=0.2','c0=0.3','c0=0.4'),col=c('black','blue','red','yellow','green'),lty =1)











plot(seq(1,n),z[1:n,11],type = 'l',ylim = c(0.063,0.076),ylab = 'efficency',xlab = 'n1',main = "p1=0.1,p2=0.5")
#lines(seq(1,n),z[1:n,12],type = 'l',col='blue')
lines(seq(1,n),z[1:n,13],type = 'l',col='red')
#lines(seq(1,n),z[1:n,14],type = 'l',col='yellow')
lines(seq(1,n),z[1:n,15],type = 'l',col='green')
legend('topright',c('c0=0','c0=0.1','c0=0.2','c0=0.3','c0=0.4'),col=c('black','blue','red','yellow','green'),lty =1)


plot(seq(1,n),z[1:n,26],type = 'l',ylim = c(0.069,0.079),ylab = 'efficency',xlab = 'n1',main = "p1=0.3,p2=0.5")
#lines(seq(1,n),z[1:n,27],type = 'l',col='blue')
lines(seq(1,n),z[1:n,28],type = 'l',col='red')
#lines(seq(1,n),z[1:n,29],type = 'l',col='yellow')
lines(seq(1,n),z[1:n,30],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)

n=5
N=15
z1=apply(scenarioz,1,function(x) eff_stra3(x[1],x[2],c0=x[3]))
c1=eff_stra1(.1,0)+eff_stra1(.5,0)
c2=eff_stra1(.1,0.2)+eff_stra1(.5,0.2)
c3=eff_stra1(.1,0.4)+eff_stra1(.5,0.4)
plot(seq(1,n),c1[,2]/c1[,3],type = 'l',ylim = c(0.039,0.074),ylab = 'efficency',xlab = 'n1',main = "p1=0.1,p2=0.5")
lines(seq(1,n),c2[,2]/c2[,3],type = 'l',col='red')
lines(seq(1,n),c3[,2]/c3[,3],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)

plot(seq(1,n),z1[1:n,11],type = 'l',ylim = c(0.049,0.074),ylab = 'efficency',xlab = 'n1',main = "p1=0.1,p2=0.5")
#lines(seq(1,n),z[1:n,27],type = 'l',col='blue')
lines(seq(1,n),z1[1:n,13],type = 'l',col='red')
#lines(seq(1,n),z[1:n,29],type = 'l',col='yellow')
lines(seq(1,n),z1[1:n,15],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)



colors=rainbow(10)
plot(seq(1,n),eff_baye(.1,.5,.1,.04,.5,.04,c0=0)[,1],type = 'l',col=colors[2],
     ylab = "efficiency",xlab = "n1",ylim = c(0.064,0.074),main = "few c0=0")
lines(seq(1,n),eff_baye(0.1,0.5,0.1,0.01,0.5,0.01,c0=0)[,1],type = 'l',col=colors[1])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.04,0.1,0.04,c0=0)[,1],type = 'l',col=colors[8])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.01,0.1,0.01,c0=0)[,1],type = 'l',col=colors[7])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0)[,1],type = 'l',col=colors[5])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0,baye=F)[,1],type = 'l',col='black')
legend('top',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],colors[c(5)],'black'),lty = 1)#,'white',4,3,6,9


plot(seq(1,n),eff_baye(.1,.5,.1,.04,.5,.04,c0=0.4)[,1],type = 'l',col=colors[2],
     ylab = "efficiency",xlab = "n1",ylim = c(0.058,0.065),main = "few c0=0.4")
lines(seq(1,n),eff_baye(0.1,0.5,0.1,0.01,0.5,0.01,c0=0.4)[,1],type = 'l',col=colors[1])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.04,0.1,0.04,c0=0.4)[,1],type = 'l',col=colors[8])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.01,0.1,0.01,c0=0.4)[,1],type = 'l',col=colors[7])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0.4)[,1],type = 'l',col=colors[5])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0.4,baye=F)[,1],type = 'l',col='black')
legend('top',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],colors[c(5)],'black'),lty = 1)#,'white',4,3,6,9


n=10
N=30
z=apply(scenarioz,1,function(x) eff_stra3(x[1],x[2],c0=x[3]))
c1=eff_stra1(.1,0)+eff_stra1(.5,0)
c2=eff_stra1(.1,0.2)+eff_stra1(.5,0.2)
c3=eff_stra1(.1,0.4)+eff_stra1(.5,0.4)
plot(seq(1,n),c1[,2]/c1[,3],type = 'l',ylim = c(0.049,0.074),ylab = 'efficency',xlab = 'n1',main = "p1=0.1,p2=0.5")
lines(seq(1,n),c2[,2]/c2[,3],type = 'l',col='red')
lines(seq(1,n),c3[,2]/c3[,3],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)

plot(seq(1,n),z[1:n,11],type = 'l',ylim = c(0.055,0.075),ylab = 'efficency',xlab = 'n1',main = "p1=0.1,p2=0.5")
#lines(seq(1,n),z[1:n,27],type = 'l',col='blue')
lines(seq(1,n),z[1:n,13],type = 'l',col='red')
#lines(seq(1,n),z[1:n,29],type = 'l',col='yellow')
lines(seq(1,n),z[1:n,15],type = 'l',col='green')
legend('bottomright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)




colors=rainbow(10)
plot(seq(1,n),eff_baye(.1,.5,.1,.04,.5,.04,c0=0)[,1],type = 'l',col=colors[2],
     ylab = "efficiency",xlab = "n1",ylim = c(0.064,0.074),main = "few c0=0")
lines(seq(1,n),eff_baye(0.1,0.5,0.1,0.01,0.5,0.01,c0=0)[,1],type = 'l',col=colors[1])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.04,0.1,0.04,c0=0)[,1],type = 'l',col=colors[8])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.01,0.1,0.01,c0=0)[,1],type = 'l',col=colors[7])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0)[,1],type = 'l',col=colors[5])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0,baye=F)[,1],type = 'l',col='black')
legend('top',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],colors[c(5)],'black'),lty = 1)#,'white',4,3,6,9


plot(seq(1,n),eff_baye(.1,.5,.1,.04,.5,.04,c0=0.4)[,1],type = 'l',col=colors[2],
     ylab = "efficiency",xlab = "n1",ylim = c(0.058,0.065),main = "few c0=0.4")
lines(seq(1,n),eff_baye(0.1,0.5,0.1,0.01,0.5,0.01,c0=0.4)[,1],type = 'l',col=colors[1])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.04,0.1,0.04,c0=0.4)[,1],type = 'l',col=colors[8])
lines(seq(1,n),eff_baye(0.1,0.5,0.5,0.01,0.1,0.01,c0=0.4)[,1],type = 'l',col=colors[7])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0.4)[,1],type = 'l',col=colors[5])
lines(seq(1,n),eff_baye(0.1,0.5,0.3,0.04,0.3,0.04,c0=0.4,baye=F)[,1],type = 'l',col='black')
legend('top',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],colors[c(5)],'black'),lty = 1)#,'white',4,3,6,9



write.csv(z1,'fplot12.csv')
write.csv(z,'fplot1212.csv')





#########   heatmaps
par(mfrow=c(1,1))
#hmplot(0.5,F)
hmplot(0.3,F)
#hmplot(0.1,F)
#hmplot(0.5,T)
hmplot(0.3,T)
#hmplot(0.1,T)



par(mfrow=c(1,1))
hmplot(0.3,T,0.01,baye=T)
hmplot(0.3,T,0.01,baye=T,rev = T)

#########   heatmaps (bayesian rule)
bayehmplot(0.3,.004,T)
bayehmplot(0.3,.004,F)




###### plot legends
par(mfrow=c(1,1))
plot(c(1:3),c(1:3))
colors=rainbow(11)
hmp=legend("topright",c("con1","con2","2","3","4","5","6","7","8","9","agg"),col = colors,pch = 15)


par(mfrow=c(1,1))
plot(c(1:3),c(1:3))
colors=rainbow(10)
legend('top',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                  colors[c(8,7)],colors[c(5)],'black'),lty = 1)#,'white',4,3,6,9
legend('top',c("bayesian rule with true prior",'variance = 0.04','variance = 0.01','variance = 0.004',
               "bayesian rule with same prior",'variance = 0.04','variance = 0.01','variance = 0.004',
               "bayesian rule with wrong prior",'variance = 0.04','variance = 0.01','variance = 0.004',
               "frequent rule"),col=c('white',colors[c(2,1,3)],'white',colors[c(5,4,6)],'white',
                                  colors[c(8,7,9)],'black'),lty = 1)

par(mfrow=c(1,1))
plot(c(1:3),c(1:3))
colors=rainbow(11)
legend('topright',c('prior','0 outof 10','1 outof 10','2 outof 10','3 outof 10','4 outof 10','5 outof 10','6 outof 10','7 outof 10','8 outof 10','9 outof 10','10 outof 10'),col=c('black',col),lty=1)
