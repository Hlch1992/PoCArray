source('multiurn.R')

pfix=c(.1,.3,.5)
pfix=c(.3,.4,.5)
### frequntist rule with resource reallocation
scenario=matrix(c(rep(c(0.1,0.3,0.5),each=3),rep(seq(0,0.4,0.2),3)),ncol=2)
x=apply(scenario,1,function(x) eff_multiurn(p=c(pfix,x[1]),c0=x[2],nsim=10))
write.csv(x,paste0(pfix[1],pfix[2],pfix[3],'frequentistrule.csv'))
x=read.csv(paste0(pfix[1],pfix[2],pfix[3],'frequentistrule.csv'))[,-1]

### bayesian rule with true prior
y=apply(scenario,1,function(x) eff_multiurn(p=c(pfix,x[1]),c0=x[2],pest = c(pfix,x[1]),pvar = rep(.01,4),baye = T,nsim=10))
write.csv(y,paste0(pfix[1],pfix[2],pfix[3],'bayesianruletrueprior.01.csv'))
y=read.csv(paste0(pfix[1],pfix[2],pfix[3],'bayesianruletrueprior.01.csv'))[,-1]

### bayesian rule with wrong prior
z=apply(scenario,1,function(x) eff_multiurn(p=c(pfix,x[1]),c0=x[2],pest = (.6-c(pfix,x[1])),pvar = rep(.01,4),baye = T,nsim=10))
write.csv(z,paste0(pfix[1],pfix[2],pfix[3],'bayesianrulewrongprior.01.csv'))
z=read.csv(paste0(pfix[1],pfix[2],pfix[3],'bayesianrulewrongprior.01.csv'))[,-1]

### bayesian rule with true prior
w=apply(scenario,1,function(x) eff_multiurn(p=c(pfix,x[1]),c0=x[2],pest = c(pfix,x[1]),pvar = rep(.04,4),baye = T,nsim=10))
write.csv(w,paste0(pfix[1],pfix[2],pfix[3],'bayesianruletrueprior.04.csv'))
w=read.csv(paste0(pfix[1],pfix[2],pfix[3],'bayesianruletrueprior.04.csv'))[,-1]

### bayesian rule with wrong prior
v=apply(scenario,1,function(x) eff_multiurn(p=c(pfix,x[1]),c0=x[2],pest = (.6-c(pfix,x[1])),pvar = rep(.04,4),baye = T,nsim=10))
write.csv(v,paste0(pfix[1],pfix[2],pfix[3],'bayesianrulewrongprior.04.csv'))
v=read.csv(paste0(pfix[1],pfix[2],pfix[3],'bayesianrulewrongprior.04.csv'))[,-1]

### frequent rule summary

#### p = .1 .3 .5 .1, c0 = 0 .2 .4
plot(seq(1:5),x[1:5,1],ylim = range(x[1:5,1:3]),xlab = 'n1',main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.1)'),ylab = 'efficency',
     xlim = c(1,6),col='black',type = 'l')
lines(seq(1:5),x[1:5,2],type = 'l',col='red')
lines(seq(1:5),x[1:5,3],type = 'l',col='green')
legend('topright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)

#### p = .1 .3 .5 .3, c0 = 0 .2 .4
plot(seq(1:5),x[1:5,4],ylim = range(x[1:5,4:6]),xlab = 'n1',main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.3)'),ylab = 'efficency',
     xlim = c(1,6),col='black',type = 'l')
lines(seq(1:5),x[1:5,5],type = 'l',col='red')
lines(seq(1:5),x[1:5,6],type = 'l',col='green')
legend('topright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)

#### p = .1 .3 .5 .5, c0 = 0 .2 .4
plot(seq(1:5),x[1:5,7],ylim = range(x[1:5,7:9]),xlab = 'n1',main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.5)'),ylab = 'efficency',
     xlim = c(1,6),col='black',type = 'l')
lines(seq(1:5),x[1:5,8],type = 'l',col='red')
lines(seq(1:5),x[1:5,9],type = 'l',col='green')
legend('topright',c('c0=0','c0=0.2','c0=0.4'),col=c('black','red','green'),lty =1)



colors=rainbow(10)
plot(seq(1,n),w[1:5,1],type = 'l',col=colors[2], ylab = "efficiency",xlab = "n1",xlim = c(1,7),
     ylim = range(c(w[1:5,1],v[1:5,1],y[1:5,1],z[1:5,1])),main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.1)'))
lines(seq(1,n),y[1:5,1],type = 'l',col=colors[1])
lines(seq(1,n),v[1:5,1],type = 'l',col=colors[8])
lines(seq(1,n),z[1:5,1],type = 'l',col=colors[7])
lines(seq(1,n),x[1:5,1],type = 'l',col='black')
legend('topright',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',#"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],'black'),lty = 1)#,colors[c(5)],'white',4,3,6,9


plot(seq(1,n),w[1:5,4],type = 'l',col=colors[2], ylab = "efficiency",xlab = "n1",xlim = c(1,7),
     ylim = range(c(w[1:5,4],v[1:5,4],y[1:5,4],z[1:5,4])),main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.3)'))
lines(seq(1,n),y[1:5,4],type = 'l',col=colors[1])
lines(seq(1,n),v[1:5,4],type = 'l',col=colors[8])
lines(seq(1,n),z[1:5,4],type = 'l',col=colors[7])
lines(seq(1,n),x[1:5,4],type = 'l',col='black')
legend('topright',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',#"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],'black'),lty = 1)#,colors[c(5)],'white',4,3,6,9


plot(seq(1,n),w[1:5,7],type = 'l',col=colors[2], ylab = "efficiency",xlab = "n1",xlim = c(1,7),
     ylim = range(c(w[1:5,7],v[1:5,7],y[1:5,7],z[1:5,7])),main = paste('p=(',pfix[1],pfix[2],pfix[3], '0.5)'))
lines(seq(1,n),y[1:5,7],type = 'l',col=colors[1])
lines(seq(1,n),v[1:5,7],type = 'l',col=colors[8])
lines(seq(1,n),z[1:5,7],type = 'l',col=colors[7])
lines(seq(1,n),x[1:5,7],type = 'l',col='black')
legend('topright',c("true prior",'variance = 0.04','variance = 0.01',#'variance = 0.004',
               #'variance = 0.04','variance = 0.01','variance = 0.004',
               "wrong prior",'variance = 0.04','variance = 0.01',#"neutral prior",#'variance = 0.004',
               "frequentist rule"),col=c('white',colors[c(2,1)],'white',
                                         colors[c(8,7)],'black'),lty = 1)#,colors[c(5)],'white',4,3,6,9





#########   heatmaps
par(mfrow=c(1,1))
#hmplot(0.5,F)
hmplot()
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
