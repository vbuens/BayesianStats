# Hypothesis: 
# presence of dwarf mutation has a significant effect on leaf unrolling compared to WT
library(rjags)
library(runjags)

#Load data
load("datasets/ibs4b.RData")


model<-"model {
  for(i in 1:25){
    unrolling[i] ~ dnorm(mu[i],tau)
    mu[i] <- mu.wt*x.wt[i] + mu.dm*x.dm[i]
  }
  mu.wt ~ dnorm(0,1e-6); mu.dm ~ dnorm(0,1e-6)
  tau ~ dgamma(0.001,0.001); sigma <- 1/sqrt(tau)
}"

muta<-function(x)  return(list(
    diff.wt.dm = x[,'mu.wt'] - x[,'mu.dm'],
    ratio.wt.dm = x[,'mu.wt'] / x[,'mu.dm']))

# Initialisation
in1<-list(mu.wt=0,mu.dm=4,tau=1)
in2<-list(mu.wt=3,mu.dm=1,tau=20)
ini<-list(in1,in2)
nodes<-c('mu.wt','mu.dm','sigma')

# Extract posteriors
post<-run.jags(model,nodes,p9e1data,n.chains=2,ini,mutate=muta)

# Diagnostic plots
xyplot(as.mcmc.list(post)) # Trace plots
acfplot(as.mcmc.list(post),auto.key=T) #Autocorrelation plots
densityplot(as.mcmc.list(post)) #Posterior densities
# Print posterior summaries and densities
post; densityplot(as.mcmc(post))

# Chance that ratio.wt.dm is larger than 1.2 ?
mean(as.mcmc(post)[,'ratio.wt.dm']>=1.2)


##
## BAYESIAN HYPOTHESIS TESTING
##

# Option A. Use 95% HPD credible interval
post$HPD 
# Result: CI: 0.02146, 0.57068). Mean: 0.299

# Option BCompute a quantity analogous to the frequentist p-value, called contour probability (pb)
# Exploiting symmetric shape, determining Pr(diff.wt.dm <= 0) and then multiply by 2
mean(as.mcmc(post)[,'diff.wt.dm']<=0)*2 
# Result: 0.0345

## EXPRESING GROUP MEANS IN A LINEAR MODEL FRAMEWORK

# Get mut variety only (line B), to have 2 leaf sections and 3 treatments
p9e2sub<-subset(p9e2data,line=='B' & mutation=='Mut')
# Create boxplot with treatments:
boxplot(unrolling~Treatment, main='Leaf unrolling by treatment',data=p9e2sub)


# Evaluate the following combinations:
# ratio of means of t1-control, t2-control, t2-t1
# difference of average effect of t2 and t1 combined vs control

muta<-function(x)  return(list(
    ratio.t1.c =x[,'mu.t1']/x[,'mu.c'],  #t1/control
    ratio.t2.c =x[,'mu.t2']/x[,'mu.c'],  #t2/control
    ratio.t2.t1=x[,'mu.t2']/x[,'mu.t1'],
    diff.t1t2.c=(x[,'mu.t1']+x[,'mu.t2'])/2-x[,'mu.c']))
muta


model3<- "model {
for(i in 1:30){
unrolling[i] ~ dnorm(mu[i],tau)
mu[i] <- mu.c*x.c[i] + mu.t1*x.t1[i] + mu.t2*x.t2[i]
}
tau ~dgamma(0.001,0.001); sigma<- 1/sqrt(tau)
mu.c ~dnorm(0,1e-6)
mu.t1~dnorm(0,1e-6)
mu.t2~dnorm(0,1e-6)
}"

in1<-list(mu.c=0,mu.t1=4,mu.t2=-1,tau=1)
in2<-list(mu.c=3,mu.t1=0,mu.t2= 5,tau=20)
ini<-list(in1,in2)
#List the nodes to monitor.
nodes<-c('mu.c','mu.t1','mu.t2','sigma')
post2<-run.jags(model3,nodes,p9e2sub,inits=ini,n.chains=2,mutate=muta, burnin=0)

# Diagnostic plots
library(lattice)
post2m<-as.mcmc.list(post2)
xyplot(post2m,layout=c(2,4),as.table=T,main='trace plots')
acfplot(post2m,layout=c(4,2),type='l',as.table=T, main='autocorrelation plots')
densityplot(post2m,layout=c(4,2),as.table=T,aspect=1,main='density plots')
gelman.plot(post2m,lwd=2,auto.layout=F) #auto.layout=T

# From the plots, we decide a burn-in of 6k is needed, but to be conservative, we use 10k:

burned<-run.jags(model3,nodes,p9e2sub,inits=ini,n.chains=2, mutate=muta, burnin=10000)
burnedm<-as.mcmc.list(burned)
# Diagnostic plots
xyplot(burnedm,layout=c(2,4),as.table=T,main='trace plots after a burn-in of 10,000')
acfplot(burnedm,layout=c(4,2),type='l',as.table=T, main='autocorrelation plots after a burn-in of 10,000')
densityplot(burnedm,layout=c(4,2),as.table=T,aspect=1,main='density plots')
gelman.plot(burnedm,lwd=2,auto.layout=F) #auto.layout=T

# Posterior summaries
burnedmc<-as.mcmc(burned)
densityplot(burnedmc,layout=c(4,2),as.table=T,aspect=1, main='density plots after a burn-in of 10,000')

burned
# Simmulation accuracy: (MCerr/SD) 
(burned$summaries[,7]/burned$summaries[,5])*100
burned$summaries[,8]
# All values are less than 0.8 %. Based on this statistic, we don't need more accuracy

# Probability that mean of t2 is at least 1/5 larger than t1:
mean(burnedmc[,'ratio.t2.t1']>=1.2)
#  0.10945 - i.e. 1 in 9

# Assessing Null Hypothesis: Average effect of t1 + t2 combined is no different to control
mean(burnedmc[,'diff.t1t2.c']<=0)*2
#Result: 0.0036 - i.e. reject H0. There is a 95% chance that the average combined effect is
# between 0.17 and 0.73 mm larger than that of control. (see diff.t1t2.c)


