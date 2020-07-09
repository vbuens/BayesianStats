library(rjags)
library(runjags)

#Load data
load("datasets/ibs4b.RData")

#Inspect data
head(p10e1data)

meanlogdose<-mean(p10e1data$logdose.c) # result= -2.218278e-17
model<-"model {
for(i in 1:5) {
y[i] ~ dbin(pi[i],n[i])
logit(pi[i]) <- b0 + b1*logdose.c[i]
}
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6)
LD50<-exp(-b0/b1+-2.218278e-17) # FUNCTION OF NODES
}"

# Initialisation
in1<-list(b0=1,b1=4)
in2<-list(b0=3,b1=1)
ini<-list(in1,in2)
nodes<-c('b0','b1','LD50')
post<-run.jags(model,nodes,p10e1data,inits=ini,n.chains=2,burnin=0)

library(lattice)
postm<-as.mcmc.list(post)
xyplot(postm,main='trace plots',as.table=T)
acfplot(postm,type='l')
densityplot(postm,as.table=T,main='density plots',layout=c(3,1), aspect=1)
par(mfrow=c(2,2))
gelman.plot(postm,auto.layout=F,lwd=2,ylim=c(1,1.2))
mtext('Gelman-Rubin-Brooks plots',line=-2,outer=T,cex=1.2)

# Looking at the Gelman-Rubin-Brook plots, a burn-in of 1000 might be sufficient
# However, we can take a burn-in of 5000 just to be safe
burned<-run.jags(model,nodes,p10e1data,inits=ini,n.chains=2, burnin=5000)
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
# Accuracy: all are less than 1% (0.9). Good accuracy :)

###
## BINOMIAL REGRESSION IN A BIOASSAY STUDY
###
# Four bettles were sprayed with one of 2 insecticides (A or B), then reaction was measured

#Inspect the data
head(p10e2data)

# First, calculate logarithmic transform of deposit and centre it (logdepo.c), add to data
p10e2data[,'logdepo'] <-log(p10e2data$deposit)
p10e2data[,'logdepo.c'] <-log(p10e2data$deposit)-mean(log(p10e2data$deposit))
# Then, add binary indication variable for each insecticide (A or B)
p10e2data[,'A'] <-(p10e2data$insecticide=='A')*1
p10e2data[,'B'] <-(p10e2data$insecticide=='B')*1
# Add proportion of beetles reacting (pob) and its logic (logitpob):
p10e2data$pob <-p10e2data$reacted/p10e2data$trials
p10e2data[,'logitpob'] <-qlogis(p10e2data$pob)

#Check columns have been added correctly
head(p10e2data)

# Plot logitpob against log(deposit) using the insecticide letter as symbol
plot(p10e2data$logitpob~p10e2data$logdepo.c, pch=as.numeric(p10e2data$insecticide)+64)
# Looking at the graph, A appears to be more potent


# A "parallel lines" logistic regression model seems appropriate for modelling the data
model<-"model {
for(i in 1:12) {
reacted[i] ~ dbin(pi[i],trials[i])
logit(pi[i]) <- b0A*A[i] + b0B*B[i] + b1*logdepo.c[i]
}
b0A ~ dnorm(0,1e-6); b0B ~ dnorm(0,1e-6); b1 ~ dnorm(0,1e-6)
}"

# To prevent warnings from run.jags() about "unused variable", we will get a subset of the data
# keeping reacted, trials, logdepo.c, A, B
p10e2dat<-p10e2data[,c(3:4,5,7:8)]; colnames(p10e2dat)

# Now we are going to calculate:
# 1. Difference between insecticide A and B (b0A - b0B)
# 2. Dose that makes 1 out of 2 beetles react (for each insecticide, i.e. ED50)
# 3. Relative potency of the 2 insecticides - RP_AB = ED50_B/ED50_A


#Define functions of stochastic nodes to pass to run.jags()
muta<-function(x) return(list(
  diffAB=x[,'b0A']-x[,'b0B'],
  ED50A=exp(-x[,'b0A']/x[,'b1']+1.385998),  # mean(log(deposit) = 1.3859..
  ED50B=exp(-x[,'b0B']/x[,'b1']+1.385998),
  rp.AB=exp((x[,'b0A']-x[,'b0B'])/x[,'b1'])))


#Init and list the nodes to monitor
in1<-list(b0A=-2,b0B= 1,b1=-1)
in2<-list(b0A= 2,b0B=-3,b1= 5)
ini<-list(in1,in2)
nodes<-c('b0A','b0B','b1')
# SET BURNIN=0 TO DECIDE HOW LONG IT MUST BE FOR CONVERGENCE TO BE REACHED
# By default 10,000 MCMC samples are drawn:
post2<-run.jags(model,nodes,p10e2dat,burnin=0,inits=ini,n.chains=2, mutate=muta)
#Plots before selecting burn-in using lattice and coda packages:
post2m<-as.mcmc.list(post2)
library(lattice)
xyplot(post2m,as.table=T,main='trace plots')
acfplot(post2m,main='autocorrelation plots',type='l',aspect=1,as.table=T, layout=c(4,2))
densityplot(post2m,main='density plots',aspect=1,as.table=T,layout=c(4,2))
# par(mfrow=c(2,4))
gelman.plot(post2m,auto.layout=F,lwd=2)
mtext('Gelman-Rubin-Brooks plots',line=-2,col=4,outer=T)

# Burn-in period of 4000 is the minimum, let's get 10k
burned2<-run.jags(model,nodes,p10e2dat,burnin=10000,inits=ini,n.chains=2, mutate=muta)
burned2m<-as.mcmc.list(burned2)
# Plotting
xyplot(burned2m,as.table=T,main='trace plots after a burn-in of 10,000')
acfplot(burned2m,main='autocorrelation plots after a burn-in of 10,000',type='l',aspect=1,as.table=T,layout=c(4,2))
burned2mc<-as.mcmc(burned2)
densityplot(burned2mc,main='density plots after a burn-in of 10,000',aspect=1,as.table=T,layout=c(4,2))
burned2

# Accuracy: none of them are greater than 1%
burned2$summaries[,'MC%ofSD'] > 1 

# Get 95% HPD credible interval for all 4 parameters
round(burned2$summaries[4:7,c(4,1,3)],2)

# What is the probability of A being 1.4 times larger than B?
mean(burned2mc[,'rp.AB']>=1.4)
# Result: 0.28
