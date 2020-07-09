library(rjags)

# EX. A
model <- "model {
y ~ dbin(theta,n) # likelihood
theta ~ dbeta(alfa,beta) # prior
ratio <- (1-theta)/theta # odds of no lactase
}"
cat(model)
dat<-data.frame(n=98,y=43,alfa=1,beta=1)
init<-list(theta=0.5,.RNG.name='base::Mersenne-Twister',.RNG.seed=1)
mod.a<-jags.model(textConnection(model),data=dat,inits=init)

post.a<-coda.samples(mod.a, variable.names=c('theta','ratio')
                     ,n.iter=1e4)
plot(post.a); summary(post.a)

qbeta(c(0.025,0.5,0.975),43+1,98-43+1)

# EX. B
dat<-data.frame(n=98,y=43,alfa=48.5,beta=51.5)
init<-list(theta=0.5,.RNG.name='base::Mersenne-Twister',.RNG.seed=1)
mod.b<-jags.model(textConnection(model),data=dat,inits=init)
post.b<-coda.samples(mod.b, variable.names=c('theta','ratio')
                     ,n.iter=1e4)
plot(post.b); summary(post.b)

# EX. C
modelc <- "model {
y ~ dbin(theta,n) # likelihood
theta1 ~ dunif(0.2,0.3)
theta2 ~ dunif(0.2,0.3)
theta <- theta1+theta2 # prior
ratio <- (1-theta)/theta # odds of no lactase
}"

init<-list(theta1=0.21,theta2=0.29,.RNG.name='base::Mersenne-Twister',.RNG.seed=1)
dat[,1:2]
mod.c<-jags.model(textConnection(modelc),data=dat,inits=init)
post.c<-coda.samples(mod.c, variable.names=c('theta','ratio')
                     ,n.iter=1e4)
plot(post.c); summary(post.c)

# EX. D
head(post.a)
mean(as.matrix(post.a)[,2]<0.485)
head(post.b)
mean(as.matrix(post.b)[,2]<0.485)
head(post.c)
mean(as.matrix(post.c)[,2]<0.485)

######
# EX 6


beta <- 2
curve(dgamma(x,2*beta,beta),0,8)
qgamma(c(0.025,0.975),2*beta,beta)  #alpha=2*beta
vis.gamma()

beta <- 1.5
qgamma(c(0.025,0.975),2*beta,beta)
alpha<-2*beta

model<-"model{
for(i in 1:20){
logsurv[i]~dnorm(mu,tau)
}
beta<- 1.5 # parameter of prior for tau
tau ~ dgamma(2*beta,beta) # prior for tau
mu ~ dnorm(log(30),tau) # prior for mu, given tau
sigma<-1/sqrt(tau)
logsurv.new~dnorm(mu,tau) # new predicted log surv time
}"
load("datasets/ibs4b.RData")

#Check model specification
cat(model)
#Check data
p6e1data
#Initial values
init<-list(mu=0, tau=1)
# To init values using a specific random number generator (for reproducibility):
# Options: "base::Wichmann-Hill" , "base:: Marsaglia-Multicarry" , "base::Super-Duper"
# If you're running several chains, then change the RNG.seed for each chain!!
init<-list(mu=0, tau=1,.RNG.name='base::Mersenne-Twister',.RNG.seed=1)
mod<-jags.model(textConnection(model),data=p6e1data[2], inits=init)
#Nodes to be monitored
nodes<-c('mu','sigma','logsurv.new'); nodes
#Draw 10,000 samples from the posterior distributions
post<-coda.samples(mod,variable.names=nodes,n.iter=1e4)
plot(post); summary(post)
library(lattice)
densityplot(post)

head(post)
mean(as.matrix(post)[,1]>log(150))
mean(as.matrix(post)[,'logsurv.new']>log(150))

#Burn-in and simulation accuracy
list.samplers(mod)
#Burn-in period ends at iteration 5000
plot(window(post,5001))
summary(window(post,5001))
# Add more iterations as we only have 5000 now
nodes<-c('mu','sigma','logsurv.new')
more<-coda.samples(mod,nodes,5000)
summary(more)
combo<- combo<-as.mcmc(rbind(as.mcmc(post),as.mcmc(more)))
# Check the number of iterations we have now and drop the first 5000
niter(combo)
summary(window(combo,5001))

#######
## PRACTICAL 7 - MCMC DIAGNOSTICS
library(rjags)

model<-"model{
y~dbin(theta,n)
theta~dunif(0,1)
ratio<-(1-theta)/theta
}"
# Create data frame
dat<-data.frame(n=98,y=43); dat
# Create initial values. Note initialisation by a list of lists
init<-list(list(theta=0.5),list(theta=1e-5),list(theta=0.99999))
init

mod7<-jags.model(textConnection(model),data=dat,inits=init, n.chains=3)
nodes<-c('theta','ratio')
post<-coda.samples(mod7,variable.names=nodes,n.iter=1e4)
summary(post) # Pooled
lapply(post,summary) #Separate

library(lattice)
# TRaceL plot of the iteration history of selected nodes
#Trace for all parameters; chains superimposed:
xyplot(post)
#Trace for all parameters, chains in separate panels:
xyplot(post,outer=T,layout=c(3,2))
#Trace for a selected parameter, chains in separate panels:
xyplot(post[,'theta'],outer=T,layout=c(1,3))
# Densities: plot of posterior distributions for the separate chains
densityplot(post,outer=T)
# chains can be superimposed:
densityplot(post)
# Quantiles: plots the evolution of sample quantiles as a function of the nb of iterations
# Default quantiles are median and 2.5 , 97.5 % quantile
# for the first chain:
cumuplot(post[[1]],col=2:1)
cumuplot(post[[1]],probs=0.5)
# for the third...
cumuplot(post[[3]],col=2:1)
# Autocorrelation plots
acfplot(post,lwd=2,auto.key=T)
acfplot(post[,'theta'],outer=T)
gelman.plot(post,lwd=2)

# EX 7 - c
burned<-window(post,8001)
head(burned)
densityplot(burned,groups=F)
summary(burned)
# Time-series SE / SD *100 < 2 -> reached burn-in
(0.0035847/0.26498)*100 < 2
(0.0006654/0.04921)*100 < 2
# examine posterior summaries separately for each chain with:
densityplot(burned)
# To examine the chains separetly:
lapply(burned,summary)
# To see if all chains are similar and consistent with convergence

# EX 7 - D
# To increase accuracy, we can run more iterations
more<-coda.samples(mod7,variable.names=nodes,n.iter=1000)
summary(more)
# append the new iterations in more to the ones in burned
library(runjags)
combo<-combine.mcmc(list(burned,more))
summary(combo)



