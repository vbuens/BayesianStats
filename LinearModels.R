library(rjags)
library(runjags)

m2<-"model {
for (i in 1:11) {
y[i]~dnorm(mu[i],tau)
mu[i]<-b0+b1*x.c[i]+b2*x.c2[i]
}
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6); b2~dnorm(0,1e-6)
tau~dgamma(1e-6,1e-6); sigma<-1/sqrt(tau)
}"

#Load data
load("datasets/ibs4b.RData")

# Add xolumn x.c (centred) and x.c2 (centred and squared)
p8e1data[,'x.c']<-p8e1data[,'x']-mean(p8e1data[,'x'])
p8e1data[,'x.c2']<-p8e1data[,'x.c']**2

# For two chains, create initial values for stochastic nodes:
ini1<-list(b0=-1,b1=0, b2=-10,tau=1)
ini2<-list(b0= 1,b1=10,b2= 10,tau=20000)
ini<-list(ini1,ini2)
# List the nodes
nodes<-c('b0','b1','b2','sigma','dic')
# Fit quadratic regression model (default: burn-in 4000, draw 10k samples)
post2<-run.jags(model=m2,monitor=nodes,data=p8e1data,inits=ini, n.chains=2,burnin=0)
# Obtain summaries of posteriors and model fit assessment
post2



m1<-"model {
for (i in 1:11) {
y[i]~dnorm(mu[i],tau)
mu[i]<-b0+b1*x.c[i]
}
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6); 
tau~dgamma(1e-6,1e-6); sigma<-1/sqrt(tau)
}"
ini1<-list(b0=-1,b1=0, tau=1)
ini2<-list(b0= 1,b1=10,tau=20000)
inim1<-list(ini1,ini2)
nodes<-c('b0','b1','sigma','dic')
post1<-run.jags(m1,nodes,p8e1data,inits=inim1,n.chains=2)

m0<-"model {
for (i in 1:11) {
y[i]~dnorm(mu[i],tau)
mu[i]<-b0
}
b0~dnorm(0,1e-6)
tau~dgamma(1e-6,1e-6); sigma<-1/sqrt(tau)
}"
ini1<-list(b0=-1,tau=1)
ini2<-list(b0= 1,tau=20000)
inim0<-list(ini1,ini2)
nodes0<-c('b0','sigma','dic')
post0<-run.jags(m0,nodes0,p8e1data,inits=inim0,n.chains=2)
post0

###
## CHOICE OF NON-NESTED REGRESSION MODELS USING DIC
##

load("datasets/ibs4b.RData")
p8e2<-p8e2data
p8e2[,-1]<-scale(p8e2[,-1],scale=F)
summary(p8e2); attach(p8e2)
# Produce pairwise scatterplots
pairs(p8e2,panel=panel.smooth,lwd=2,lower.panel=NULL)
# Fit a frequentist MLR model with 3 x-variables and then reduced models for each pair of x-variables
f123<-lm(y~x1+x2+x3)
f12 <-lm(y~x1+x2)
f23 <-lm(y~x2+x3)
f13 <-lm(y~x1+x3)
AIC(f123,f12,f23,f13)

model<-"model{
for (i in 1:100) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- b0+b1*x1[i]+b2*x2[i]+b3*x3[i]
}
tau ~ dgamma(0.001,0.001); sigma <- 1/sqrt(tau)
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6)
b2~dnorm(0,1e-6); b3~dnorm(0,1e-6)
}"

# Fit Bayesian regression models
ini1<-list(b0=-200,b1=0, b2= 0,b3= 0,tau=0.001)
ini2<-list(b0= 0,b1=20,b2=-10,b3=-5,tau=1)
ini<-list(ini1,ini2)
nodes<-c('b0','b1','b2','b3','sigma','dic')
b123<-run.jags(model,nodes,p8e2,inits=ini,n.chain=2)
b123

# Fit the three models each with 2 x-variables.
# Assign zero to the regression coefficient and comment out its prior:
model<-"model{
for (i in 1:100) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- b0+b1*x1[i]+b2*x2[i]+b3*x3[i]
}
tau ~ dgamma(0.001,0.001); sigma <- 1/sqrt(tau)
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6)
b2~dnorm(0,1e-6); b3<-0
}"
# Drop the initial value of the removed coefficient:
ini12<-list(ini1[-4],ini2[-4]) # b3 is the 4th element in each list
b12<-run.jags(model,nodes,p8e2,inits=ini12,n.chain=2)
b12


model23<-"model{
for (i in 1:100) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- b0+b1*x1[i]+b2*x2[i]+b3*x3[i]
}
tau ~ dgamma(0.001,0.001); sigma <- 1/sqrt(tau)
b0~dnorm(0,1e-6); b1<-0
b2~dnorm(0,1e-6); b3~dnorm(0,1e-6)
}"
ini1<-list(b0=-200,b1=0, b2= 0,b3= 0,tau=0.001)
ini2<-list(b0= 0,b1=20,b2=-10,b3=-5,tau=1)
ini23<-list(ini1[-2],ini2[-2]) # b1 is the 2th element in each list
b23<-run.jags(model23,nodes,p8e2,inits=ini23,n.chain=2)
b23


model13<-"model{
for (i in 1:100) {
y[i] ~ dnorm(mu[i], tau)
mu[i] <- b0+b1*x1[i]+b2*x2[i]+b3*x3[i]
}
tau ~ dgamma(0.001,0.001); sigma <- 1/sqrt(tau)
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6)
b2<-0; b3~dnorm(0,1e-6)
}"
ini1<-list(b0=-200,b1=0, b2= 0,b3= 0,tau=0.001)
ini2<-list(b0= 0,b1=20,b2=-10,b3=-5,tau=1)
ini13<-list(ini1[-3],ini2[-3]) # b2 is the 3th element in each list
b13<-run.jags(model13,nodes,p8e2,inits=ini13,n.chain=2)
b13

##
## CHANGE POINT REGRESSION
##

model <- "model{
for (i in 1:41) {
y[i] ~ dnorm(mu[i],tau)
mu[i] <- b0+b1*(x[i]-cp)*step(x[i]-cp)
}
tau~dgamma(0.001,0.001); sigma<-1/sqrt(tau)
b0~dnorm(0,1e-6); b1~dnorm(0,1e-6)
cp ~ dunif(1,15)
#change.point ~ dunif(1,15)
}"
# Run 2 chains with initial values for stochastic nodes
in1<-list(b0= 0,b1= 2,tau= 0.1,change.point=2)
in2<-list(b0=10,b1=-2,tau=20, change.point=14)
in1<-list(b0= 0,b1= 2,tau= 0.1,cp=2)
in2<-list(b0=10,b1=-2,tau=20, cp=14)
ini<-list(in1,in2)
nodes<-c('b0','b1','cp','sigma','dic')

# Set burn-in at zero and dra 10k samples
post3<-run.jags(model,nodes,p8e3data,n.chains=2,inits=ini,burnin=0)

library(lattice)
xyplot(as.mcmc.list(post3)) # Trace plots
acfplot(as.mcmc.list(post3),auto.key=T) #Autocorrelation plots
densityplot(as.mcmc.list(post3)) #Posterior densities
gelman.plot(as.mcmc.list(post3))

# With 10k burn-in:
burned<-run.jags(model,nodes,data=p8e3data,n.chains=2,inits=ini,burnin=10000)
xyplot(as.mcmc.list(burned)) # Trace plots
acfplot(as.mcmc.list(burned),auto.key=T) #Autocorrelation plots
densityplot(as.mcmc.list(burned)) #Posterior densities
gelman.plot(as.mcmc.list(burned))
burned

# Simulation accuracy. Largest MC error/sd is 2.7% for CP
(0.022088/0.8059)* 100 > 2  # TRUE

# Values are lower than 5 %, so the accuracy is good enough :)

