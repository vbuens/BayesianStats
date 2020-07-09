# First download JAGS-4.3.0 from: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
install.packages("rjags")
install.packages("runjags")
install.packages("modeest")
install.packages("car")

library(rjags)
library(runjags)
library(modeest)
library(car)


### PRACTICAL 3
curve(log(5+x),0,40)
curve(dbeta(x,2,4),0,1)
qbeta(c(0.025,0.975),9,17)

# Ex. 1 A
a<-8  # 7+1
b<-14  # 20-7 +1
curve(dbeta(x,a,b),0,1,xlab = expression(theta), ylab = 'Prob. Density', main='Posterior distribution')
qbeta(c(0.025,0.975),a,b)
#Solutions: 0.1810716 0.5696755
qtl<-c(0.025,0.975)
qbeta(qtl,8,14)

# Ex. 1 B
#v = a +b = 5 +2 = 7
#a = v* µ =7 * 0.2 = 1.4
#b = v(1-µ)= 7(1-0.2) = 5.6
# PRIOR: Beta(1.4,5.6)   POSTERIOR = Beta(7+1.4,20-7 + 5.6)
a<-8.4
b<-13+5.6 
curve(dbeta(x,a,b),0,1,xlab = expression(theta), ylab = 'Prob. Density', main='Prior and Posterior')
curve(dbeta(x,1.4,5.6),0,1, add=T, lty=2)
legend(0.6,3,c("prior","posterior"),lty=c(2,1),bty="n")
qbeta(qtl,a,b)
#[1] 0.1546589 0.4940403

# Ex. 2
# mode = 10, theta < 50
# IN gamma distribution: mode(expectation)= (a-1)/beta -> a = 10B +1
Beta <- c(0.01,0.1,1,10,100)
Alpha <- 10*Beta+1  
qgamma(0.95,Alpha,Beta)
# 318.66641  47.43865  16.96222  11.80793  10.53603
Beta<-0.094
Alpha <- 10*Beta+1  
qgamma(0.95,Alpha,Beta)
curve(dgamma(x,Alpha,Beta),0,80,ylab="Prob.density")
abline(v=10,lty=2)

y<-c(24,25,31,31,22,21,26,20,16,22); y
sum(y) #238  n=10 
# Gamma dsitribution = 1.94 +238 , 0.09 +10) =(239.94,10.094)
Alpha_p <- 239.94
Beta_p  <- 10.094
curve(dgamma(x,Alpha_p,Beta_p),0,80,ylab="Prob.density")
curve(dgamma(x,Alpha,Beta), 0,80, add=T, lty=2)
legend(40,0.2,c("prior","posterior"),lty=c(2,1),bty="n")
posteriormean <- Alpha_p/Beta_p #23.77056
qgamma(qtl,Alpha_p,Beta_p) # 20.85775 26.87097
qtl <- c(0.025, 0.975)
qgamma(qtl,239.9,10.9)

