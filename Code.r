# install necessary packages

#install.packages('mvtnorm')
library(mvtnorm)
#install.packages('truncnorm')
library(truncnorm)

## QUESTION 1 ##

#importing and analysing data

cw <- read.csv('CW18.csv')

bio1<-as.numeric(cw$bio1)
x<- bio1
head(cw)
str(cw)
hist(cw$r)
plot(cw$d,cw$r)
plot(cw$d, cw$r,col=c("black","red")[as.factor(cw$bio1)],main = "Response to Dose received", sub = "hilighted red for the existing of biomarkers", xlab = "Dose Received", ylab="Measured Response",pch = 20)
legend(x = "bottomright", legend = levels(as.factor(cw$bio1)), col=c("black","red"), pch = 20)
boxplot(cw$d,cw$r)
summary(cw)

## QUESTION 2 ##

# checking helpful code

ll <- function(theta,r,d,x) # calculates the log likelihood 
{
  E0 <- theta[1] 
  Emax <- theta[2]
  ED50 <- theta[3]
  lambda <- theta[4]
  sigmasq <- exp(theta[5])
  beta1 <- theta[6]
  n <- length(r)
  
  if((ED50 + beta1) > 0)  # check to ensure parameters are positive
  {
    mu <- E0 + ((d^lambda)*Emax)/((d^lambda)+(ED50+beta1*x)^lambda)
    out <-  (-n/2)*(log(2*pi)+log(sigmasq))-(1/(2*sigmasq))*sum((r-mu)^2)
  } else{
    out <- -Inf
  }
  return(out)
}

pri <- function(theta) # gives the prior distributions and collective log likelihood   
{
  E0 <- theta[1]
  Emax <- theta[2]
  ED50 <- theta[3]
  lambda <- theta[4]
  sigmasq <- exp(theta[5]) # exponential since theta[5] = log(sigma) 
  beta1 <- theta[6]  
  
  E0pri <- dnorm(E0,0,10)
  Emaxpri <- dnorm(Emax,100,10)  
  pED50 <- ED50/200
  pED50pri <- dbeta(pED50,2.5,5)
  lambdapri <- dunif(lambda,0.5,3)
  sigmasqpri <- dtruncnorm(sigmasq,a=0,mean=3,sd=5)
  betapri <- dnorm(beta1,10,4)
  
  ind <- ((ED50 + beta1)>0) # indicator function to check positivity 
  
  return(log(E0pri*Emaxpri*pED50pri*lambdapri*sigmasqpri*betapri*ind))
}

## QUESTIONS 3, 4 & 5 ## include something about burn-in for Q4

# Metropolis-Hastings sampler

theta <- c(5, 115, 80, 3, 3.8, 8) # initial guess for parameters
sig.prop <- c(1.35, 2, 1.8, 0.03,0.055,1.5) # proposed vector of (tuned) standard deviations
nrep <- 100000 # number of repetitions
ll0 <- ll(theta, cw$r, cw$d, x) + pri(theta) # initial log likelihood
all.alpha <- rep(0, nrep) # empty vectors for acceptance probabilities and accepted proposals
all.accept <- rep(0, nrep)
th <- matrix(0, 6, nrep) # empty matrix for theta vectors
th[, 1] <- theta # first column of th is theta

#Metropolis-Hastings loop

for (i in 2:nrep){
  
  th[1, i] <- th[1, i - 1] + rnorm(1)*sig.prop[1] #proposed values for each parameter
  th[2, i] <- th[2, i - 1] + rnorm(1)*sig.prop[2]
  th[3, i] <- th[3, i - 1] + rnorm(1)*sig.prop[3]
  th[4, i] <- th[4, i - 1] + rnorm(1)*sig.prop[4] 
  th[5, i] <- th[5, i - 1] + rnorm(1)*sig.prop[5]
  th[6, i] <- th[6, i - 1] + rnorm(1)*sig.prop[6]

  # check positivity in 3rd and 4th parameters
  
  if (th[3, i] <= 0 || th[4, i] <= 0){ # reject proposal
    
    th[,i] <- th[,i-1] 
    
  }
    
  else {
    
    ll1 <- ll(th[,i], cw$r, cw$d, x) + pri(th[, i])
    acc <- exp(min(0, ll1 - ll0, na.rm = TRUE)) # exponential to cancel log
    all.alpha[i] <- acc }
    
    if (runif(1) <= acc){# accept proposal
      
      all.accept[i] <- 1
      ll0 <- ll1 # keeps ll0 in sync with th
      
    }
    
    else { ## reject proposal
      
      th[,i] <- th[,i-1]
      ll1 <- ll0 ## keeps ll1 in sync with th
      
    }
  
  }

prob <- sum(all.accept)/nrep # overall acceptance rate
prob

# trace plots and autocorrelation plots for each parameter

par(mfcol=c(2,3),ps = 12, cex.lab = 1.3, cex.main = 2)
plot(th[1, 1:i], type = 'l', main = expression("MH 1: Trace Plot for" ~~ {E0}),ylab = expression({E0^(i)}),xlab = expression({index: i})) 
abline(h=mean(th[1,1:i]),col=2)
acf(th[1, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E0}))
plot(th[2, 1:i], type = 'l', main = expression("MH 1: Trace Plot for" ~~ {E[max]}),ylab = expression({E[max]^(i)}),xlab = expression({index: i}))
abline(h=mean(th[2,1:i]),col=2)
acf(th[2, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E[max]}))
plot(th[3, 1:i], type = 'l',main = expression("MH 1: Trace Plot for" ~~ {ED[50]}),ylab = expression({ED[50]^(i)}),xlab = expression({index: i}))
abline(h=mean(th[3,1:i]),col=2)
acf(th[3, 1:i],lag.max = 300,main = expression("ACF for" ~~ {ED[50]}))
plot(th[4, 1:i], type = 'l',main = expression("MH 1: Trace Plot for" ~~ {lambda}),ylab = expression({lambda^(i)}),xlab = expression({index: i}))
abline(h=mean(th[4,1:i]),col=2)
acf(th[4, 1:i],lag.max = 300,main = expression("ACF for" ~~ {lambda}))
plot(th[5, 1:i], type = 'l',main = expression("MH 1: Trace Plot for" ~~ {sigma^2}),ylab = expression({sigma^2^(i)}),xlab = expression({index: i}))
abline(h=mean(th[5,1:i]),col=2)
acf(th[5, 1:i],lag.max = 300,main = expression("ACF for" ~~ {sigma^2}))
plot(th[6, 1:i], type = 'l',main = expression("MH 1: Trace Plot for" ~~ {beta}),ylab = expression({beta^(i)}),xlab = expression({index: i}))
abline(h=mean(th[6,1:i]),col=2)
acf(th[6, 1:i],lag.max = 300,main = expression("ACF for" ~~ {beta}))

## These plotting codes are for report writing
#par(mfcol=c(2,3), cex.axis = 1.8, cex.main = 2.2)

#  plot(th[1, 1:i], type = 'l', main = expression("Trace Plot for" ~~ {E0}),ylab = '',xlab = '') 
#  abline(h=mean(th[1,1:i]),col=2)
#  acf(th[1, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E0}),ylab = '',xlab = '')
#  plot(th[2, 1:i], type = 'l', main = expression("Trace Plot for" ~~ {E[max]}),ylab = '',xlab = '')
#  abline(h=mean(th[2,1:i]),col=2)
#  acf(th[2, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E[max]}),ylab = '',xlab = '')
#  plot(th[3, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {ED[50]}),ylab = '',xlab = '')
#  abline(h=mean(th[3,1:i]),col=2)
#  acf(th[3, 1:i],lag.max = 300,main = expression("ACF for" ~~ {ED[50]}),ylab = '',xlab = '')
#  plot(th[4, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {lambda}),ylab = '',xlab = '')
#  abline(h=mean(th[4,1:i]),col=2)
#  acf(th[4, 1:i],lag.max = 300,main = expression("ACF for" ~~ {lambda}),ylab = '',xlab = '')
#  plot(th[5, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {sigma^2}),ylab = '',xlab = '')
#  abline(h=mean(th[5,1:i]),col=2)
#  acf(th[5, 1:i],lag.max = 300,main = expression("ACF for" ~~ {sigma^2}),ylab = '',xlab = '')
#  plot(th[6, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {beta}),ylab = '',xlab = '')
#  abline(h=mean(th[6,1:i]),col=2)
#  acf(th[6, 1:i],lag.max = 300,main = expression("ACF for" ~~ {beta}),ylab = '',xlab = '')

## QUESTION 6 ##

#ccf(th[1, 1:i],th[2, 1:i] , lag.max = 300)

CorS=matrix(0,6,6) 
CorV=matrix(0,6,6) 
Me=matrix(0,6,1)
for (k in 1 : 6){
    Me[k]= mean(th[k,])
  for (j in 1 : 6){
    CorS[k,j] = cor( th[k,],th[j,] )
    CorV[k,j] = cov( th[k,],th[j,] )
  }
}
CorS # this is the correlation between parameters
CorV

## QUESTION 7 & 8 ##

theta <- Me
nrep1 <- 100000 # number of repetitions
ll0 <- ll(theta, cw$r, cw$d, x) + pri(theta) # initial log likelihood
all.alpha1 <- rep(0, nrep1) # empty vectors for acceptance probabilities and accepted proposals
all.accept1 <- rep(0, nrep1)
th <- matrix(0, 6, nrep1) # empty matrix for theta vectors
th[, 1] <- theta # first column of th is theta
k = 0.9


for (i in 2:nrep1){

th[1:6, i] = th[1:6, i - 1] + t(rmvnorm(1, c(0, 0, 0, 0, 0, 0), (k^2)*CorV))

if (th[3, i] <= 0 || th[4, i] <= 0)
  th[,i] <- th[,i-1] 
  

else{
  
  ll1 <- ll(th[,i], cw$r, cw$d, x) + pri(th[, i])
  acc <- exp(min(0, ll1 - ll0, na.rm = TRUE)) # exponential to cancel log
  all.alpha1[i] <- acc }

  if (runif(1) <= acc){# accept proposal
  
    all.accept1[i] <- 1
    ll0 <- ll1 # keeps ll0 in sync with th
  }
else { ## reject proposal
  
  th[,i] <- th[,i-1]
  ll1 <- ll0 ## keeps ll1 in sync with th
  
}

}
prob1 <- sum(all.accept1)/nrep1 # overall acceptance rate
prob1

par(mfcol=c(2,3),cex.axis = 1.5, cex.lab = 1.3, cex.main = 2.2)

plot(th[1, 1:i], type = 'l', main = expression("Trace Plot for" ~~ {E0}),ylab = expression({E0^(i)}),xlab = expression({index: i})) 
abline(h=mean(th[1,1:i]),col=2)
acf(th[1, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E0}))
plot(th[2, 1:i], type = 'l', main = expression("Trace Plot for" ~~ {E[max]}),ylab = expression({E[max]^(i)}),xlab = expression({index: i}))
abline(h=mean(th[2,1:i]),col=2)
acf(th[2, 1:i],lag.max = 300, main = expression("ACF for" ~~ {E[max]}))
plot(th[3, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {ED[50]}),ylab = expression({ED[50]^(i)}),xlab = expression({index: i}))
abline(h=mean(th[3,1:i]),col=2)
acf(th[3, 1:i],lag.max = 300,main = expression("ACF for" ~~ {ED[50]}))
#par(mfcol=c(1,1),cex.axis = 1.5, cex.lab = 1.3, cex.main = 2.2)
plot(th[4, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {lambda}),ylab = '',xlab = '')
abline(h=mean(th[4,1:i]),col=2)
acf(th[4, 1:i],lag.max = 300,main = expression("ACF for" ~~ {lambda}))
plot(th[5, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {sigma^2}),ylab = expression({sigma^2^(i)}),xlab = expression({index: i}))
abline(h=mean(th[5,1:i]),col=2)
acf(th[5, 1:i],lag.max = 300,main = expression("ACF for" ~~ {sigma^2}))
plot(th[6, 1:i], type = 'l',main = expression("Trace Plot for" ~~ {beta}),ylab = expression({beta^(i)}),xlab = expression({index: i}))
abline(h=mean(th[6,1:i]),col=2)
acf(th[6, 1:i],lag.max = 300,main = expression("ACF for" ~~ {beta}))

## QUESTION 9 ##

acl <- matrix(0, 6, 1)
n <- matrix(0,6,1)

for (i in 1:6){
    my_acf <- acf(th[i,], plot = FALSE, lag.max =500)
    acl[i,1] <- 2*sum(my_acf$acf) + 1
    n[i,1] <- nrep1/acl[i,1] # effective sample size
    
}

## QUESTION 10 ##

th1 <- matrix(0, 6, floor(nrep1/45))

for (i in 1:nrep1){
  
  if (i %% 45 == 0){
    
    th1[,i/45] <- th[,i]
  }
}
    
## QUESTION 11 ##

sample1 <- th1[, 1:1111]
sample2 <- th1[, 1112:2222]
ks.test(sample1, sample2) # not significant so samples are from the same distribution

## QUESTION 12 ##

hist(th1[1,]) # normal
hist(th1[2,]) # normal
hist(th1[3,]) # normal
hist(th1[4,]) # negatively skewed
hist(th1[5,]) # normal
hist(th1[6,]) # normal
par(mfcol=c(2,3))
qqnorm(th1[1,],main = expression("Q-Q plot for" ~~ {E0}))
qqnorm(th1[2,],main = expression("Q-Q plot for" ~~ {E[max]}))
qqnorm(th1[3,],main = expression("Q-Q plot for" ~~ {ED[50]}))
qqnorm(th1[4,],main = expression("Q-Q plot for" ~~ {lambda}))
qqnorm(th1[5,],main = expression("Q-Q plot for" ~~ {sigma^2}))
qqnorm(th1[6,],main = expression("Q-Q plot for" ~~ {beta}))


CorS1=matrix(0,6,6) 
Me=matrix(0,6,1)
for (k in 1 : 6){
  for (j in 1 : 6){
    CorS1[k,j] = cor( th[k,],th[j,] )
  }}

## Question 13 ##

quantile(th1[1,],c(0.025,0.975)) # wide intervals apart from theta 4 and 5
quantile(th1[2,],c(0.025,0.975))
quantile(th1[3,],c(0.025,0.975))
quantile(th1[4,],c(0.025,0.975))
quantile(th1[5,],c(0.025,0.975))
quantile(th1[6,],c(0.025,0.975))

## Question 14 ##

E0 <- mean(th1[1,])
Emax <- mean(th1[2,])
ED50 <- mean(th1[3,])
lambda <- mean(th1[4,])
sigmasq <- mean(exp(th1[5,]))
beta1 <- mean(th[6,])

r <- function(d, y){
  
  r = E0 + (d^lambda * Emax/(d^lambda + (ED50 + beta1 * y)^lambda))
  
}

par(mfcol=c(1,1))

d <- c(0:200)
plot(r(d, 0), type = 'l', col = "black",main = "Overal Dose Toxicity Model") # plots with/without biomarker
lines(r(d, 1), col = "red", type = 'l') 
points(cw$d, cw$r,col=c("black","red")[as.factor(cw$bio1)],pch = 20)
legend(x = "bottomright", legend = levels(as.factor(cw$bio1)), col=c("black","red"), pch = 20)
  
E30 <- th1[3,] * (0.3/0.7)^(1/th1[4,])
plot(E30, type = 'l')
abline(h=mean(E30),col=2)
quantile(E30,c(0.025,0.975))
mean(E30)

  
