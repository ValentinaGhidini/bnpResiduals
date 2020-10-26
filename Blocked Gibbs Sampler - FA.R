library(VGAM)
library(truncnorm)
library(invgamma)
library(readr)
library(matrixcalc)
library(MASS)
library(datasets)
library(matrixStats)
library(progress)
library(plotrix)

## Blocked Gibbs Sampler

# 1. Simulated data with particular residuals

n <- 1000
k <- 10
beta = rep(1,k+1)#c(50,100,500,20,200)# rpois(k+1,1)
beta_true = beta



#Generate residuals
muu = 30
tau1 = 4
tau2 = 1

gen_residuals <- function(n, mu, tau1, tau2){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(T,F),n, replace = T)
  res[bol] <- rnorm(sum(bol),mu,tau1)
  res[!bol] <- rnorm(n-sum(bol),-mu, tau2)
  return(res)
}

#dat <- rnorm(n)
dat <- gen_residuals(n,muu,tau1,tau2)


# Density of the residuals
plot(density(dat))


residuals <- dat 


#Generate data

X <- round(matrix(rnorm(n*k), nrow = n, ncol = k),2) #k
X <- cbind(1, X)#cbind(1,X)


#residuals <- gen_residuals(n,mu$mu,tau1,tau2)
gen_data <- function(X,beta,residuals){
  y <- X%*%beta+ residuals
  return(y)
}

Y <- gen_data(X,beta_true,residuals)
fit <- lm(Y~.,data = data.frame(X[,-1]))


### Real Data

#Get the data
setwd("C:\\Users\\valen\\Desktop\\Misture di Misture")

library(zoo)
load("data/FamaFrench.RData")
load("data/mfunds.RData")
GSATX <- mfunds$GSATX

data <- mfunds$GSATX
data <- as.vector(data[4153:length(data)])
GSATX <- rep(0,length(data)-1)
for(i in 2:length(data)){
  GSATX[i-1] <- 100*(data[i]-data[i-1])/data[i-1] #compute returns
}

dat <- window(FamaFrench, start = "2007-06-25") #start of GSATX
dat <- as.data.frame(dat)
dat <- data.frame(dat, GSATX = GSATX[1:2775]) #end of FamaFrench

Y <- dat$GSATX - dat$RF
X <- cbind(1,dat$Mkt.RF,dat$SMB,dat$HML)


b = coef(lm(Y~X))
beta_true = coef(lm(Y~X))




## Let's run the algorithm
# Initialize the values:

n <- nrow(X)
k <- ncol(X)-1

# Initialization

b0 <-  rep(0, k+1)
B0 <-  diag(10, k+1)
s1 <-  2

s2 <-  2

S1 <- 1#0.01
S2 <-  1#30

sigma2=10
theta = 1

b <- coef(lm(Y~X[,-1]))

epsilon = Y-X%*%b
comp <- ifelse(epsilon > 0, 1, -1)
N = 150

log_s = T


V <- rbeta(N-1,shape1 = 1,theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W




mu0 = 1
sigma0 = sqrt(25)
Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))



complete_matrix = matrix(1, nrow=n, ncol=3)

initializer <- list(Z = Z, p = p, comp = comp, b0 = b0, B0 = B0, N=N, mu0=mu0, sigma0=sigma0, complete_matrix=complete_matrix)



p_sampler <- function(K,theta,N){ # distribution of the Ks - IT DOES NOT DEPEND ON THE MODEL!
  #K = n classifications labels
  #theta = concentration parameter
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  M <- M[-N]
  #sample betas
  V <- rbeta(length(M),shape1 = 1+M,theta+sum(M)-cumsum(M))
  W <- cumprod(1-V)
  V <- c(V,1)
  W <- c(1,W)
  p <- V*W
  return(p)
}

### Sample the betas

sampler_beta <- function(mu, tau1, tau2, comp, y, B0, b0){
  n <- length(mu)
  #sigma = c()
  #M = c()
  #print(length(tau1))
   #print(length(tau2))
  #print(length(comp))
  #print("Mu:")
  #print(summary(mu))
  #print("Tau1:")
  #print(summary(tau1))
  #print("Tau2:")
  #print(summary(tau2))
  #print("comp:")
  #print(summary(as.factor(comp)))
  if (!is.null(dim(comp))){
    comp = comp[,1]
  }
  
  if (is.null(dim(B0))){
    B0 <- diag(c(B0))
  }
  
  #for (i in 1:n){
    
    #sigma[i]=ifelse(comp[i]==1, tau1[i]^2,  tau2[i]^2) 
    #M[i]= comp[i]*mu[i] # M is the array containing the real means (with the signs adjusted)
  #}
  M = comp*mu
  #print(M)
  #print(ifelse(comp==1, 1/tau1^2,  1/tau2^2))
  sigma_inv = diag(ifelse(comp==1, 1/tau1^2,  1/tau2^2))
  
  #sigma = diag(ifelse(comp==1, tau1^2,  tau2^2)[,1])
  
  
  #print(ifelse(comp==1, 1/tau1^2,  1/tau2^2))
  #print(sigma)
  #print(diag(sigma))
  #sigma_inv = solve(sigma)
  
  B = solve(B0) + t(X) %*% sigma_inv %*% X
  #print(B)
  B[lower.tri(B)] = t(B)[lower.tri(B)]
  b = solve(B0) %*% b0 + t(X) %*% sigma_inv %*%(Y-M)
  
  #print(solve(B))
  return(list(beta=mvrnorm(n = 1, mu = solve(B)%*%b, Sigma=solve(B)), 
              sigma_inv = sigma_inv, B0 = B0, X = X))
}


ss <- sampler_beta(complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], comp, Y, B0, b0)
B <- ss$beta
sigma <- ss$sigma



### Sample the gammas
sampler_gamma_prob <- function(X, Y, beta, mu, tau1, tau2, log_rate = F){ # mu is M
  
  ww <- matrix(0, nrow=n, ncol=2)
  
  w1 <- dnorm(Y, mean = X%*%beta+mu, sd = tau1, log = log_rate)
  w2 <- dnorm(Y, mean = X%*%beta-mu, sd = tau2, log = log_rate)
  
  ww <- cbind(w1, w2)
  #for (i in 1:n){
   # w1 <- dnorm(Y[i], mean=(X%*%beta+mu)[i,], sd=tau1[i], log = log_rate)
    #w2 <- dnorm(Y[i], mean=(X%*%beta-mu)[i,], sd=tau2[i], log = log_rate)
    
    #ww[i, ] <- c(w1, w2)
  #}
  # Normalize the weights - obtain a probability measure
  if (!log_rate){
    ww <- apply(ww, 1, function(x) {x/sum(x)})
    
  }
  
  if (log_rate){
    ww <- apply(ww, 1, function(x) {exp(x-logSumExp(x))})
  }
  ww <- t(ww)
  return(ww)
}

#w <- sampler_gamma_prob(X, Y, B, Z[,1], Z[,2], Z[,3], log_rate = F)
### Sample the triplet (mu, tau1, tau2)

sampler_triplet_blockedGS <- function(Y,epsilon, comp, p, Z, mu0, sigma0){
  # Let's sample KZ
  N <- length(p)
  n <- length(Y)
  K=c()
  weights_k1 <- log(p)+t(sapply(epsilon, function(x) dnorm(x, mean=Z[,1], sd=Z[,2], log=T)))
  weights_k2 <- log(p)+t(sapply(epsilon, function(x) dnorm(x, mean=-Z[,1], sd=Z[,3], log=T)))
  weights_k1[comp==-1,] = weights_k2[comp==-1,]
  weights_k <- apply(t(weights_k1), 1, function(x) exp(x-logSumExp(x)))#exp(weights_k1-logSumExp(weights_k1))
  K = apply(weights_k, MARGIN=1, function(x) sample(1:N, size=1, replace=T, prob=x))
  #print(K)
  #for (i in 1:n){
    #print(weights_k[i,])
   # K[i] = sample(1:N, size=1, prob = weights_k[i,])
  #}
  #for (i in 1:n){
   # weights_k <- log(p)+dnorm(epsilon[i], mean=comp[i]*Z[,1], sd=ifelse(rep(comp[i],N)==1, Z[,2], Z[,3]), log=T)
    #weights_k <- exp(weights_k-logSumExp(weights_k))
    #print(sum(weights_k))
    #print(weights_k)
    #K[i] = sample(1:N, size=1, prob = weights_k)
  #}
  #print("K:")
  #print(summary(as.factor(K)))

  mu <- c()
  tau1 <- c()
  tau2 <- c()
  
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #positions of atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
 # print(M)
 
  Z[M == 0,] = matrix(c(rtruncnorm(n=sum(M==0), a=0, b=Inf, mean = mu0, sd = sigma0), sqrt(rinvgamma(n=sum(M==0), shape=s1, rate = S1)), 
               sqrt(rinvgamma(n=sum(M==0), shape=s2, rate = S2))), ncol=3, byrow=F)

  
  index_pos <- which(M!=0)
  for (j in index_pos){
    m1 <- sum(K==j & comp==1)
    
    m2 <- sum(K==j & comp==-1)
    
    mu_tilda = (mu0*Z[j,2]^2*Z[j,3]^2
                +m1*mean(epsilon[comp==1])*sigma0^2*Z[j,3]^2
                +m2*mean(-epsilon[comp==-1])*sigma0^2*Z[j,2]^2)/(Z[j,2]^2*Z[j,3]^2
                                                                  +m1*sigma0^2*Z[j,3]^2+m2*sigma0^2*Z[j,2]^2)

    sigma2_tilda = (1/sigma0^2+m1/Z[j,2]^2+m2/Z[j,3]^2)^(-1)
    
    Z[j,1] <- rtruncnorm(n = 1, mean=mu_tilda, sd=sqrt(sigma2_tilda), a = 0)
    
    s1_star = s1+m1/2
    S1_star = ifelse(is.na(sum(epsilon[K==j & comp==1])), S1, S1+.5*sum((epsilon[K==j & comp==1]-Z[j,1])^2))
    
    s2_star = s2+m2/2
    #print(m2)
    S2_star = ifelse(is.na(sum(epsilon[K==j & comp==-1])), S2, S2+0.5*sum((epsilon[K==j & comp==-1]+Z[j,1])^2))
    
    
    Z[j,2] <- sqrt(rinvgamma(1, shape=s1_star, rate=S1_star)) #controlla scale-location
    Z[j,3] <- sqrt(rinvgamma(1, shape=s2_star, rate=S2_star)) #controlla scale-location
    
  }
  #print(Z)
  complete_matrix = Z[K,]
  #print(s1_star)
  #print(S1_star)
  #print(s2_star)
  #print(S2_star)
  return(list(Z=Z, K=K, complete_matrix = complete_matrix))
}



gs <- sampler_triplet_blockedGS(Y,epsilon, comp,p, Z, mu0, sigma0)



sampler <- function(X, Y, initializer, iter, burn_in = 10000, thin = 50, Verbose=T, log_rate=F){
  
  null_df <- data.frame(NULL)
  
  write.table(null_df, file = "mu.csv", row.names = F)
  write.table(null_df, file = "tau1.csv", row.names = F)
  write.table(null_df, file = "tau2.csv", row.names = F)
  write.table(null_df, file = "weights.csv", row.names = F)
  write.table(null_df, file = "betas.csv", row.names = F)
  write.table(null_df, file = "component.csv", row.names = F)
  write.table(null_df, file = "epsilon.csv", row.names = F)
  
  # grid to compute the predictive
  grid <- seq(-50, 50, by=.1)
  predictive <- matrix(NA, nrow=length(grid), ncol=iter)
  
  
  if (thin<=0){
    thin=1
  }
  
  # Hyperpriors
  t=2
  T_big=4
  
  #sigma0 = sqrt(rinvgamma(n=1, shape = t, rate = T_big))
  #mu0 = rtruncnorm(1, a=0, b=Inf, mean =0, sd = sigma0)
  mu0=initializer$mu0
  sigma0=initializer$sigma0
  
  pb <- progress_bar$new(total = iter)
  Z <- initializer$Z
  b0 <- initializer$b0
  B0 <- initializer$B0
  #w <- initializer$w
  comp = initializer$comp
  p <- initializer$p
  #mu0 = initializer$mu0
  #sigma0 = initializer$sigma0
  complete_matrix = initializer$complete_matrix
  times_beta <- c()
  times_gamma <- c()
  times_gs <- c()
  
  for (j in 1:iter){
    # print(j)
    
    if (Verbose){
      pb$tick()
      Sys.sleep(1 /iter)}
    
    
    begin_beta <- Sys.time()
    
    
    # Betasss
    #if (Verbose){print("beta")}
    ss <- sampler_beta(complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], comp, Y, B0, b0)
    
    end_beta <- Sys.time()
    times_beta[j] <- end_beta - begin_beta
    beta = ss$beta
    epsilon <- Y-X%*%beta
    
    # Gamma prob
    begin_gamma <- Sys.time()
    w <- sampler_gamma_prob(X, Y, beta, complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], log_rate = log_rate)
    
    end_gamma <- Sys.time()
    times_gamma[j] <- end_gamma - begin_gamma
    #Update comp
    comp = apply(w,1,function(x){
      return(sample(c(1,-1), size=1, prob = x))
    })
    
    
    #if (Verbose){print("Triplet - Z")}
    # mu sampler
    begin_gs <- Sys.time()
    gs <- sampler_triplet_blockedGS(Y, epsilon, comp, p, Z, mu0, sigma0)
    end_gs <- Sys.time()
    
    times_gs[j] = end_gs - begin_gs
    
    K <- gs$K
    Z <- gs$Z
    
    complete_matrix <- gs$complete_matrix
    
    p <- p_sampler(K, theta, N)
    
    # Predictive
    #print(sapply(grid, function(x) sum(p*(.5*dnorm(x, Z[,1], sd=Z[,2])+.5*dnorm(x, -Z[,1], sd=Z[,3])))))
    predictive[,j] = sapply(grid, function(x) sum(p*(.5*dnorm(x, Z[,1], sd=Z[,2])+.5*dnorm(x, -Z[,1], sd=Z[,3]))))
    #print(predictive)
    
    
    ll <- list(beta = beta, component = comp, 
               weights = w, mu = complete_matrix[,1], tau1 = complete_matrix[,2], tau2=complete_matrix[,3], epsilon = epsilon)
    
    
    
    if (j >= burn_in & j%%thin ==0){
      # Append to file
      write.table(matrix(ll$beta, nrow=1), "betas.csv", sep = ";", col.names = !file.exists("betas.csv"), append = T, row.names = F)
      write.table(matrix(ll$component, nrow=1), "component.csv", sep = ";", col.names = !file.exists("component.csv"), append = T, row.names = F)
      write.table(matrix(ll$weights, nrow=1), "weights.csv", sep = ";", col.names = !file.exists("weights.csv"), append = T, row.names = F)
      write.table(matrix(ll$mu, nrow=1), "mu.csv", sep = ";", col.names = !file.exists("mu.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau1, nrow=1), "tau1.csv", sep = ";", col.names = !file.exists("tau1.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau2, nrow=1), "tau2.csv", sep = ";", col.names = !file.exists("tau2.csv"), append = T, row.names = F)
      write.table(matrix(ll$epsilon, nrow=1), "epsilon.csv", sep = ";", col.names = !file.exists("epsilon.csv"), append = T, row.names = F)
      
      
    }
    
    
  }
  return(predictive)
}




credible_intervals <- function(parameters, level = 0.05){
  if (is.null(dim(parameters))){
    parametrs = as.matrix(parameters)
    
  }
  intervals <- matrix(NA, ncol=2, nrow=ncol(parameters))
  for (i in 1:ncol(parameters)){
    intervals[i,] = c(quantile(parameters[,i], level/2), quantile(parameters[,i], 1-level/2))
  }
  return(intervals)
}

#mu0 = 10
#sigma0 = sqrt(5)


begin = Sys.time()
print(begin)
trial <- sampler(X, Y, initializer, iter=5000, burn_in = 100, thin = 10, log_rate = log_s )
end = Sys.time()
end - begin 

predictive_estimates = apply(trial, 1, mean)


eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")

m <- read.table("betas.csv", header=F, skip=1, sep=";")


plot(density(dat), type="l", lwd=2, main = "Density of the residuals \n First and Last Iteration")

plot(density(as.numeric(eps[1,])), type="l", col=2, lty=3, lwd=2)

lines(density(as.numeric(eps[dim(eps)[1],])), type="l", col=3, lty=4, lwd=2)



tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

plot(mu[,12], type="l")

plot(tau1[,90], type="l")

plot(tau2[,90], type="l")

to_study = 3


beta_est = m[,to_study]



acf(beta_est)

plot(beta_est, type="l")
abline(h=beta_true[to_study], col=2) 


abline(h=mean(beta_est, na.rm=T), col=3, lty=2)


d <- data.frame(cbind(beta_true[to_study], mean(beta_est)))

colnames(d) <- c("True Value", "Estimate")
rownames(d) <- paste("Beta", to_study)
d


betas_estimates <- apply(m, 2, mean)

d <- cbind(betas_estimates, beta_true)
colnames(d) = c("estimates", "true values")
rownames(d) = NULL

d

cred_int = credible_intervals(m)
### OLS linear regression


ols_model <- lm(Y~X[,-1])

plot(ols_model)

coef(ols_model)

plot(beta_true, type="l", main="OLS coefficients vs Proposed coefficients", ylim=c(0,2))
# add fill
polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(cred_int[ ,2]), cred_int[ ,1]), col = rgb(1,0,0, alpha=0.2), border = NA)
# model
# intervals
lines(1:(k+1), cred_int[ ,2], lty = 'dashed', col = 'red')
lines(1:(k+1), cred_int[ ,1], lty = 'dashed', col = 'red')

lines(betas_estimates, col=2)
points(betas_estimates, col=2, pch=19)


polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(confint(ols_model)[ ,2]), confint(ols_model)[ ,1]), col = rgb(0,1,0, alpha=0.2), border = NA)
# model
# intervals
lines(1:(k+1), confint(ols_model)[ ,2], lty = 'dashed', col = 'green')
lines(1:(k+1), confint(ols_model)[ ,1], lty = 'dashed', col = 'green')

#plotCI(x=1:(k+1),y=betas_estimates,li=t(cred_int[,1]),ui=t(cred_int[,2]),col=2, add=T)

lines(coef(ols_model), col=3)
points(coef(ols_model), col=3, pch=19)
legend("bottomright", y = NULL, legend=c("OLS", "Proposed"), fill = NULL, col = c(3,2),
       pch=19)


#plotCI(x=1:(k+1),y=coef(ols_model),li=t(confint(ols_model)[,1]),ui=t(confint(ols_model)[,2]),col=3, add=T)


coverage <- function(real_value, intervals){
  cover <- sum(real_value>=intervals[,1] & real_value<=intervals[,2])/length(beta_true)
  return(cover)
}

coverage(beta_true, cred_int)
coverage(coef(ols_model), confint(ols_model))
