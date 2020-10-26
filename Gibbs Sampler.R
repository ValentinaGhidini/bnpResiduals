library(VGAM)
library(truncnorm)
library(invgamma)
library(readr)
library(matrixcalc)
library(MASS)
library(datasets)
library(matrixStats)
library(progress)




setwd("C:\\Users\\valen\\Desktop\\Misture di Misture")

set.seed(2020)

#### Sample the coefficients

## Sampler for Beta coefficients
sampler_beta <- function(mu, tau1, tau2,comp, y, B0, b0){
  sigma = c()
  M = c()
  
  if (is.null(dim(B0))){
    B0 <- diag(c(B0))
  }
  
  for (i in 1:n){
    
    sigma[i]=ifelse(comp[i]==1, tau1[i]^2,  tau2[i]^2) 
    M[i]= comp[i]*mu[i] # M is the array containing the real means (with the signs adjusted)
  }
  sigma = diag(sigma)
  sigma_inv = solve(sigma)

  B = solve(B0) + t(X) %*% solve(sigma) %*% X
  B[lower.tri(B)] = t(B)[lower.tri(B)]
  b = solve(B0) %*% b0 + t(X) %*% sigma_inv %*%(Y-M)
  
  
  return(list(beta=mvrnorm(n = 1, mu = solve(B)%*%b, Sigma=solve(B)), 
              sigma_inv = sigma_inv, sigma=sigma, M = M, B0 = B0, X = X))
}


## Sampler for gamma
sampler_gamma_prob <- function(X, Y, beta, mu, tau1, tau2, log_rate = F){ # mu is M
  ww <- matrix(0, nrow=n, ncol=2)
  
  for (i in 1:n){
    w1 <- dnorm(Y[i], mean=(X%*%beta+mu)[i,], sd=tau1[i], log = log_rate)
    w2 <- dnorm(Y[i], mean=(X%*%beta-mu)[i,], sd=tau2[i], log = log_rate)
    
    ww[i, ] <- c(w1, w2)
  }
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


## Sampler for mu
mu_sampler <- function(M,tau1, tau2, epsilon, sigma, comp, log_rate = F){
  
  mu = abs(M)
  
  if (!log_rate){
    for (i in 1:n){
      
      weights <- dnorm(epsilon[i], mean=comp[i]*mu, sd=diag(sigma)[i])
      weights[i] <- theta*2*dnorm(epsilon[i], mean = 0, sd = sqrt(sigma2+diag(sigma)[i]^2))*
        pnorm(0, mean=comp[i]*sigma2/(sigma2+diag(sigma)[i]^2)*epsilon[i], sd= sqrt(sigma2*diag(sigma)[i]^2/(sigma2+diag(sigma)[i]^2)), 
              lower.tail = F)
      weights[tau1 != tau1[i] | tau2 != tau2[i]] <- 0
      weights <- weights/sum(weights)
      
      
      mixture_component <- sample(1:n, prob=weights, size=1) 

      
      if (mixture_component==i){
        mu[i] <- rtruncnorm(1, a=0, b=Inf, mean = (comp[i]*epsilon[i]*(sigma2/(sigma2+diag(sigma)[i]^2))), sd = sqrt(sigma2*diag(sigma)[i]^2/(sigma2+diag(sigma)[i]^2)))
        
      }
      
      else{
        mu[i] <- mu[mixture_component]
        
      }
      
    }
    
  }
  if(log_rate){
    for (i in 1:n){
      
      weights <- dnorm(epsilon[i], mean=comp[i]*mu, sd=diag(sigma)[i], log = T)
      weights[i] <- log(theta)+log(2)+dnorm(epsilon[i], mean = 0, sd = sqrt(sigma2+diag(sigma)[i]^2),log=T)+
        pnorm(0, mean=comp[i]*sigma2/(sigma2+diag(sigma)[i]^2)*epsilon[i], sd= sqrt(sigma2*diag(sigma)[i]^2/(sigma2+diag(sigma)[i]^2)), 
              lower.tail = F, log=T)
      
      weights[tau1 != tau1[i] | tau2 != tau2[i]] <- -Inf
      weights <- exp(weights-logSumExp(weights))
      mixture_component <- sample(1:n, prob=weights, size=1) 
      
      
      if (mixture_component==i){
        mu[i] <- rtruncnorm(1, a=0, b=Inf, mean = (comp[i]*epsilon[i]*(sigma2/(sigma2+diag(sigma)[i]^2))), sd = sqrt(sigma2*diag(sigma)[i]^2/(sigma2+diag(sigma)[i]^2)))
        
      }
      
      else{
        mu[i] <- mu[mixture_component]
        
      }
      
      
      
    }
  }
  #print(mu)
  return(list(mu=mu))
  
}



### Sampler for tau
tau_sampler <- function(tau1, tau2, M, epsilon, s1, S1, s2, S2, comp, log_rate = F){
  
  tau_1 <- tau1
  tau_2 <- tau2
  mu = abs(M)
  for (i in 1:n){
    s1_star <- s1 + .5
    S1_star <- S1 + 0.5*(epsilon[i]-M[i])^2
    
    s2_star <- s2 + .5
    S2_star <- S2 + 0.5*(epsilon[i]-M[i])^2
  if (!log_rate){
    
    weights11 <- dnorm(epsilon[i], mean=M[i], sd=tau1)
    weights11[i] <- theta/sqrt(2*pi)*(S1)^s1/gamma(s1)*gamma(s1_star)/(S1_star^s1_star)
    weights11[mu != mu[i] | tau_2 != tau_2[i]] = 0
    weights11 <- weights11/sum(weights11)
    
    weights12 <- rep(1,n)
    weights12[i] <- theta
    weights12[mu != mu[i] | tau_2 != tau_2[i]] = 0
    weights12 <- weights12/sum(weights12)
    
    weights21 <- rep(1,n)
    weights21[i] <- theta
    weights21[mu != mu[i] | tau_1 != tau_1[i]] = 0
    weights21 <- weights21/sum(weights21)
    
    
    weights22 <- dnorm(epsilon[i], mean=M[i], sd=tau2)
    weights22[i] <- theta/sqrt(2*pi)*(S2)^s2/gamma(s2)*gamma(s2_star)/(S2_star^s2_star)
    weights22[mu != mu[i] | tau_1 != tau_1[i]] = 0
    weights22 <- weights22/sum(weights22)
    }
  
    if (log_rate){
      weights11 <- dnorm(epsilon[i], mean=M[i], sd=tau1, log=T)
      weights11[i] <- log(theta)-log(sqrt(2*pi))+s1*log(S1)-log(gamma(s1))+log(gamma(s1_star))-s1_star*log(S1_star)
      weights11[mu != mu[i] | tau_2 != tau_2[i]] = -Inf
      weights11 <- exp(weights11-logSumExp(weights11))
      
      weights12 <- rep(0,n)
      weights12[i] <- log(theta)
      weights12[mu != mu[i] | tau_2 != tau_2[i]] = -Inf
      weights12 <- exp(weights12-logSumExp(weights12))
      
      weights21 <- rep(0,n)
      weights21[i] <- log(theta)
      weights21[mu != mu[i] | tau_1 != tau_1[i]] = -Inf
      weights21 <- exp(weights21-logSumExp(weights21))
      
      
      weights22 <- dnorm(epsilon[i], mean=M[i], sd=tau2, log=T)
      weights22[i] <- log(theta)-log(sqrt(2*pi))+s2*log(S2)-log(gamma(s2))+log(gamma(s2_star))-s2_star*log(S2_star)
      weights22[mu != mu[i] | tau_1 != tau_1[i]] = -Inf
      weights22 <- exp(weights22-logSumExp(weights22))
      
    }
    
    mixture_component11 <- sample(1:n,prob=weights11,size=1)
    mixture_component22 <- sample(1:n,prob=weights22,size=1)
    mixture_component12 <- sample(1:n,prob=weights12,size=1)
    mixture_component21 <- sample(1:n,prob=weights21,size=1)
    
    if (comp[i]==1){ 
      if (mixture_component11==i){
        # sample from a new value
        tau_1[i] <- sqrt(rinvgamma(1, shape=s1_star, rate = S1_star))
      }
      else{
        # sample from the already observed values
        tau_1[i] <- tau_1[mixture_component11]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
        
      }
      
      if (mixture_component12==i){
        # sample from a new value
        tau_2[i] <- sqrt(rinvgamma(1, shape=s2, rate = S2))
      }
      else{
        # sample from the already observed values
        tau_2[i] <- tau_2[mixture_component12]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
        
      }
      
      }
    
    
    if (comp[i]==-1){ 
      
      if (mixture_component22==i){
        # sample from a new value
        tau_2[i] <- sqrt(rinvgamma(1, shape=s2_star, rate = S2_star))
      }
      else{
        # sample from the already observed values
        tau_2[i] <- tau_2[mixture_component22]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
        
      }
      
      if (mixture_component21==i){
        # sample from a new value
        tau_1[i] <- sqrt(rinvgamma(1, shape=s1, rate = S1))
      }
      else{
        # sample from the already observed values
        tau_1[i] <- tau_1[mixture_component21]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
        
      }
       }
    
  }
  return(list(tau1=tau_1, tau2=tau_2))
}



## Overall Sampler
sampler <- function(X, Y, initializer, iter, burn_in = 10000, thin = 50, Verbose=T, log_rate=F){
  mu <- initializer$mu
  tau1 <- initializer$tau1
  tau2 <- initializer$tau2
  b0 <- initializer$b0
  B0 <- initializer$B0
  #w <- initializer$w
  comp = initializer$comp
  s1 <- initializer$s1
  S1 <- initializer$S1
  s2 <- initializer$s2
  S2 <- initializer$S2
  pb <- progress_bar$new(total = iter)
  prec_sigma <- NULL
  
  for (j in 1:iter){
    
    if (Verbose){
      pb$tick()
      Sys.sleep(1 /iter)}
    
    # Betasss
    #if (Verbose){print("beta")}
    skip_to_next <- FALSE
    ss <- sampler_beta(mu$mu, tau1, tau2, comp, Y, B0, b0)

    
    if(skip_to_next) { next } 
    
    B <- ss$beta
    sigma <- ss$sigma
    M <- ss$M

    
    # Sampling Residuals 
    #if (Verbose) {print("Epsilon")}
    epsilon <- Y-X%*%B
    
    #if (Verbose){print("gamma")}
    
    # Gamma prob
    w <- sampler_gamma_prob(X, Y, B, mu$mu, tau1, tau2, log_rate = log_rate)
    
    #Update comp
    comp = apply(w,1,function(x){
      return(sample(c(1,-1), size=1, prob = x))
      })
    
    
    #if (Verbose){print("mu")}
    # mu sampler
    mu <- mu_sampler(M,tau1, tau2, epsilon, sqrt(sigma), comp, log_rate = log_rate)
    
    # tau sampler
    #if (Verbose) {print("tau")}
    
    tau <- tau_sampler(tau1, tau2, M,epsilon, s1, S1, s2, S2, comp, log_rate = log_rate)
    tau1<- tau$tau1 #sd!
    tau2<- tau$tau2

    
    ll <- list(beta = B, component = comp, 
               weights = w, mu = mu$mu, tau1 = tau$tau1, tau2=tau$tau2, epsilon = epsilon)
     
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
}



##################### DATA - choose the dataset to apply the algorithm ########################?


# 1. Simulated data with particular residuals

n <- 1000
k <- 5
beta = rep(1,k+1)#c(50,100,500,20,200)# rpois(k+1,1)
beta_true = beta


#Generate residuals
muu = 30
tau1 = 1
tau2 = 1
gen_residuals <- function(n, mu, tau1, tau2){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(T,F),n, replace = T)
  res[bol] <- rnorm(sum(bol),mu,sqrt(tau1))
  res[!bol] <- rnorm(n-sum(bol),-mu, sqrt(tau2))
  return(res)
}

dat <- gen_residuals(100,muu,tau1,tau2)

# Density of the residuals
plot(density(dat))

residuals <- dat 
#Generate data

X <- round(matrix(rnorm(n*k), nrow = n, ncol = k),2)
X <- cbind(1,X)


#residuals <- gen_residuals(n,mu$mu,tau1,tau2)
gen_data <- function(X,beta,residuals){
  y <- X%*%beta+ residuals
  return(y)
}

Y <- gen_data(X,beta_true,residuals)
fit <- lm(Y~.,data = data.frame(X[,-1]))


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

sigma2=25
theta = 1


b <- coef(fit)

epsilon = Y-X%*%b
comp <- ifelse(epsilon > 0, 1, -1)


mu <- rep(0,n)
mu[comp == 1] <- mean(epsilon[epsilon > 0])
mu[comp == -1] <- abs(mean(epsilon[epsilon < 0]))


tau1 <- 1# sqrt(rinvgamma(n,shape = s1, rate = S1)) #standard deviations
tau2 <- 1 #sqrt(rinvgamma(n,shape = s2, rate = S2)) #standard deviations
tau <- rep(0,n)
s1_star <- s1 + .5
S1_star <- S1 + 0.5*(epsilon[1]-comp[1]*mu[1])^2
s2_star <- s2 + .5
S2_star <- S2 + 0.5*(epsilon[1]-comp[1]*mu[1])^2

if(comp[1] == 1){
  tau1[1] <- sqrt(rinvgamma(1, shape=s1_star, rate = S1_star))
}else{
  tau2[1] <- sqrt(rinvgamma(1, shape=s2_star, rate = S2_star))
}
for(i in 2:n){
  s1_star <- s1 + .5
  S1_star <- S1 + 0.5*(epsilon[i]-comp[i]*mu[i])^2
  s2_star <- s2 + .5
  S2_star <- S2 + 0.5*(epsilon[i]-comp[i]*mu[i])^2
  
  #comp == 1
  weights11 <- c(dnorm(epsilon[i], mean=comp[i]*mu[i], sd=tau1[1:(i-1)], log=T),0)
  weights11[i] <- log(theta)-log(sqrt(2*pi))+s1*log(S1)-log(gamma(s1))+log(gamma(s1_star))-s1_star*log(S1_star)
  weights11 <- exp(weights11-logSumExp(weights11))
  
  weights12 <- rep(0,i)
  weights12[i] <- log(theta)
  weights12 <- exp(weights12-logSumExp(weights12))
  
  #comp == -1
  weights21 <- rep(0,i)
  weights21[i] <- log(theta)
  weights21 <- exp(weights21-logSumExp(weights21))
  
  
  weights22 <- c(dnorm(epsilon[i], mean=-mu[i], sd=tau2[1:(i-1)], log=T),0)
  weights22[i] <- log(theta)-log(sqrt(2*pi))+s2*log(S2)-log(gamma(s2))+log(gamma(s2_star))-s2_star*log(S2_star)
  weights22 <- exp(weights22-logSumExp(weights22))
  
  mixture_component11 <- sample(1:i,prob=weights11,size=1)
  mixture_component22 <- sample(1:i,prob=weights22,size=1)
  mixture_component12 <- sample(1:i,prob=weights12,size=1)
  mixture_component21 <- sample(1:i,prob=weights21,size=1)
  
  if (comp[i]==1){ 
    if (mixture_component11==i){
      # sample from a new value
      tau1[i] <- sqrt(rinvgamma(1, shape=s1_star, rate = S1_star))
    }
    else{
      # sample from the already observed values
      tau1[i] <- tau1[mixture_component11]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
      
    }
    
    if (mixture_component12==i){
      # sample from a new value
      tau2[i] <- sqrt(rinvgamma(1, shape=s2, rate = S2))
    }
    else{
      # sample from the already observed values
      tau2[i] <- tau2[mixture_component12]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
      
    }
    
  }
  
  
  if (comp[i]==-1){ 
    
    if (mixture_component22==i){
      # sample from a new value
      tau2[i] <- sqrt(rinvgamma(1, shape=s2_star, rate = S2_star))
    }
    else{
      # sample from the already observed values
      tau2[i] <- tau2[mixture_component22]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
      
    }
    
    if (mixture_component21==i){
      # sample from a new value
      tau1[i] <- sqrt(rinvgamma(1, shape=s1, rate = S1))
    }
    else{
      # sample from the already observed values
      tau1[i] <- tau1[mixture_component21]#sample(x=tau1, size=1, prob= weights1[,2]) mu[i] <- mu[mixture_component]
      
    }
  }
  
}



mu <- list(mu = mu)
log_s = T


initializer <- list(mu = mu, comp = comp, tau1=tau1, tau2=tau2, b0 = b0, B0 = B0, s1=s1, S1=S1, s2=s2, S2=S2)


# Let's create empty dataframes to save the values of the sampled coefficients

null_df <- data.frame(NULL)

write.table(null_df, file = "mu.csv", row.names = F)
write.table(null_df, file = "tau1.csv", row.names = F)
write.table(null_df, file = "tau2.csv", row.names = F)
write.table(null_df, file = "weights.csv", row.names = F)
write.table(null_df, file = "betas.csv", row.names = F)
write.table(null_df, file = "component.csv", row.names = F)
write.table(null_df, file = "epsilon.csv", row.names = F)

# Let's finally run the algorithm
log_s = T
iter = 11000
begin = Sys.time()
res <- sampler(X, Y, initializer, iter, log_rate = log_s, Verbose = T, burn_in = 1000, thin = 2)

end = Sys.time()
end-begin
### Now, let's check the results

# Let's read the datasets created
eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
m <- read.table("betas.csv", header=F, skip=1, sep=";")


plot(density(dat), type="l", lwd=2)
lines(density(as.numeric(eps[10,])), type="l", col=2)

lines(density(as.numeric(eps[50,])), type="l", col=3)

lines(density(as.numeric(eps[100,])), type="l", col=4)
lines(density(as.numeric(eps[400,])), type="l", col=5)


tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")
plot(mu[,1], type="l")


# Let's study the estimate of a coefficient
# Choose here what coefficient to study
to_study = 2

beta_est = m[,to_study]

acf(beta_est)

plot(beta_est, type="l")
abline(h=beta_true[to_study], col=2) 
abline(h=mean(beta_est, na.rm=T), col=3, lty=2)


d <- data.frame(cbind(beta_true[to_study], mean(beta_est)))
colnames(d) <- c("True Value", "Estimate")
rownames(d) <- paste("Beta", to_study)
d

## Let's see the estimate for all the coefficients

# Beta stimate
#m <- data.frame(read_table("betas.csv", col_names = T))

est_b= apply(m, mean, MARGIN=2)
est_b

acf(beta_est)
fit <- lm(Y~X)

Y_tilda = X%*%betas_est
Y_tilda
mse_new = sqrt(sum((Y-Y_tilda)^2))
mse_old = sqrt(sum((Y-predict(fit, data.frame(X[,-1])))^2))

# Let's study the estimate of the errors


### Burn in & Thinning

plot(beta_est[5000:length(beta_est)], type="l")
# Maybe the burn in is not needed?
beta_burned <- beta_est[10000:length(beta_est)]
plot(beta_burned, type="l")
# Thinning

beta_thin <- c()
i = 1
thin = 1000
j = 1
while (i <= length(beta_burned)){
  beta_thin[j] <- beta_burned[i]
  i = i+thin
  j = j+1
  
}





plot(beta_est, type="l", col="grey")
plot(beta_thin, type="l")
abline(h=beta_true[to_study], col=2) 
abline(h=mean(beta_thin), col=3, lty=2)

acf(beta_thin
  )

d <- data.frame(cbind(beta_true[to_study], mean(beta_thin)))

colnames(d) <- c("True Value", "Estimate with Thinning")
rownames(d) <- paste("Beta", to_study)

d


# Studying the residuals

betas_estimates <- apply(m, MARGIN=2,mean)


# Each residual has one mean and two variances associated. Let's study them

# Mean



means_mu <- apply(mu, MARGIN=2, mean)


acf(mu_df[,1])

for (i in 1:nrow(mu)){
  
  if (i%%100==0){
    print(mean(as.numeric(mu[i,])))
  }
  
}
#acf(mu_df[,1])

plot(density(as.numeric(mu[1000,])), type="l")
# Variance1

tau1 <- read.table("C:/Users/valen/Desktop/Misture di Misture/tau1.csv", quote="\n", comment.char="", skip=1)


tau1_df <- data.frame(NULL)
el = 1

for (j in 1:n){
  i = 1
  while (i <= iter){
    
    tau1_df[j,i] = tau1[el,1] # ogni iterazione ha una colonna
    i = i+1
    el = el+1
    
  }
  
}

means_tau1 <- apply(tau1_df, MARGIN=2, mean)



# Variance2

tau2 <- read.table("C:/Users/valen/Desktop/Misture di Misture/tau2.csv", quote="\n", comment.char="", skip=1)

head(tau2)

tau2_df <- data.frame(NULL)
el = 1

for (j in 1:n){
  i = 1
  while (i <= iter){
    
    tau2_df[j,i] = tau2[el,1] # ogni iterazione ha una colonna
    i = i+1
    el = el+1
    
  }
  
}



means_tau2 <- apply(tau2_df, MARGIN=2, mean) #1


# Resulting epsilon
eps_sim <- c()

for (i in 1:n){
  k = sample(1,c(1,2))
  if (k ==1){
    eps_sim[i] <- rnorm(1, means_mu[i], sd=means_tau1[i])
  }
  else{
    eps_sim[i] <- rnorm(1,-means_mu[i], sd=means_tau2[i])
  }
}

par(mfrow=c(1,2)) 
plot(density(residuals), main="Residuals") 
plot(density(eps_sim), main="Simulated Residuals")


