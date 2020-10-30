library(VGAM)
library(truncnorm)
library(invgamma)
library(readr)
library(matrixcalc)
library(MASS)
library(datasets)
library(matrixStats)
library(progress)


initializer <- function(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix){
  return(list(Z = Z, p = p, comp = comp, b0 = b0, B0 = B0, N=N, mu0=mu0, sigma0=sigma0, complete_matrix=complete_matrix))
  
}


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

sampler_beta <- function(mu, tau1, tau2, comp, y, B0, b0){
  n <- length(mu)
  if (!is.null(dim(comp))){
    comp = comp[,1]
  }
  
  if (is.null(dim(B0))){
    B0 <- diag(c(B0))
  }
 
  M = comp*mu
  sigma_inv = diag(ifelse(comp==1, 1/tau1^2,  1/tau2^2))
  
  B = solve(B0) + t(X) %*% sigma_inv %*% X
  B[lower.tri(B)] = t(B)[lower.tri(B)]
  b = solve(B0) %*% b0 + t(X) %*% sigma_inv %*%(Y-M)

  return(list(beta=mvrnorm(n = 1, mu = solve(B)%*%b, Sigma=solve(B)), 
              sigma_inv = sigma_inv, B0 = B0, X = X))
}

sampler_gamma_prob <- function(X, Y, beta, mu, tau1, tau2, log_rate = F){ # mu is M
  
  ww <- matrix(0, nrow=n, ncol=2)
  
  w1 <- dnorm(Y, mean = X%*%beta+mu, sd = tau1, log = log_rate)
  w2 <- dnorm(Y, mean = X%*%beta-mu, sd = tau2, log = log_rate)
  
  ww <- cbind(w1, w2)
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

  mu <- c()
  tau1 <- c()
  tau2 <- c()
  
  #Find multiplicities
  res <- table(K) #summary of labels
  pos <- as.numeric(names(res)) #positions of atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res

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
    S2_star = ifelse(is.na(sum(epsilon[K==j & comp==-1])), S2, S2+0.5*sum((epsilon[K==j & comp==-1]+Z[j,1])^2))
    
    
    Z[j,2] <- sqrt(rinvgamma(1, shape=s1_star, rate=S1_star))
    Z[j,3] <- sqrt(rinvgamma(1, shape=s2_star, rate=S2_star)) 
    
  }
  complete_matrix = Z[K,]
  return(list(Z=Z, K=K, complete_matrix = complete_matrix))
}


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
    
    predictive[,j] = sapply(grid, function(x) sum(p*(.5*dnorm(x, Z[,1], sd=Z[,2])+.5*dnorm(x, -Z[,1], sd=Z[,3]))))
    
    
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


# Diagnostic and Inference
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

coverage <- function(real_value, intervals){
  cover <- sum(real_value>=intervals[,1] & real_value<=intervals[,2])/length(beta_true)
  return(cover)
}


model_selection <- function(X, Y, sigma2, cutoff_level = 0.75){
  n = nrow(X)
  t_cutoff <- qt(cutoff_level, df=n-1)
  index <- matrix(coef(lm(Y~X[,-1])), nrow=1)%*%(t(X)%*%X)/sigma2
  selected <- index>t_cutoff
  return(X[,selected])
}

