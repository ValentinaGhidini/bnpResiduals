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


p_sampler <- function(K,theta,N){
  # distribution of the Ks -
  # K = n classifications labels
  # theta = concentration parameter
  # Find multiplicities
  res <- table(K) # summary of labels
  pos <- as.numeric(names(res)) # atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  M <- M[-N]
  # sample betas
  V <- rbeta(length(M),shape1 = 1+M,theta+sum(M)-cumsum(M))
  W <- cumprod(1-V)
  V <- c(V,1)
  W <- c(1,W)
  p <- V*W
  return(p)
}

sampler_beta <- function(mu, tau1, tau2, comp, Y, X, B0, b0){
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
  for (i in 1:n){
    weights_k <- log(p)+dnorm(epsilon[i], mean=comp[i]*Z[,1], sd=ifelse(rep(comp[i],N)==1, Z[,2], Z[,3]), log=T)
    weights_k <- exp(weights_k-logSumExp(weights_k))
    K[i] = sample(1:N, size=1, prob = weights_k)
  }
  #K = apply(weights_k, MARGIN=1, function(x) sample(1:N, size=1, replace=T, prob=x))

  mu <- c()
  tau1 <- c()
  tau2 <- c()

  #Find multiplicities
  res <- table(K) # summary of labels
  pos <- as.numeric(names(res)) # positions of atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res

  Z[M == 0,] = matrix(c(rtruncnorm(n=sum(M==0), a=0, b=Inf, mean = mu0, sd = sigma0), sqrt(rinvgamma(n=sum(M==0), shape=s1, rate = S1)),
                        sqrt(rinvgamma(n=sum(M==0), shape=s2, rate = S2))), ncol=3, byrow=F)


  index_pos <- which(M!=0)
  for (j in index_pos){
    m1 <- sum(K==j & comp==1)

    m2 <- sum(K==j & comp==-1)


    mu1 <- ifelse(m1 == 0,0,mean(epsilon[K == j & comp==1]))
    mu2 <- ifelse(m2 == 0,0,mean(-epsilon[K == j & comp==-1]))

    mu_tilda = (mu0*Z[j,2]^2*Z[j,3]^2
                +m1*mu1*sigma0^2*Z[j,3]^2
                +m2*mu2*sigma0^2*Z[j,2]^2)/(Z[j,2]^2*Z[j,3]^2+m1*sigma0^2*Z[j,3]^2+m2*sigma0^2*Z[j,2]^2)


    sigma2_tilda = (1/sigma0^2+m1/Z[j,2]^2+m2/Z[j,3]^2)^(-1)

    Z[j,1] <- rtruncnorm(n = 1, mean=mu_tilda, sd=sqrt(sigma2_tilda), a = 0, b=Inf)

    s1_star = s1+m1/2
    S1_star = ifelse(m1 == 0, S1, S1+.5*sum((epsilon[K==j & comp==1]-Z[j,1])^2))




    Z[j,2] <- sqrt(rinvgamma(1, shape=s1_star, rate=S1_star))
    Z[j,3] <- Z[j,2]

  }
  complete_matrix = Z[K,]
  #print("Zeta")
  #print(Z)
  return(list(Z=Z, K=K, complete_matrix = complete_matrix))
}

bnp.lm_symmetric <- function(X, Y, initializer, iter, burn_in = 10000, thin = 50, Verbose=T, log_rate=T){
  n <- nrow(X)
  k <- ncol(X)
  # Write on file
  null_df <- data.frame(NULL)
  Grid <- seq(-50, 50, by=.1)
  predictive <- matrix(0, ncol=length(Grid), nrow=iter)
  betas <- matrix(0, ncol=iter, nrow=k)
  #write.table(null_df, file = "mu.csv", row.names = F)
  #write.table(null_df, file = "tau1.csv", row.names = F)
  #write.table(null_df, file = "tau2.csv", row.names = F)
  #write.table(null_df, file = "weights.csv", row.names = F)
  #write.table(null_df, file = "betas.csv", row.names = F)
  #write.table(null_df, file = "component.csv", row.names = F)
  #write.table(null_df, file = "epsilon.csv", row.names = F)
  if (thin<=0){
    thin=1
  }

  # Hyperpriors
  t=2
  T_big=4

  mu0=0 # To ensure semiconjugacy
  sigma0=initializer$sigma0

  pb <- progress_bar$new(total = iter)
  Z <- initializer$Z
  b0 <- initializer$b0
  B0 <- initializer$B0
  comp <- initializer$comp
  p <- initializer$p

  sigma0_saved <- c()
  complete_matrix = initializer$complete_matrix
  # Tests
  p0 <- c()
  p1 <- c()
  p0_single <- c()
  p1_single <- c()
  tests <- matrix(0, nrow = iter, ncol = k)

  # Cluster analysis for outliers
  clusters <- matrix(0, nrow = iter, ncol=n)


  for (j in 1:iter){
    if (Verbose){
      pb$tick()
      Sys.sleep(1 /iter)}

    # Sample beta

    ss <- sampler_beta(complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], comp, Y, X, B0, b0)
    beta <- ss$beta

    betas[,j] <- beta
    # Update Residuals
    epsilon <- Y-X%*%beta
    #print("residuals")
    #print(epsilon[index_to_study])

    # Gamma prob
    w <- sampler_gamma_prob(X, Y, beta, complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], log_rate = log_rate)

    #Update Gamma (comp)
    comp = apply(w,1,function(x){
      return(sample(c(1,-1), size=1, prob = x))
    })

    # Sample mu

    gs <- sampler_triplet_blockedGS(Y, epsilon, comp, p, Z, mu0, sigma0)

    # Labels - they identify the cluster of belonging for each observation
    K <- gs$K
    clusters[j,] <- K
    Z <- gs$Z
    complete_matrix <- gs$complete_matrix
    #print(summary(K))
    #print(summary(complete_matrix[,1]))
    # Update p - stick breaking probabilities
    p <- p_sampler(K, theta, N)
    #print(K[index_to_study])
    #print(p)
    #print(complete_matrix[index_to_study,])


    K_unique <- unique(K)

    # Update sigma0 - thanks to semiconjugacy
    sigma0 = sqrt(rinvgamma(n=1, t + length(K_unique)/2, T_big+0.5*sum(Z[K_unique,1]^2)))


    # Compute empirical predictive distribution on a grid of points
    predictive[j,] = sapply(Grid, function(x) sum(p*(.5*dnorm(x, Z[,1], sd=Z[,2])+.5*dnorm(x, -Z[,1], sd=Z[,3]))))
    ll <- list(beta = beta, component = comp,
               weights = w, mu = complete_matrix[,1], tau1 = complete_matrix[,2], tau2=complete_matrix[,3], epsilon = epsilon)

    # Hypothesis testig: H0: beta=0
    p0[j] <- sum(dnorm(Y, comp*Z[K,1], sd=ifelse(comp==1, Z[K,2], Z[K,3]), log = T))
    p1[j] <- sum(dnorm(Y, comp*Z[K,1]+X%*%beta, sd=ifelse(comp==1, Z[K,2], Z[K,3]), log = T))

    # Hypothesis testing for each beta:
    if (j == 1){
      for (kk in 1:k){
        beta_h0 <- beta
        beta_h0[kk] <- 0
        p0_single[kk] <- sum(dnorm(Y, comp*Z[K,1]+X%*%beta_h0 , sd=ifelse(comp==1, Z[K,2], Z[K,3]), log = T))
        p1_single[kk] <- logSumExp(p1)
      }
    }

    if (j != 1){
      for (kk in 1:k){
        beta_h0 <- beta
        beta_h0[kk] <- 0
        p0_single[kk] <- logSumExp(c(p0_single[kk], sum(dnorm(Y, comp*Z[K,1]+X%*%beta_h0 , sd=ifelse(comp==1, Z[K,2], Z[K,3]), log = T))))
        p1_single[kk] <- logSumExp(p1)
      }
    }


    if (j >= burn_in & j%%thin ==0){
      # Append to file
      #write.table(matrix(ll$beta, nrow=1), "betas.csv", sep = ";", col.names = !file.exists("betas.csv"), append = T, row.names = F)
      #write.table(matrix(ll$component, nrow=1), "component.csv", sep = ";", col.names = !file.exists("component.csv"), append = T, row.names = F)
      #write.table(matrix(ll$weights, nrow=1), "weights.csv", sep = ";", col.names = !file.exists("weights.csv"), append = T, row.names = F)
      #write.table(matrix(ll$mu, nrow=1), "mu.csv", sep = ";", col.names = !file.exists("mu.csv"), append = T, row.names = F)
      #write.table(matrix(ll$tau1, nrow=1), "tau1.csv", sep = ";", col.names = !file.exists("tau1.csv"), append = T, row.names = F)
      #write.table(matrix(ll$tau2, nrow=1), "tau2.csv", sep = ";", col.names = !file.exists("tau2.csv"), append = T, row.names = F)
      #write.table(matrix(ll$epsilon, nrow=1), "epsilon.csv", sep = ";", col.names = !file.exists("epsilon.csv"), append = T, row.names = F)

      sigma0_saved = c(sigma0_saved, sigma0)
    }

  }
  return(list(coef = betas))
}


## Final model

bnp.lm_symmetric_resultis <- function(X, Y, initializer, cluster = F, iter=10000, burn_in = 5000, thin = 10, conf_level=0.05){
  samples <- sampler(X, Y, initializer_, iter=iter, burn_in = burn_in, thin = thin, log_rate = T )
  #eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
  #m <- read.table("betas.csv", header=F, skip=1, sep=";")
  #tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
  #tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
  #mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

  # Coefficients
  betas_estimates <- apply(m, 2, mean)

  # Empirical Credible Intervals
  #cred_int = credible_intervals(m, level = conf_level)

  #d <- cbind(betas_estimates, cred_int)
  #colnames(d) <- c("Estimate", paste0(conf_level*100/2, "%"), paste0(100-conf_level*100/2, "%"))
  #Plot Residual Density
  #plot(density(as.numeric(eps[dim(eps)[1],])), type="l", col=1, lty=4, lwd=2, main = "Density of the Empirical Residuals")

  # Clusters
  if (cluster){#
    residuals <- as.numeric(eps[dim(eps)[1],]) #apply(eps, 2, mean)
    freq_clusters <- apply(samples$clusters, 2,table) # list
    label <- c()

    for (i in 1:length(freq_clusters)){
      label[i] = min(as.numeric(names(which(freq_clusters[[i]]==max(freq_clusters[[i]])))))
    }
    #x11()
    #par(mfrow=c(1,4))
    #plot(residuals, col=label , pch=19, ylab="Residual",
         #xlab="", main="Clustered residuals", cex=min(label, 3))#ifelse(label==label[which.min(residuals)] | label==label[which.max(residuals)], 1.5,1))

    #n_clusters <- length(unique(label)) # number of clusters
    #means_per_obs <- apply(mu, 2, mean)
    #tau1_per_obs <- apply(tau1, 2, mean)
    #tau2_per_obs <- apply(tau2, 2, mean)


    #plot(means_per_obs, pch=19, col=label,#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 2, 1),
    #     xlab="Observation", ylab=expression(mu), main = "Mean", cex=min(label, 3))#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 1.5,.6))


    #plot(tau1_per_obs, pch=19, col=label,#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 2, 1),
    #     xlab="Observation", ylab="", main = "First Variance", cex=min(label, 3))#felse(label==label[label[which.min(residuals)]]| label==label[label[which.max(residuals)]], 1.5,.6))


    #plot(tau2_per_obs, pch=19, col=label,#ifelse(label==label[which.min(residuals)] | label==label[label[which.max(residuals)]], 2, 1),
       #  xlab="Observation", ylab="", main = "Second Variance",  cex=min(label,3))#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 1.5,.6))

    return(list(coef = samples$coef, point_estimates = betas_estimates))


  }

  return(list(coef = d))


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

