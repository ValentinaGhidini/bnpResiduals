#' @import progress
#' @import VGAM
#' @import truncnorm
#' @import invgamma
#' @import readr
#' @import matrixcalc
#' @import MASS
#' @import matrixStats
NULL

#' Initialize the blocked Gibbs sampler algorithm
#'
#' This function allows you to define the initialization values for the blocked Gibbs sampler algorithm, with the quantities illustrated in Ascolani, Ghidini (2021).
#' All the quantities in input can be random, since they should not impact the final linear model estimation.
#' @param Z matrix containing the unique values of the triplets (Z1, Z2, Z3), each one associated with a label in 1:N
#' @param p distribution of the N labels
#' @param comp vector of dimension n (= number of observations), with possible entries +1 or -1, denoting the element of the mixture from which the residual has been sampled from
#' @param b0 mean vector for the gaussian prior of beta - to specify just in case of linear models
#' @param B0 covariance matrix for the gaussian prior of beta - to specify just in case of linear models
#' @param N number of elements in the truncation of the Pitman-Yor process
#' @param mu0 hyperparameter for the baseline measure (mean of the Truncated Normal)
#' @param sigma0 hyperparameter for the baseline measure (variance of the Truncated Normal)
#' @param complete_matrix matrix of dimension nx3 containing all the triplets associated to each observation
#' @keywords initializer
#' @examples
#' initializer <- initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)
#' @export

initializer <- function(Z, p, comp, N, mu0, sigma0, complete_matrix, b0=NULL, B0=NULL ){
  return(list(Z = Z, p = p, comp = comp, b0 = b0, B0 = B0, N=N, mu0=mu0, sigma0=sigma0, complete_matrix=complete_matrix))

}

#' Sample the distribution of the labels
#'
#' This function allows you to sample labels k=1, ..., N from the distribution defined in Ishwaran and James (2001)
#' @param K set of classifications labels.
#' @param theta concentration parameter.
#' @param N number of elements in the truncation of the PY process.
#' @param sigma parameter to regulate the Pitman-Yor process - default: sigma=0 (Dirichlet process case).
#' @keywords p_sampler
p_sampler <- function(K,theta,N, sigma=0){
  res <- table(K) # summary of labels
  pos <- as.numeric(names(res)) # atoms with at least one label
  M <- rep(0,N)
  M[pos] <- res
  M <- M[-N]
  # sample betas
  V <- rbeta(length(M), shape1 = 1+M-sigma,  shape2 = theta+sum(M)-cumsum(M) + sigma*c(1:length(M))) ## To check!
  W <- cumprod(1-V)
  V <- c(V,1)
  W <- c(1,W)
  p <- V*W
  return(p)
}

#' Sample the coefficients of the linear model
#'
#' This function allows you to sample from the posterior distribution of the parameters beta, as explained in Ascolani, Ghidini (2021)
#' @param mu absolute value of the mean of the components of the residual mixture
#' @param tau1 variance of the first component of the residual mixture
#' @param tau2 variance of the second component of the residual mixture
#' @param comp vector of dimension n, with possible entries +1 or -1, denoting the element of the mixture from which the residual has been sampled from
#' @param Y response variable
#' @param X design matrix
#' @param b0 mean vector for the gaussian prior of beta
#' @param B0 covariance matrix for the gaussian prior of beta
#' @keywords sampler_beta, beta_sampler
sampler_beta <- function(mu, tau1, tau2, comp, Y, X, B0, b0){
  n <- length(Y)
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
              sigma_inv = sigma_inv))
}

#' Compute the distribution of gamma
#'
#' This function allows you to compute the probability distribution of gamma, as explained in Ascolani, Ghidini (2021)
#' @param X design matrix
#' @param Y response variable
#' @param beta coefficients of the linear model
#' @param mu absolute value of the mean of the components of the residual mixture
#' @param tau1 variance of the first component of the residual mixture
#' @param tau2 variance of the second component of the residual mixture
#' @param log_rate compute all the probabilities in log-scale
#' @keywords sampler_gamma_prob
sampler_gamma_prob <- function(epsilon, mu, tau1, tau2, log_rate = T){ # mu is M

  ww <- matrix(0, nrow=n, ncol=2)

  w1 <- dnorm(epsilon, mean = mu, sd = tau1, log = log_rate)
  w2 <- dnorm(epsilon, mean = -mu, sd = tau2, log = log_rate)

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


#' Blocked Gibbs sampler algorithm
#'
#' This function allows you to sample an instance of (Z1, Z2, Z3) from the posterior distribution,
#' exploiting the Blocked Gibbs Sampler algorithm as explained in Ascolani, Ghidini (3032)
#' @param Y response variable
#' @param epsilon residuals
#' @param comp vector of dimension n, with possible entries +1 or -1, denoting the element of the mixture from which the residual has been sampled from
#' @param Z matrix containing the unique values of the triplets (Z1, Z2, Z3), each one associated with a label in 1:N
#' @param p distribution of the N labels
#' @param mu0 hyperparameter for the baseline measure (mean of the Truncated Normal)
#' @param sigma0 hyperparameter for the baseline measure (variance of the Truncated Normal)
#' @keywords sampler_triplet_blockedGS, blocked gibbs sampler

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

    s2_star = s2+m2/2
    S2_star = ifelse(m2 == 0, S2, S2+0.5*sum((epsilon[K==j & comp==-1]+Z[j,1])^2))


    Z[j,2] <- sqrt(rinvgamma(1, shape=s1_star, rate=S1_star))
    Z[j,3] <- sqrt(rinvgamma(1, shape=s2_star, rate=S2_star))

  }
  complete_matrix = Z[K,]
  return(list(Z=Z, K=K, complete_matrix = complete_matrix))
}




#' Sample all the parameters defining the model
#'
#' This function allows you to obtain a Markov Chain with stationary distribution equal to the posterior of the parameters sampled, as explained in Ascolani, Ghidini (year)
#' @param X design matrix
#' @param Y response variable
#' @param initializer list of initial values - can be obtained as the output of the function initializer()
#' @param iter number of iterations to consider
#' @param burn_in number of initial sampled valued to discard
#' @param thin thinning value
#' @param Verbose boolean value which allows to print information about the running of the algorithm
#' @param log_rate compute all the probabilities in log-scale (default: log_rate=T)
#' @param sigma parameter of the Pitman Yor process - default: sigma=0 (Dirichlet process case)
#' @keywords sampler
sampler <- function(X, Y, initializer, iter, burn_in = 10000, thin = 50, Verbose=T, log_rate=T, sigma=0){
  n <- nrow(X)
  k <- ncol(X)
  # Write on file
  null_df <- data.frame(NULL)
  Grid <- seq(-50, 50, by=.1)
  predictive <- matrix(0, ncol=length(Grid), nrow=iter)
  write.table(null_df, file = "mu.csv", row.names = F)
  write.table(null_df, file = "tau1.csv", row.names = F)
  write.table(null_df, file = "tau2.csv", row.names = F)
  write.table(null_df, file = "weights.csv", row.names = F)
  write.table(null_df, file = "betas.csv", row.names = F)
  write.table(null_df, file = "component.csv", row.names = F)
  write.table(null_df, file = "epsilon.csv", row.names = F)
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

    # Update Residuals
    epsilon <- Y-X%*%beta

    # Gamma prob
    w <- sampler_gamma_prob(epsilon, complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], log_rate = log_rate)

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

    # Update p - stick breaking probabilities
    p <- p_sampler(K, theta, N, sigma)
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
      write.table(matrix(ll$beta, nrow=1), "betas.csv", sep = ";", col.names = !file.exists("betas.csv"), append = T, row.names = F)
      write.table(matrix(ll$component, nrow=1), "component.csv", sep = ";", col.names = !file.exists("component.csv"), append = T, row.names = F)
      write.table(matrix(ll$weights, nrow=1), "weights.csv", sep = ";", col.names = !file.exists("weights.csv"), append = T, row.names = F)
      write.table(matrix(ll$mu, nrow=1), "mu.csv", sep = ";", col.names = !file.exists("mu.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau1, nrow=1), "tau1.csv", sep = ";", col.names = !file.exists("tau1.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau2, nrow=1), "tau2.csv", sep = ";", col.names = !file.exists("tau2.csv"), append = T, row.names = F)
      write.table(matrix(ll$epsilon, nrow=1), "epsilon.csv", sep = ";", col.names = !file.exists("epsilon.csv"), append = T, row.names = F)

      sigma0_saved = c(sigma0_saved, sigma0)
    }

    }
  return(list(predictive=predictive, sigma0=sigma0_saved, clusters=clusters, global_test = cbind(p0, p1),
              single_test_avg = cbind(p0_single, p1_single)))
}

#' Fitting a Bayesian Nonparametric linear model
#'
#' This function allows you to ofit the linear model defined in Ascolani, Ghidini (year)
#' @param X design matrix
#' @param Y response variable
#' @param initializer list of initial values - can be obtained as the output of the function initializer()
#' @param iter number of iterations to consider
#' @param burn_in number of initial sampled valued to discard
#' @param thin thinning value
#' @param conf_level confidence level of the credible intervals
#' @param cluster boolean value indicating whether the cluster analysis is desired
#' @param sigma parameter to regulate the Pitman-Yor process - default: sigma=0 (Dirichlet process case).
#' @keywords model, linear model, lm
#' @export

bnp.lm <- function(X, Y, initializer, cluster = F, iter=10000, burn_in = 5000, thin = 10, conf_level=0.05, sigma=0){
  samples <- sampler(X, Y, initializer_, iter=iter, burn_in = burn_in, thin = thin, log_rate = T, , sigma=sigma)
  eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
  m <- read.table("betas.csv", header=F, skip=1, sep=";")
  tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
  tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
  mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

  # Coefficients
  betas_estimates <- apply(m, 2, mean)

  # Empirical Credible Intervals
  cred_int = credible_intervals(m, level = conf_level)

  d <- cbind(betas_estimates, cred_int)
  colnames(d) <- c("Estimate", paste0(conf_level*100/2, "%"), paste0(100-conf_level*100/2, "%"))

  #Plot Residual Density
  x11()
  plot(density(as.numeric(eps[dim(eps)[1],])), type="l", col=1, lty=4, lwd=2, main = "Density of the Empirical Residuals")


  # Clusters
  if (cluster){
    residuals <- apply(eps, 2, mean)
    freq_clusters <- apply(samples$clusters, 2,table) # list
    label <- c()

    for (i in 1:length(freq_clusters)){
      label[i] = min(as.numeric(names(which(freq_clusters[[i]]==max(freq_clusters[[i]])))))
    }
    x11()
    par(mfrow=c(1,4))
    plot(residuals, col=label , pch=19, ylab="Residual",
         xlab="", main="Clustered residuals", cex=min(label, 3))#ifelse(label==label[which.min(residuals)] | label==label[which.max(residuals)], 1.5,1))

    n_clusters <- length(unique(label)) # number of clusters
    means_per_obs <- apply(mu, 2, mean)
    tau1_per_obs <- apply(tau1, 2, mean)
    tau2_per_obs <- apply(tau2, 2, mean)


    plot(means_per_obs, pch=19, col=label,#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 2, 1),
         xlab="Observation", ylab=expression(mu), main = "Mean", cex=min(label, 3))#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 1.5,.6))


    plot(tau1_per_obs, pch=19, col=label,#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 2, 1),
         xlab="Observation", ylab="", main = "First Variance", cex=min(label, 3))#felse(label==label[label[which.min(residuals)]]| label==label[label[which.max(residuals)]], 1.5,.6))


    plot(tau2_per_obs, pch=19, col=label,#ifelse(label==label[which.min(residuals)] | label==label[label[which.max(residuals)]], 2, 1),
         xlab="Observation", ylab="", main = "Second Variance",  cex=min(label,3))#ifelse(label==label[label[which.min(residuals)]] | label==label[label[which.max(residuals)]], 1.5,.6))

    return(list(coef = d, number_of_clusters = n_clusters, clusters=label))


  }

  return(list(coef = d))


}

# AutoRegressive model for time series

#' This function defines a design matrix, to estimate the AR(q) model
#'
#' @param X Covariates of interest
#' @param Y Time Series of interest
#' @param q Order of the AutoRegressive model AR(q)
designmatrix.ar <- function(X, Y, q){
  k <- ncol(X)
  # Let's define the new design matrix
  Y_new <- Y[q:length(Y)]
  matrix1 <- NULL
  matrix2 <- NULL
  for (i in 0:(nrow(X)-q - 1)){
    ind1 <- sort((1+i):(i+q), decreasing = T)
    matrix1 <- rbind(matrix1, Y[ind1])

    if (dim(X)[1] != 0){
      matrix2 <- rbind(matrix2, X[q+i+1,])
      design <- cbind(matrix1, matrix2)
    }
    else{
      design <- matrix1    }

  }

  return(design)

}


#' Fitting a Bayesian Nonparametric AutoRegressive model of order q - AR(q)
#'
#' @param Y Time series
#' @param q Order of the AutoRegressive model AR(q)
#' @param X Covariates of interest (by default, the function fits only the autoregressive part)
#' @param initializer list of initial values - can be obtained as the output of the function initializer()
#' @param iter number of iterations to consider
#' @param burn_in number of initial sampled valued to discard
#' @param thin thinning value
#' @param conf_level confidence level of the credible intervals
#' @param cluster boolean value indicating whether the cluster analysis is desired
#' @param sigma parameter to regulate the Pitman-Yor process - default: sigma=0 (Dirichlet process case).
#' @keywords model, ar model, ar, autoregressive model
#' @example armodel <- bnp.ar(Y, qq, initializer, X, iter = 20000, burn_in = 5000,thin = 50, conf_level = 0.05)
#' @export
bnp.ar <- function(Y, q, initializer,X = matrix(0, ncol=0, nrow=0), iter = 10000, burn_in = 5000,
                   thin = 10, conf_level = 0.05, sigma = 0, plot=F){
  n = length(Y)
  k = ncol(X)
  X_design <- designmatrix.ar(X, Y, q)
  Y_design = Y[(q+1):length(Y)]

  samples <- sampler(X_design, Y_design, initializer_, iter = iter, burn_in = burn_in,
                                    thin = thin, log_rate = T, sigma = sigma)
  eps <- read.table("epsilon.csv", header = F, skip = 1,
                    sep = ";")
  m <- read.table("betas.csv", header = F, skip = 1,
                  sep = ";")
  tau1 <- read.table("tau1.csv", header = F, skip = 1,
                     sep = ";")
  tau2 <- read.table("tau2.csv", header = F, skip = 1,
                     sep = ";")
  mu <- read.table("mu.csv", header = F, skip = 1, sep = ";")
  betas_estimates <- apply(m, 2, mean)
  cred_int = bnpResiduals:::credible_intervals(m, level = conf_level)
  d <- cbind(betas_estimates, cred_int)
  colnames(d) <- c("Estimate", paste0(conf_level * 100/2,
                                      "%"), paste0(100 - conf_level * 100/2, "%"))
  if (plot){
    x11()
    plot(density(as.numeric(eps[dim(eps)[1], ])), type = "l",
         col = 1, lty = 4, lwd = 2, main = "Density of the Empirical Residuals")
  }

  return(list(coef = d))

}



# Gaussian Process with Bayesian Nonparametrics residuals

#' Fitting a Gaussian Process with the Residual Structure explicited in Ascolani, Ghidini (2021+)
#'
#' @param X1 first matrix for the computation of the kernel
#' @param X2 second matrix for the computation of the kernel
#' @param l parameter of the kernel
#' @keywords kernel, covariance
#'
#'
calcSigma <- function(X1,X2,l=1) {
  Sigma <- matrix(rep(0, dim(X1)[1]*dim(X2)[1]), nrow=dim(X1)[1])
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*sum(abs(X1[i,]-X2[j,])/l)^2)
      # here you need to substitute with the appropriate k(.,.)
    }
  }
  return(Sigma)
}




#' Fitting a Gaussian Process with the Residual Structure explicited in Ascolani, Ghidini (2021+)
#'
#' @param X design matrix
#' @param Y response variable
#' @param initializer list of initial values - can be obtained as the output of the function initializer()
#' @param iter number of iterations to consider
#' @param burn_in number of initial sampled valued to discard
#' @param thin thinning value
#' @param Xstar grid to predict
#' @param sigma parameter to regulate the Pitman-Yor process - default: sigma=0 (Dirichlet process case).
#' @keywords model, gaussian process, gp
#' @export


bnp.GP <- function(X, Y, initializer, iter, Xstar, burn_in = 0, thin = 1, Verbose=T, log_rate=T, sigma=0){
  n <- nrow(X)
  k <- ncol(X)
  if (is.null(dim(Xstar))){
    Xstar = matrix(Xstar, ncol=1)
  }
  # Write on file
  null_df <- data.frame(NULL)


  write.table(null_df, file = "mu.csv", row.names = F)
  write.table(null_df, file = "tau1.csv", row.names = F)
  write.table(null_df, file = "tau2.csv", row.names = F)
  write.table(null_df, file = "weights.csv", row.names = F)
  write.table(null_df, file = "component.csv", row.names = F)
  write.table(null_df, file = "epsilon.csv", row.names = F)
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
  comp <- initializer$comp
  p <- initializer$p

  eta <- matrix(0, nrow = n, ncol=iter)
  preds <- matrix(0, nrow = nrow(Xstar), ncol=iter)

  for (j in 1:iter){

    #print(j)

    if (Verbose){
      pb$tick()
      Sys.sleep(1 /iter)}

    # Gaussian Process
    if (is.null(dim(X))){
      X <- as.matrix(X, ncol=1)
    }

    # Calculate the covariance matrix
    sigma_bnp <- diag(c(ifelse(comp==1, complete_matrix[,2], complete_matrix[,3])))

    k.xx <- calcSigma(X,X)


    B = solve(sigma_bnp) + solve(k.xx)
    b = solve(sigma_bnp)%*%(Y-comp*complete_matrix[,1])

    # Save eta
    sampled_eta <- mvrnorm(1, solve(B)%*%b, solve(B))
    eta[,j] <- sampled_eta

    # Update Residuals
    epsilon <- Y-sampled_eta

    # Gamma prob
    w <- sampler_gamma_prob(epsilon, complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], log_rate = log_rate)

    #Update Gamma (comp)
    # print(dim(w))
    comp = apply(w,1,function(x){
      return(sample(c(1,-1), size=1, prob = x))
    })

    # Sample mu

    gs <- sampler_triplet_blockedGS(Y, epsilon, comp, p, Z, mu0, sigma0)

    # Labels - they identify the cluster of belonging for each observation
    K <- gs$K

    Z <- gs$Z
    complete_matrix <- gs$complete_matrix

    # Update p - stick breaking probabilities
    p <- p_sampler(K, theta, N, sigma)
    K_unique <- unique(K)

    # Update sigma0 - thanks to semiconjugacy
    sigma0 = sqrt(rinvgamma(n=1, t + length(K_unique)/2, T_big+0.5*sum(Z[K_unique,1]^2)))


    # Compute empirical predictive distribution on a grid of points
    ll <- list(component = comp,
               weights = w, mu = complete_matrix[,1], tau1 = complete_matrix[,2], tau2=complete_matrix[,3], epsilon = epsilon)


    # Predictions on Xstar (formulae by Rasmussen)

    k.xxs <- calcSigma(X,Xstar)
    k.xsx <- calcSigma(Xstar,X)
    k.xsxs <- calcSigma(Xstar,Xstar)

    mu.star <- k.xsx%*%solve(k.xx + sigma_bnp)%*%(Y-comp*complete_matrix[,1])
    sigma.star <- k.xsxs - k.xsx%*%solve(k.xx + sigma_bnp)%*%k.xxs


    preds[,j] <- mvrnorm(1, mu.star, sigma.star)


    if (j >= burn_in & j%%thin ==0){
      # Append to file
      write.table(matrix(ll$component, nrow=1), "component.csv", sep = ";", col.names = !file.exists("component.csv"), append = T, row.names = F)
      write.table(matrix(ll$weights, nrow=1), "weights.csv", sep = ";", col.names = !file.exists("weights.csv"), append = T, row.names = F)
      write.table(matrix(ll$mu, nrow=1), "mu.csv", sep = ";", col.names = !file.exists("mu.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau1, nrow=1), "tau1.csv", sep = ";", col.names = !file.exists("tau1.csv"), append = T, row.names = F)
      write.table(matrix(ll$tau2, nrow=1), "tau2.csv", sep = ";", col.names = !file.exists("tau2.csv"), append = T, row.names = F)
      write.table(matrix(ll$epsilon, nrow=1), "epsilon.csv", sep = ";", col.names = !file.exists("epsilon.csv"), append = T, row.names = F)

    }



  }
  return(list(preds=preds, mu.star = mu.star, sigma.star = sigma.star, Xstar = Xstar,
              eta=eta))

}


# Diagnostic and Inference

#' Credible Intervals
#'
#' This function allows you to compute credible intervals of the linear regression coefficients
#' @param parameters empirical posterior distribution of the parameters of interest
#' @param level credible level (default: level=0.05)
#' @keywords credible intervals, diagnostic

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

#' Coverage
#'
#' This function allows you to compute the coverage of the credible (or confidence) intervals
#' @param real_value true value of the parameters of interest
#' @param intervals credible or confidence intervals
#' @keywords coverage, diagnostic
coverage <- function(real_value, intervals){
  cover <- sum(real_value>=intervals[,1] & real_value<=intervals[,2])/length(beta_true)
  return(cover)
}

#' Model Selection
#'
#' This function allows you to perform model selection
#' @param X design matrix
#' @param Y response variable
#' @param cutoff_level threshold to include a variable (default: cutoff_level = 0.75)
#' @keywords model selection, diagnostic
model_selection <- function(X, Y, sigma2, cutoff_level = 0.75){
  n = nrow(X)
  t_cutoff <- qt(cutoff_level, df=n-1)
  index <- matrix(coef(lm(Y~X[,-1])), nrow=1)%*%(t(X)%*%X)/sigma2
  selected <- index>t_cutoff
  return(X[,selected])
}





