## Final Application BNP ##
rm(list=ls())
library(orthopolynom)
library(MASS)
library(invgamma)
library(truncnorm)
library(invgamma)
library(bnpResiduals)
library(ggplot2)
library(reshape2)

# Set working directory
setwd("C:\\Users\\valen\\Desktop\\Misture di Misture")
source("code complete/bnp_residuals_blockedGS.R")


######### PARAMETRIC MODEL ##################

posterior_sigma2 <- function(c0, d0, X, Y, beta, groups){
  N <- nrow(X)
  new_shape <- c0 + N/2
  residuals <- rep(0, N)
  L <- length(unique(groups))
  for (l in 1:L){
    residuals[groups==l] <- Y[groups==l] - X[groups==l,] %*% beta[,l]
  }
  new_scale <- d0 + .5*sum(residuals^2)
  sigma2_post <- rinvgamma(1, shape=new_shape, rate = new_scale)
  return(sigma2_post)
}


sample_posterior_beta0 <- function(beta, b0, B, B0){
  if (is.null(dim(b0))){b0 <- matrix(b0, ncol=1)}
  if (is.null(dim(beta))){beta <- matrix(beta, ncol=1)}
  mean_vector <- apply(beta, 1, mean)
  L <- ncol(beta)
  k <- length(beta)
  new_mean <- solve(L*solve(B) + solve(B0))%*%(L*solve(B)%*%mean_vector+solve(B0)%*%b0)
  new_var <- solve(solve(B) + solve(B0))
  beta0_post <- mvrnorm(1, mu= new_mean, Sigma = new_var)
  return(beta0_post)
}



posterior_beta <- function(X, Y, sigma2, B, beta0){
  if (is.null(dim(b0))){b0 <- matrix(b0, ncol=1)}
  if (is.null(dim(Y))){ Y <- matrix(Y, ncol=1)}
  
  Bp = t(X)%*%X/sigma2+solve(B)
  bp = t(X)%*%Y/sigma2 + solve(B)%*%beta0
  post <- mvrnorm(1, mu = solve(Bp)%*% bp ,  Sigma= solve(Bp))
  return(post)
}


loglik <- function(Y, sigma2, beta, X, groups){
  L <- length(unique(groups))
  llik <- 0
  for (l in 1:L){
    llik <- llik + sum(dnorm(Y[groups==l], X[groups==l,]%*%beta[,l], sd = sqrt(sigma2), log = T))
  }
  return(llik)
}


GS_parametric <- function(X, Y, b0, B0,  B, c0, d0, sigma2, iters, groups){
  # L = number of groups
  # gruppo = gruppo di appartenenza
  p <- ncol(X)
  if (is.null(dim(Y))){ Y <- matrix(Y, ncol=1)}
  N <- nrow(X)
  beta0_post <- matrix(0, ncol=iters, nrow = p)
  beta_post_full <- list()#matrix(0, ncol=iters, nrow = p)
  sigma2_post <- c()#matrix(0, ncol=iters, nrow = N)
  L <- length(unique(groups))
  beta0 <- rep(0,p)
  for (it in 1:iters){
    beta_post <- matrix(0, ncol=L, nrow = p)
    
    # sample betas for each group
    for (l in 1:L){
      
      beta_post[,l] <- posterior_beta(X[groups==l, ], Y[groups==l], sigma2, B, beta0)
      
    }
    beta0 <- sample_posterior_beta0(beta_post, b0, B, B0)
    sigma2 <- posterior_sigma2(c0, d0, X, Y, beta_post, groups )
    
    beta0_post[,it] <- beta0
    beta_post_full[[it]] <- beta_post
    sigma2_post <- c(sigma2_post, sigma2)
  }
  
  return(list(beta0_post = beta0_post, beta_post_full = beta_post_full, sigma2_post = sigma2_post))  
}


#### NONPARAMETRIC MODEL ####


GS_nonparametric <- function(X, Y, groups, initializer, iter, Verbose=T, 
                              log_rate=T, sigma=0, c0=1, d0=1){
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
  
  
  # Cluster analysis for outliers
  clusters <- matrix(0, nrow = iter, ncol=n)
  
  B <- matrix(0, ncol=ncol(X), nrow=ncol(X))
  diag(B) <- 1
  beta_post_full <- list()
  ZZ <- list()
  sigma2_post <- c()
  p_weights <- matrix(0, nrow=N, ncol=iter)
  for (j in 1:iter){
    if (Verbose){
      pb$tick()
      Sys.sleep(1 /iter)}
    
    L <- length(unique(groups))
    
    beta_post <- matrix(0, ncol=L, nrow = ncol(X))
    beta0_post <- matrix(0, ncol=iter, nrow = p)
    # sample betas for each group
    
    for (l in 1:L){
      
      ss <- sampler_beta(complete_matrix[groups==l,1], complete_matrix[groups==l,2], complete_matrix[groups==l,3], 
                         comp[groups==l], Y[groups==l], X[groups==l,], B0, matrix(b0, ncol=1))
      
      beta_post[,l] <- ss$beta
      #print("Finished betas")
      # Update Residuals
      epsilon[groups==l] <- Y[groups==l]-X[groups==l,]%*%beta_post[,l]
      #print("Finished epsilon")
    }
    beta0 <- sample_posterior_beta0(beta_post, b0, B, B0)
    sigma2 <- posterior_sigma2(c0, d0, X, Y, beta_post, groups)
    
    beta0_post[,j] <- beta0
    beta_post_full[[j]] <- beta_post
    sigma2_post <- c(sigma2_post, sigma2)
    # Sample beta - for each group
    #ss <- sampler_beta(complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], comp, Y, X, B0, b0)
    #beta <- ss$beta
    
    
    
    # Gamma prob
    w <- bnpResiduals:::sampler_gamma_prob(epsilon, complete_matrix[,1], complete_matrix[,2], complete_matrix[,3], log_rate = log_rate)
    
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
    ZZ[[j]] <- Z
    complete_matrix <- gs$complete_matrix
    
    # Update p - stick breaking probabilities
    p <- bnpResiduals:::p_sampler(K, theta, N, sigma)
    
    p_weights[,j] <- p
    
    K_unique <- unique(K)
    
    # Update sigma0 - thanks to semiconjugacy
    sigma0 = sqrt(rinvgamma(n=1, t + length(K_unique)/2, T_big+0.5*sum(Z[K_unique,1]^2)))
    
    
    # Compute empirical predictive distribution on a grid of points
    predictive[j,] = sapply(Grid, function(x) sum(p*(.5*dnorm(x, Z[,1], sd=Z[,2])+.5*dnorm(x, -Z[,1], sd=Z[,3]))))
    ll <- list(beta = beta_post_full, component = comp, beta0 = beta0_post,
               weights = w, mu = complete_matrix[,1], tau1 = complete_matrix[,2], tau2=complete_matrix[,3], epsilon = epsilon)
    
    
    # Append to file
    #write.table(matrix(ll$beta, nrow=1), "betas.csv", sep = ";", col.names = !file.exists("betas.csv"), append = T, row.names = F)
    
    write.table(matrix(ll$component, nrow=1), "component.csv", sep = ";", col.names = !file.exists("component.csv"), append = T, row.names = F)
    write.table(matrix(ll$weights, nrow=1), "weights.csv", sep = ";", col.names = !file.exists("weights.csv"), append = T, row.names = F)
    write.table(matrix(ll$mu, nrow=1), "mu.csv", sep = ";", col.names = !file.exists("mu.csv"), append = T, row.names = F)
    write.table(matrix(ll$tau1, nrow=1), "tau1.csv", sep = ";", col.names = !file.exists("tau1.csv"), append = T, row.names = F)
    write.table(matrix(ll$tau2, nrow=1), "tau2.csv", sep = ";", col.names = !file.exists("tau2.csv"), append = T, row.names = F)
    write.table(matrix(ll$epsilon, nrow=1), "epsilon.csv", sep = ";", col.names = !file.exists("epsilon.csv"), append = T, row.names = F)
    write.table(matrix(ll$beta0, nrow=1), "beta0.csv", sep = ";", col.names = !file.exists("beta0.csv"), append = T, row.names = F)
    
    sigma0_saved = c(sigma0_saved, sigma0)
    
  }
  saveRDS(ll$beta, "betas.RData")
  return(list(predictive=predictive, sigma0=sigma0_saved, clusters=clusters,
              sigma2 = sigma2_post, betas = ll$beta, p = p_weights, Z = ZZ))
}



##### Let's try it on the data #####

# Data
insurance <- read.csv("C:/Users/valen/Desktop/insurance.csv")
groups <- factor(insurance$region)
levels(groups) <- c(1,2,3,4) 
insurance <- insurance[,-6]
insurance$sex <- ifelse(insurance$sex == 'female', 1, 0) # 1 female, 0 male
insurance$smoker <- ifelse(insurance$smoker == 'yes', 1, 0)


Y <- log(insurance$charges)
X <- as.matrix(insurance[,1:5])
X <- cbind(rep(1, nrow(X)), X)


# Initialization + hyperpars for the BNP Model
n <- nrow(X)
k <- ncol(X)-1

b0 <-  rep(0, k+1)

B0 <-  diag(10, k+1)
s1 <-  2
s2 <-  2
S1 <- 1#0.01
S2 <-  1#30

sigma2=10

theta = 1

ols_model <- lm(Y~X[,-1])
b <- coef(ols_model)

epsilon = Y-X%*%b # ols residuals

comp <- ifelse(epsilon > 0, 1, -1)

N = 150

log_s = T
V <- rbeta(N-1,shape1 = 1,theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)


p <- rep(1/N, N)#V*W

mu0 = 1
sigma0 = sqrt(25)



K <- sample(1:N, size = n, replace = T, prob = p)
Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
temp <- sampler_triplet_blockedGS(Y,epsilon, comp, p, Z, mu0, sigma0)
K <- temp$K
Z <- temp$Z
complete_matrix = Z[K,]


initializer_ <- initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)

ITERS <- 10000

#### Standard OLS Model ###

insurance$sex <- relevel(insurance$sex, "male")
ols_model <- lm(log(charges)~ ., data=insurance[,-6])
summary(ols_model)


### Parametric Model ###


parametric_result<- GS_parametric(X, Y, b0=rnorm(k+1), B0=diag(1, ncol=k+1, nrow=k+1),  
                                  B=diag(2, ncol=k+1, nrow=k+1), 
               c0=1, d0=1, sigma2=4, iters=ITERS, groups=groups)

# Coefficients
parametric_estimates <- matrix(0, nrow=k+1, ncol=length(unique(groups)))
for (el in parametric_result$beta_post_full){
  for (i in 1:nrow(el)){
    for (j in 1:ncol(el)){
      parametric_estimates[i, j] <- parametric_estimates[i, j] + el[i,j]
    }
  }
}
parametric_estimates <- parametric_estimates/ITERS
parametric_estimates

# Predictions
preds_parametric <- rep(0, nrow(X))

for (l in 1:length(unique(groups))){
  preds_parametric[groups==l] <- X[groups==l,]%*%parametric_estimates[,l]
}
plot(Y, preds_parametric, pch=19, main="Residuals")
abline(a=0,b=1,col=2)


# Traceplot
trace <- c()

for (it in 1:ITERS){
  tt <- loglik(Y, X = X , beta = parametric_result$beta_post_full[[it]], 
               sigma2 = parametric_result$sigma2_post[it], groups=groups)
  trace <- c(trace, tt)
}

plot(trace, type="l", main="Traceplot")


### Nonparametric model ###

nonparametric_result <- GS_nonparametric(X, Y, groups=groups, initializer_, 
                                         iter=ITERS, 
                         Verbose=T, log_rate=T, sigma=0,  c0=1, d0=1)


# Coefficients

nonparametric_estimates <- matrix(0, nrow=ncol(X), ncol=length(unique(groups)))

for (el in nonparametric_result$betas){
  for (i in 1:nrow(el)){
    for (j in 1:ncol(el)){
      nonparametric_estimates[i, j] <- nonparametric_estimates[i, j] + el[i,j]
    }
  }
}
nonparametric_estimates <- nonparametric_estimates/ITERS
nonparametric_estimates

# Let's read the results
eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")
m <- nonparametric_result$betas

# Traceplots
to_study = c(3,2) # third par of group 2

beta_est = c()
for (el in m){
  beta_est <- c(beta_est, el[to_study[1], to_study[2]])
}
# Autocorrelation
acf(beta_est)
# Traceplots
plot(beta_est, type="l")
abline(h=beta_true[to_study], col=2)
abline(h=mean(beta_est, na.rm=T), col=3, lty=2)



############# Comparison - three methods ##############

# 1) Coefficients
credible_intervals_pergroup <- list()
beta_est <- list()
coefficients <- data.frame(betaOLS = coef(ols_model),
                           lowerOLS = confint(ols_model)[ ,1], upperOLS = confint(ols_model)[ ,2])

for(l in 1:ncol(nonparametric_estimates)){
  for (i in 1:length(m)){
    el <- m[[i]]
    if (i == 1){
      beta_est[[l]] <- el[,l]
    }
    else{
      beta_est[[l]] <- rbind(beta_est[[l]], el[,l])
    }
  }
  credible_intervals_pergroup[[l]] <- credible_intervals(beta_est[[l]])
  
}




for(l in 1:ncol(nonparametric_estimates)){
  coefficients[paste("lowerEstimate_Group",l, sep="" )] = credible_intervals_pergroup[[l]][,1]
  coefficients[paste("upperEstimate_Group",l, sep="" )] = credible_intervals_pergroup[[l]][,2]
}

# 2) Residuals

preds_nonparametric <- rep(0, nrow(X))

for (l in 1:length(unique(groups))){
  preds_nonparametric[groups==l] <- X[groups==l,]%*%nonparametric_estimates[,l]
}



# let's plot

x11()
p <- ggplot(coefficients, aes(x = 1:length(betaOLS)))  +
  
  xlab("") + ylab("Coefficients") +
  ggtitle("Proposed VS OLS Estimates") + 
  geom_line(aes(x=1:length(betaOLS), y = nonparametric_estimates[,1]) , col=1) + 
  geom_point(aes(x=1:length(betaOLS), y = nonparametric_estimates[,1] ) , col=1) +
  geom_ribbon(aes(x = c(1:length(betaOLS)), ymin = credible_intervals_pergroup[[1]][,1],
                  ymax = credible_intervals_pergroup[[1]][,2]), alpha = 0.2, fill=1) + 
  
  geom_line(aes(x=1:length(betaOLS), y = nonparametric_estimates[,2]) , col=2) + 
  geom_point(aes(x=1:length(betaOLS), y = nonparametric_estimates[,2] ) , col=2) +
  geom_ribbon(aes(x = c(1:length(betaOLS)), ymin = credible_intervals_pergroup[[2]][,1],
                  ymax = credible_intervals_pergroup[[2]][,2]), alpha = 0.2, fill=2)+
  
  geom_line(aes(x=1:length(betaOLS), y = nonparametric_estimates[,3]) , col=3) + 
  geom_point(aes(x=1:length(betaOLS), y = nonparametric_estimates[,3] ) , col=3) +
  geom_ribbon(aes(x = c(1:length(betaOLS)), ymin = credible_intervals_pergroup[[3]][,1],
                  ymax = credible_intervals_pergroup[[3]][,2]), alpha = 0.2, fill=3) +
  
  geom_line(aes(x=1:length(betaOLS), y = nonparametric_estimates[,4]) , col=4) + 
  geom_point(aes(x=1:length(betaOLS), y = nonparametric_estimates[,4] ) , col=4) +
  geom_ribbon(aes(x = c(1:length(betaOLS)), ymin = credible_intervals_pergroup[[4]][,1],
                  ymax = credible_intervals_pergroup[[4]][,2]), alpha = 0.2, fill=4) +
  xlim(2,6) + ylim(0, 2)
p


epsi <- apply(eps,2, mean)


## Empirical residuals
x_empirical <- data.frame(proposedResiduals = Y - preds_nonparametric,
                          parametricResiduals = Y - preds_parametric,#as.numeric(eps[dim(eps)[1],]),
                          OLSresiduals=ols_model$residuals)


data_empirical<- melt(x_empirical)
x11()
emp <- ggplot(data_empirical,aes(x = value, fill=variable)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Residuals", labels = c("Proposed", "Estimated", "Parametric", "OLS")) + 
  ggtitle("Residual densities (empirical)") + xlab("Residual") + ylab("Density") 

emp


## Theoretical results
Grid <-seq(-2, 2, length = 200)

dens_proposedResiduals <- matrix(0, ncol=ITERS, nrow=length(Grid))
dens_parametricResiduals <- matrix(0, ncol=ITERS, nrow=length(Grid))

for (i in 1:length(Grid)){
  for (j in 1:ITERS){
    mus <- nonparametric_result$Z[[j]][,1]
    tau1s <- nonparametric_result$Z[[j]][,2]
    tau2s <- nonparametric_result$Z[[j]][,3]
    ww <- nonparametric_result$p[,j]
    dens_proposedResiduals[i,j] <- sum(ww*(0.5*dnorm(Grid[i], mus, tau1s) + 0.5*dnorm(Grid[i], -mus, tau2s)))
    dens_parametricResiduals[i,j] <- dnorm(Grid[i], 0, sqrt(parametric_result$sigma2_post[j]))
  }
}

mean_density_proposal <- apply(dens_proposedResiduals, 1, mean)
mean_density_parametric <- apply(dens_parametricResiduals, 1, mean)
mean_density_ols <- dnorm(Grid, mean=0, sd=summary(ols_model)$sigma)


x_theoretical <- data.frame(proposedResiduals = mean_density_proposal, x =Grid,
                parametricResiduals =mean_density_parametric,
                OLSresiduals=mean_density_ols)


data_theor<- melt(x_theoretical, id.var="x")
x11()

theor <- ggplot(data,aes(x= x, y = value, col=variable)) + 
  geom_line() +
  geom_area(aes(fill = variable), alpha = 0.3, position = 'identity') +
  ggtitle("Residual densities (theoretical)") + xlab("Residual") + ylab("Density") +
  scale_fill_manual(name = 'Residuals', values = c("deepskyblue4","darkviolet", "springgreen4" ),
                    labels = c("Proposed", "Parametric", "OLS")) +
  scale_color_manual(name = 'Residuals', values = c("deepskyblue4","darkviolet", "springgreen4"),
                     labels = c("Proposed", "Parametric", "OLS"))

theor



### Study of clusters
p_last <- nonparametric_result$p[,dim(nonparametric_result$p)[2]]
z_last <- nonparametric_result$Z[[length(nonparametric_result$Z)]]

bigger <- order(p_last, decreasing = T)
z_last[bigger[1:10],]
p_last[bigger[1:10]]


# Let's study clusters 2 and 3
cluster2 <- nonparametric_result$clusters[ITERS,] == bigger[2]
cluster3 <- nonparametric_result$clusters[ITERS,] == bigger[3]
cluster4 <- nonparametric_result$clusters[ITERS,] == bigger[4]
cluster5 <- nonparametric_result$clusters[ITERS,] == bigger[5]


plot(Y, Y, col=cluster3+1, pch=19, cex = 0.2 + cluster3)
plot(Y, Y, col=cluster2+1, pch=19, cex = 0.2 + cluster2)
plot(Y, Y, col=cluster4+1, pch=19, cex = 0.2 + cluster4) # good
plot(Y, Y, col=cluster5+1, pch=19, cex = 0.2 + cluster5)
x11()
par(mfrow=c(2,5))
for (k in 1:10){
  cluster5 <- nonparametric_result$clusters[ITERS,] == bigger[k]
  plot(Y, Y, col=cluster5+1, pch=19, cex = 0.2 + cluster5, 
       main = paste("Cluster: ", k))
}

summary(insurance[cluster2, ])

summary(insurance[cluster3, ])

summary(insurance[cluster4, ])


x11()
par(mfrow=c(2,5))
for (k in 1:10){
  cluster5 <- nonparametric_result$clusters[ITERS,] == bigger[k]
  plot(insurance$bmi, Y, col=cluster5+1, pch=19, cex = 0.2 + cluster5, 
       main = paste("Cluster: ", k))
}
x11()
par(mfrow=c(2,5))
for (k in 1:10){
  cluster5 <- nonparametric_result$clusters[ITERS,] == bigger[k]
  plot(insurance$bmi, Y, col=cluster5+1, pch=19, cex = 0.2 + cluster5, 
       main = paste("Cluster: ", k), type=".")
}

library(GGally)

ggpairs(insurance, title="Insurance data - Cluster 4",
        mapping=ggplot2::aes(colour = cluster4))

x11()
hist(exp(Y))
summary(Y)

Y[Y>quantile(Y,.9)]

outliers_right <- which(Y> mean(Y) + 1.96*sd(Y))
outliers_left <- which(Y< mean(Y) - 1.96*sd(Y))

sort(nonparametric_result$clusters[ITERS,outliers_left])
sort(nonparametric_result$clusters[ITERS,outliers_right])


# Residuals
x11()
plot(x_empirical$proposedResiduals, col=1+cluster2, pch=19)

left_tail_residuals <- which(x_empirical$proposedResiduals<quantile(x_empirical$proposedResiduals, .01))
right_tail_residuals <- which(x_empirical$proposedResiduals>quantile(x_empirical$proposedResiduals, .99))

sort(nonparametric_result$clusters[ITERS,left_tail_residuals])
sort(nonparametric_result$clusters[ITERS,right_tail_residuals])

plot(x_empirical$proposedResiduals, col=1+(x_empirical$proposedResiduals>quantile(x_empirical$proposedResiduals, .99)), 
pch=19)


plot(x_empirical$proposedResiduals, col=1+(x_empirical$proposedResiduals>quantile(x_empirical$proposedResiduals, .99)), 
     pch=19)

plot(x_empirical$proposedResiduals, col=1+(x_empirical$proposedResiduals>quantile(x_empirical$proposedResiduals, .99)), 
     pch=19)
