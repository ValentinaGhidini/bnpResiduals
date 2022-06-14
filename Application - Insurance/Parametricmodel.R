library(orthopolynom)
library(MASS)
library(invgamma)
# Model Y = Xbeta + epsi ---> Y | beta, sigma2 ~ N(Xbeta, sigma2) 
# epsi | SIGMA2 ~ N(0, sigma2)
# beta | beta0 ~ N(beta0, sigma2b)
# beta0 ~ N(0, sigma20)
# sigma2 ~ IG(a, b)

######### NEW MODEL ##################

posterior_sigma2 <- function(c0, d0, X, Y, beta, groups){
  N <- nrow(X)
  new_shape <- c0 + N/2
  residuals <- rep(0, N)
  L <- length(unique(groups))
  for (l in 1:L){
    residuals[groups==l] <- Y[groups==l] - X[groups==l,] %*% beta[,l]
  }
  new_scale <- d0 + .5*sum(residuals^2)
  sigma2_post <- rinvgamma(1, shape=new_shape, scale = new_scale)
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

######## SIMULATION ######

#### HYPERPARS ####
p <- 5
# beta0 ~ N(b0, B0)
b0 <- rep(0, p)
B0 <- matrix(0.1, ncol=p, nrow=p)
diag(B0) <- 1


# beta|beta0 ~ N(beta0, B)
B <- matrix(0, ncol=p, nrow=p)
diag(B) <- 1

# sigma2 ~ IG(c0, d0)
c0 <- 1
d0 <- 1


### INITIALIZATION ###
sigma2 <- 1

#### DATA ###
n <- 100


GRUPPI = 4

beta0 <- mvrnorm(1, b0, B0)
beta_true <- matrix(0, nrow=p, ncol=GRUPPI)#sample(1:10, replace=T, size=p)
for (g in 1:GRUPPI) {
  beta_true[,g]<- mvrnorm(1, beta0, B)
}

X = matrix(rnorm(n*p), ncol=p, nrow=n)
groups <- sample (1:GRUPPI, size=n, replace=T)

Y <- rep(0,n)
for (g in 1:GRUPPI) {
  Y[groups==g] = X[groups==g,]%*% beta_true[,g] + rnorm(sum(groups==g))
}



### GIBBS SAMPLER
GS_new <- function(X, Y, beta0, B0,  B, c0, d0, sigma2, iters, groups){
  # L = number of groups
  # gruppo = gruppo di appartenenza
  if (is.null(dim(Y))){ Y <- matrix(Y, ncol=1)}
  N <- nrow(X)
  beta0_post <- matrix(0, ncol=iters, nrow = p)
  beta_post_full <- list()#matrix(0, ncol=iters, nrow = p)
  sigma2_post <- c()#matrix(0, ncol=iters, nrow = N)
  L <- length(unique(groups))
  
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


final<- GS_new(X, Y, beta0, B0,  B, c0, d0, sigma2, iters=10000, groups=groups)


final$beta_post_full[[10000]]
beta_true

final_estimates_MAP <- matrix(0, nrow=p, ncol=GRUPPI)
for (el in final$beta_post_full){
  for (i in 1:nrow(el)){
    for (j in 1:ncol(el)){
      final_estimates_MAP[i, j] <- final_estimates_MAP[i, j] + el[i,j]
    }
  }
}
final_estimates_MAP <- final_estimates_MAP/10000
final_estimates_MAP
beta_true
#beta_post <- apply(final$beta_post[,5000:10000], 1, mean)

# Traceplot
trace <- c()
iters = 10000
for (it in 1:iters){
  tt <- loglik(Y, X = X , beta = final$beta_post_full[[it]], sigma2 = final$sigma2_post[it], groups=groups)
  trace <- c(trace, tt)
}

plot(trace, type="l", main="Traceplot")

# predictions
preds <- rep(0, n)

for (l in 1:GRUPPI){
  preds[groups==l] <-  X[groups==l,]%*%final_estimates_MAP[,l]
}

# Residuals
plot(c(Y), preds, pch=19, main="Residuals")
abline(a=0,b=1,col=2)

# Density of residuals
plot(density(Y-preds), main="Density of Residuals")
lines(seq(-3,3, length=100), dnorm(seq(-3,3, length=100)), col=2)

             
