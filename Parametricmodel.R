# Model Y = Xbeta + epsi ---> Y | beta, sigma2 ~ N(Xbeta, sigma2) 
# epsi | SIGMA2 ~ N(0, sigma2)
# beta | beta0 ~ N(beta0, sigma2b)
# beta0 ~ N(0, sigma20)
# sigma2 ~ IG(a, b)




sample_posterior_beta0 <- function(beta, sigma2beta, sigma20){
   k <- length(beta)
   beta0_post <- rnorm(n = 1, mean = mean(beta)*(sigma2beta*sigma20*k)/(k*sigma20+sigma2beta),
                       sd= sqrt((sigma20*sigma2beta)/(k*sigma20+sigma2beta)))
  return(beta0_post)
}



loglik <- function(Y, sigma2, beta, X){
  return(sum(dnorm(Y, X%*%beta, sd = sqrt(sigma2), log = T)))
}

posterior_beta <- function(X, Y, sigma2, sigma2beta, beta0){
  k = ncol(X)
  beta0 = rep(beta0, k)
  B = diag(rep(sigma2beta, k))
  Bp = t(X)%*%X/sigma2+solve(B)
  bp = t(X)%*%Y/sigma2 + solve(B)%*%beta0
  post <- mvrnorm(1, mu = solve(Bp)%*% bp ,  Sigma= solve(Bp))
  return(post)
}


posterior_sigma2 <- function(c0, d0, X, Y, beta){
  N <- nrow(X)
  new_shape <- c0 + N/2
  new_scale <- d0 + .5*sum((Y-X%*%beta)^2)
  sigma2_post <- rinvgamma(1, shape=new_shape, scale = new_scale)
  return(sigma2_post)
}
#### DATA ###
n <- 100
p <- 10
beta_true <- rnorm(n=p, beta0, sigma20)#sample(1:10, replace=T, size=p)

X = matrix(rnorm(n*p), ncol=p, nrow=n)
Y = X%*% beta_true + rnorm(n)


#### HYPERPARS ####
sigma20 <- 1
beta0 <- 0
sigma2beta <- 1
c0 <- 1
d0 <- 1

### INITIALIZATION ###
sigma2 <- 1

### GIBBS SAMPLER

GS <- function(X, Y, beta, beta0, sigma2, sigma2beta, c0, d0, iters){
  N <- nrow(X)
  beta0_post <- c()
  beta_post <- matrix(0, ncol=iters, nrow = p)
  sigma2_post <- c()#matrix(0, ncol=iters, nrow = N)
  
  for (it in 1:iters){
    beta0 <- sample_posterior_beta0(beta, sigma2beta, sigma20)
    beta <- posterior_beta(X, Y, sigma2, sigma2beta, beta0)
    sigma2 <- posterior_sigma2(c0, d0, X, Y, beta)
    
    beta0_post <- c(beta0_post, beta0)
    beta_post[,it] <- beta
    sigma2_post <- c(sigma2_post, sigma2)
  }
 
 return(list(beta0_post = beta0_post, beta_post = beta_post, sigma2_post = sigma2_post))  
}



final <- GS(X, Y, beta_true, beta0, sigma2, sigma2beta, c0, d0, iters=10000)

final$beta0_post
beta_post <- apply(final$beta_post[,5000:10000], 1, mean)
trace <- c()
iters = 10000
for (it in 1:iters){
  tt <- loglik(Y, X = X , beta = final$beta_post[,it], sigma2 = final$sigma2_post[it])
  trace <- c(trace, tt)
}

plot(trace, type="l")
sum(Y-X%*%beta_post)/n
plot(c(Y), c(X%*%beta_post))
abline(a=1,b=1,col=2)

beta_true
beta_post

plot(density(Y-X%*%beta_post))
lines(seq(-3,3, length=100), dnorm(seq(-3,3, length=100)), col=2)

             
