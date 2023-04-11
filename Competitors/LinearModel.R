library(mvtnorm)
library(invgamma)
library(truncnorm)
## Linear model


sampler_beta_linearmodel <- function(X, Y, Lambda0, mu0, sigma2){
  Lambdan = t(X)%*%X + Lambda0
  mun = solve(Lambdan)%*%(t(X)%*%Y + Lambda0 %*%mu0)
  beta_post <- rmvnorm(1, mean= mun, sigma= sigma2*solve(Lambdan))
  return(beta_post)
}



sampler_sigma2_linearmodel <- function(X, Y, Lambda0, mu0, a0, b0){
  n = nrow(X)
  Lambdan = t(X)%*%X + Lambda0
  mun = solve(Lambdan)%*%(t(X)%*%Y + Lambda0 %*%mu0)
  an = a0 + n/2
  bn = b0 + 0.5*(t(Y)%*%Y + t(mu0)%*%Lambda0%*%mu0 - t(mun)%*%Lambdan%*%mun)
  sigma2_post <- rinvgamma(1, an, bn)
  return(sigma2_post)
}




linear_model <- function(X, Y, iters, Lambda0, mu0, a0, b0){
  p = ncol(X)
  sigma2_post <- rep(0, iters)
  beta_post <- matrix(0, nrow=p, ncol=iters)
  
  for (it in 1:iters){
    beta <- sampler_beta_linearmodel(X, Y, Lambda0, mu0, sigma2)
    sigma2 <- sampler_sigma2_linearmodel(X, Y, Lambda0, mu0, a0, b0)
    sigma2_post[it] <- sigma2
    beta_post[,it] <- beta
  }
  
  return(list(beta=beta_post, sigma2=sigma2_post))
}
