### Competitors - BNP Resisuals

library(mvtnorm)
library(invgamma)
library(truncnorm)
### Skew-T (Sahu et al.)

# Part 1 : linear model
beta_sampler <- function(sigma2, w, delta, z, X, Y, beta0, Lambda0){
  Q_minus1 <- diag(w/sigma2)
  Lambdan <- solve(Lambda0) + t(X)%*%Q_minus1 %*% X
  mun <- solve(Lambdan)%*%(solve(Lambda0)%*% beta0+
                             t(X)%*%Q_minus1%*%Y + delta*t(X)%*%Q_minus1%*%z)

  beta_post <- rmvnorm(n=1, mean=mun, sigma=Lambdan)
  return(matrix(beta_post, ncol=1))
}

#beta = beta_sampler(sigma2, w, delta, z, X, Y, beta0, Lambda0)

sigma2_sampler <- function(beta, X, Y, w, delta, z, a0, b0){
  n = length(Y)
  an = a0 + n/2
  Q_minus1 = diag(w)
  bn = b0 + .5*(t(Y- X%*%beta - delta*z)%*%Q_minus1%*%(Y- X%*%beta - delta*z))

  sigma2_post <- rinvgamma(1, an, bn)
  return(sigma2_post)
}
#sigma2_sampler(beta, X, Y, w, delta, z, a0=1, b0=1)
# Part 2: new

z_sampler <- function(Y, X, beta, delta, sigma2, w){

  A = 1 + w*delta^2/sigma2 # vector nx1
  a = w*delta/sigma2*(Y-X%*%beta)

  z_post <- rtruncnorm(1, mean=a/A, sd=sqrt(1/A), a=0)
  return(z_post)

}

#z_sampler(Y, X, beta, delta, sigma2, w)

delta_sampler <- function(Y, X, beta, z, w, sigma2, gamma2){
  B = 1/gamma2+sum(z^2)/sigma2
  b = sum(z*(Y-X%*%beta))/sigma2

  delta_post <- rnorm(n=1, mean=1/B*b, sd=sqrt(1/B))
  return(delta_post)
}

#delta_sampler(Y, X, beta, z, w, sigma2, gamma2=1)


nu_sampler <- function(w, nu_old, sigma_old, iter, acc_old, acc_star=0.23){
  # MH with adaptive variance
  sigma_new <- exp(log(sigma_old) + iter^(-0.8)*(acc_old - acc_star))
  nu_new <- rnorm(n=1, mean= nu_old, sd=sigma_new)
  if (nu_new <= 2){
    lacc_new <- -Inf
  }
  else{
    lacc_new <- min(1,
                 dnorm(nu_old, mean=nu_new, sd=sigma_new, log = T) +
                   sum(dgamma(w, nu_new/2, nu_new/2, log = T)) +
                   dgamma(nu_new, 1, 0.1) -
                   (pnorm(nu_new, mean=nu_old, sd=sigma_new, log = T) +
                      sum(dgamma(w, nu_old/2, nu_old/2, log = T)) +
                      dgamma(nu_old, 1, 0.1)))
  }

  if (runif(1)< exp(lacc_new)){
    return(list(nu_post=nu_new, sigma_post = sigma_new, acc_post = exp(lacc_new)))
  }

  else{
    return(list(nu_post=nu_old, sigma_post = sigma_old, acc_post = acc_old))
  }
}


# nu_sampler(w, nu_old=30, sigma_old=.00001, iter=10000, acc_old=.5, acc_star=0.23)


sampler_w <- function(X, Y, beta, nu, z, delta){
  n <- length(Y)
  w_post <- rgamma(n, nu/2 + 1/2, 1/2*(nu + (Y-X%*%beta - delta*z)^2/sigma2))
  return(w_post)

}

#sampler_w(X, Y, beta, nu=.8, z, delta)



skewt_sampler <- function(X, Y, iters, sigma2_init, delta_init, w_init, z_init, nu_init,
                          beta0, Lambda0, a0, b0, gamma2
                          ){
  n <- nrow(X)
  p <- ncol(X)
  beta_post <- matrix(0, ncol=iters, nrow=p)
  sigma2_post <- rep(0, iters)
  nu_post <- rep(0, iters)
  z_post <- matrix(0, ncol=iters, nrow=n)
  w_post <- matrix(0, ncol=iters, nrow=n)
  delta_post <- rep(0, iters)
  # Init
  sigma2 <- sigma2_init
  delta <- delta_init
  w <- w_init
  z <- z_init
  nu <- nu_init
  acc_old <- 1
  sigma_old <- 1
  for (t in 1:iters){

    # part 1 - linear model
    beta  <- beta_sampler(sigma2, w, delta, z, X, Y, beta0, Lambda0)
    sigma2 <- sigma2_sampler(beta, X, Y, w, delta, z, a0, b0)

    # part 2 - residual part
    z <- z_sampler(Y, X, beta, delta, sigma2, w)
    delta <- delta_sampler(Y, X, beta, z, w, sigma2, gamma2)
    nu_post <- nu_sampler(w, nu, sigma_old, iter=t, acc_old, acc_star=0.23)
    nu <- nu_post$nu_post
    sigma_old <- nu_post$sigma_post
    acc_old <- nu_post$acc_post
    w <- sampler_w(X, Y, beta, nu, z, delta)


    beta_post[,t] <- beta
    sigma2_post[t] <- sigma2
    z_post[,t] <- z
    delta_post[t] <- delta
    nu_post[t] <- nu
    w_post[,t] <- w
  }

  return(list(beta = beta_post, sigma2=sigma2_post, z=z_post, delta=delta_post,
              nu = nu_post, w = w_post))
}



skewt_sampler(X, Y_gaussian, iters=5000, sigma2_init=1, delta_init=1, w_init=rep(1,n),
              z_init=rep(1,n), nu_init=1,
              beta0=b0, Lambda0=B0, a0=1, b0=1, gamma2=1)



