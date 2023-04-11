#### SImulated data ####


# Base scenario: Linear regression #

n <- 1000
p <- 10

X <- matrix(rnorm(n*p), ncol=p)

# Add the intercept to the design matrix
X <- cbind(rep(1,n), X)
true_beta <- rpois(p+1, 3)

######################################
# SCENARIO I: Gaussian residuals
######################################
true_sigma <- 3
epsilon_gaussian <- rnorm(n, mean=0, sd=true_sigma)

Y_gaussian <- X%*%true_beta + epsilon_gaussian

plot(density(Y_gaussian), main="Gaussian Residuals")
######################################
# SCENARIO II: Heavy-tailed residuals
######################################

df_t <- 1.1
epsilon_heavytailed <- rt(n, df_t)

Y_heavytailed <- X%*%true_beta + epsilon_heavytailed

plot(density(Y_heavytailed), main="Heavy-Tailed Residuals")

######################################
# SCENARIO III: asymmetric residuals (skew normal)
######################################
library(sn)

alpha = 4
epsilon_skewnormal <- rsn(n=n, xi=0, omega=1, alpha)
Y_skewnormal <- X%*%true_beta + epsilon_skewnormal

plot(density(Y_skewnormal), main="Skew-normal Residuals")

######################################
# SCENARIO IV: asymmetric residuals (mixture model)
######################################


p1 <- .5
p2 <- .5

mu1 <- -3
sigma1 <- 1


mu2 <- 100
sigma2 <- 30


epsilon_mixture <- c()

for (i in 1:n){
  if (runif(1)<p1){
    epsilon_mixture <- c(epsilon_mixture, rnorm(1, -mu2, sigma1))
  }
  else{
    epsilon_mixture <- c(epsilon_mixture, rnorm(1, mu2, sigma2))
  }
  
}
Y_mixture <- X%*%true_beta + epsilon_mixture

plot(density(Y_mixture), main="Mixture Residuals")

######################################
# SCENARIO V: outliers
######################################

mean_outlier = 100

epsilon_outlier <- epsilon_gaussian
epsilon_outlier[sample(1:n, size=n/20, replace=F)] <- rnorm(n/20, mean_outlier, sd=sqrt(mean_outlier))

Y_outlier <- X%*%true_beta + epsilon_outlier

plot(density(Y_outlier), main="Residuals - with outlier")


######################################
# SCENARIO VI: SkewT residuals
######################################

epsilon_skewT <- rst(n)

Y_skewT <- X%*%true_beta + epsilon_skewT

plot(density(Y_skewT), main="Skew-T Residuals")

#############################
######## MODELS #############
#############################



setwd("C:\\Users\\valen\\Desktop\\bnpResiduals\\Competitors")



#########################################
### MODEL 1: our model, but symmetric ###
#########################################
source("Competitor_symmetricmodel.R")


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

b <- coef(lm(Y_mixture~X[,-1]))

epsilon = Y_mixture-X%*%b # ols residuals

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
library(truncnorm)
library(invgamma)


K <- sample(1:N, size = n, replace = T, prob = p)
Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
temp <- sampler_triplet_blockedGS(Y_mixture,epsilon, comp, p, Z, mu0, sigma0)
K <- temp$K
Z <- temp$Z
complete_matrix = Z[K,]



#attach(bnpResiduals)
initializer_ <- initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)

symmetric_model <- bnp.lm_symmetric(X, Y_mixture, initializer_, iter=5000, burn_in = 500, thin = 10, cluster = T)

symmetric_model$coef

#########################################
### MODEL 2: Standard Linear Model ###
#########################################

source("LinearModel.R")

linear_model <- linear_model(X, Y_mixture, iters=5000, Lambda0=B0, mu0=b0, a0=1, b0=1)

apply(linear_model$beta, 1, mean)
true_beta


#########################################
### MODEL 3: Skew-T model ###
#########################################

source("skewT.R")
skewt_model <- skewt_sampler(X, Y_skewT, iters=5000, sigma2_init=1, delta_init=1, w_init=rep(1,n), 
                          z_init=rep(1,n), nu_init=1,
                          beta0=b0, Lambda0=B0, a0=1, b0=1, gamma2=1)




dim(skewt_model$beta)

apply(skewt_model$beta, 1, mean)
true_beta



