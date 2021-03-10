n <- 100 # length of the time serie
k <- 2 # let's include k covariates

# Simulate the starting data
Y <- rnorm(n)
X <- matrix(rnorm(n*k), ncol=k)



qq = 2 # AR(qq)
# Let's compute the exact data
for (t in qq:(length(Y)-1)){
  Y[t+1] =  .5*Y[t] + .5*Y[t-1] + X[t+1,1] + X[t+1,2]+rnorm(1, 0, 5)#X[t,]%*%delta +Y[(t-count):(t-count+qq-1)]%*%gamma + rnorm(1)#gen_residuals(1, 100, 30, 3) #+
 }


# Define the new design matrix, in order to be able to treat the AR model as a linear one
designmatrix.ar <- function(X, Y, q){
  k <- ncol(X)

  #if (n+q > length(Y)){
  # warning("You are trying to predict something too far in the future")
  # }
  
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
X_design <- designmatrix.ar(X, Y, qq)
Y_design <- Y[(qq+1):length(Y)]

# OLS coef
b <- coef(lm(Y_design~ X_design -1))
b

# Initialization
b0 <-  rep(0, length(b))

B0 <-  diag(10, length(b))

s1 <-  2
s2 <-  2
S1 <- 1#0.01
S2 <-  1#30

sigma2=10

theta = 1
epsilon = Y_design-X_design%*%b # ols residuals

comp <- ifelse(epsilon > 0, 1, -1)

N = 150

log_s = T
V <- rbeta(N-1,shape1 = 1,theta)


W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W

mu0 = 1
sigma0 = sqrt(25)
library(truncnorm)
library(invgamma)
Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
complete_matrix = matrix(1, nrow=nrow(X_design), ncol=3)


initializer_ <- bnpResiduals:::initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)
# model
library(bnpResiduals)


#model <- bnpResiduals:::bnp.lm(X_design, Y_design, initializer_, cluster = F, iter=5000, burn_in = 500, thin=50, sigma=0.8)


bnp.ar <- function(Y, q, initializer,X = matrix(0, ncol=0, nrow=0), iter = 10000, burn_in = 5000, 
                   thin = 10, conf_level = 0.05, sigma = 0, plot=F){
  n = length(Y)
  k = ncol(X)
  X_design <- designmatrix.ar(X, Y, q)  
  Y_design = Y[(q+1):length(Y)]
  
  samples <- bnpResiduals:::sampler(X_design, Y_design, initializer_, iter = iter, burn_in = burn_in, 
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

armodel <- bnp.ar(Y, qq, initializer, X, iter = 20000, burn_in = 5000, 
                  thin = 50, conf_level = 0.05)

armodel


## Diagnostic

eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
m <- read.table("betas.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")


betas_estimates <- apply(m, 2, mean)
betas_estimates

# Empirical Credible Intervals
cred_int = bnpResiduals:::credible_intervals(m)

ind=2
plot(m[,ind], type="l")



