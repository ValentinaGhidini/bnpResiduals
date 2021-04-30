require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)
library(bnpResiduals)



X <- matrix(seq(-100, 100), by=.1, ncol=1)
n = length(X)
Y = sin(X) + rnorm(n, 0, 0.2)

## Initialization
n <- nrow(X)
k <- ncol(X)

b0 <-  rep(0, k+1)

B0 <-  diag(10, k+1)

s1 <-  2
s2 <-  2
S1 <- 1#0.01
S2 <-  1#30

sigma2=10

theta = 1

#b <- coef(lm(Y~X))

epsilon = rnorm(n)#Y-X%*%b # ols residuals

comp <- ifelse(epsilon > 0, 1, -1)

N = 200

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
complete_matrix = matrix(1, nrow=n, ncol=3)


initializer_ <- bnpResiduals:::initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)


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


#X <- matrix(seq(-20, 20), by=.1, ncol=1)
n = length(X)

#Y = X+ rnorm(n)
#x.star <- apply(X, 2, function(x) {seq(min(x)-.5, max(x)+.5, len=n.samples)}) # infamous grid
x.star <- seq(-20, 20, by = 0.5)

bnp.GP <- function(X, Y, initializer, iter, x.star, burn_in = 0, thin = 1, Verbose=T, log_rate=T, sigma=0){
  n <- nrow(X)
  k <- ncol(X)
  if (is.null(dim(x.star))){
    x.star = matrix(x.star, col=1)
  }
  # Write on file
  null_df <- data.frame(NULL)
  Grid <- seq(-50, 50, by=.1)
  predictive <- matrix(0, ncol=length(Grid), nrow=iter)
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
  preds <- matrix(0, nrow = nrow(x.star), ncol=iter)

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

    
    # Predictions on x.star (formulae by Rasmussen)

    k.xxs <- calcSigma(X,x.star)
    k.xsx <- calcSigma(x.star,X)
    k.xsxs <- calcSigma(x.star,x.star)
    
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
  return(list(preds=preds, mu.star = mu.star, sigma.star = sigma.star, x.star = x.star, 
              eta=eta))
  
}
#x.star <- seq(-5, 5,len=50) # infamous grid
x.star <- matrix(x.star, ncol=1)
ss <- bnp.GP(X, Y, initializer_, x.star, iter=5000)
values <- cbind(x = x.star, ss$preds)
prediction = apply(ss$preds, 1, mean)
head(values)
values <- melt(values,id="x")

df <- data.frame(cbind(x.star = x.star[,1], f.star.bar=prediction))
colnames(df) <- c("x.star", "f.star.bar")
figg <-  ggplot(values,aes(x=Var1,y=value)) +
  #geom_line(aes(group=Var1), colour="grey80") +
  geom_line(data=data.frame(X, Y), aes(x=X, y=Y), colour="blue")+
  geom_line(data=df,aes(x=x.star,y=prediction), colour="red", size=1) + 
  geom_point(data=data.frame(x=X, y=Y),aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(min(Y)-1,max(Y)+1), name="output, f(x)") +
  xlab("input, x")
figg



eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui Gaussiani\\Results\\
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui Mistura di T\\Results\\
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui T1\\Results\\
#m <- read.table("betas.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

head(mu)

plot(mu[,10], type="l")
plot(tau1[,40], type="l")
plot(tau2[,40], type="l")

plot(apply(ss$eta, 1, mean), pch=19)
