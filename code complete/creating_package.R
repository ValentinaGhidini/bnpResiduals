library(devtools)
library(roxygen2)
setwd("C:\\Users\\valen\\Desktop")
create("bnpResiduals")

setwd("./bnpResiduals")
document()

setwd("..")
install("bnpResiduals")
library(bnpResiduals)
?initializer
?p_sampler
?sampler_beta
?sampler
?coverage


## Trying analysis with package functions

n <- 100
k <- 9
beta = c(rep(-1, 3), rep(0,1), rep(1,6))
beta_true = beta


#Generate residuals
muu = 100
tau1 = 2
tau2 = 3

gen_residuals <- function(n, mu, tau1, tau2){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(T,F),n, replace = T)
  res[bol] <- rnorm(sum(bol),mu,sqrt(tau1))
  res[!bol] <- rnorm(n-sum(bol),-mu, sqrt(tau2))
  return(res)
}

gen_residuals_T <- function(n, df1, df2, muu){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(T,F),n, replace = T)
  res[bol] <- rt(sum(bol),df1)+muu
  res[!bol] <- rt(n-sum(bol),df2)-muu
  return(res)
}

#dat <- rnorm(n)
#dat <- gen_residuals_T(n,tau1,tau2, muu)
dat <- rt(n, df=2.1)
#dat <- gen_residuals_T(n, 2, 3, 100)
# Density of the residuals

# Insert an outlier
#dat[c(1,length(dat))] <- c(-30, 100)
dat[sample(1:n, 15)] <- -20:-6
dat[sample(1:n, 15)] <- 6:20

plot(density(dat))
residuals <- dat

#Generate data
#residuals <- gen_residuals(n,mu$mu,tau1,tau2)
gen_data <- function(X,beta,residuals){
  y <- X%*%beta+ residuals
  return(y)
}

X <- round(matrix(rnorm(n*k), nrow = n, ncol = k),2) #k
X <- cbind(1, X)#cbind(1,X)
Y <- gen_data(X,beta_true,residuals)
fit <- lm(Y~.,data = data.frame(X[,-1]))
## Initialization
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

b <- coef(lm(Y~X[,-1]))

epsilon = Y-X%*%b # ols residuals

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
complete_matrix = matrix(1, nrow=n, ncol=3)

initializer_ <- bnpResiduals:::initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)

## Algorithm
begin = Sys.time()
print(begin)
trial <- bnpResiduals:::sampler(X, Y, initializer_, iter=10000, burn_in = 5000, thin = 10, log_rate = log_s )
end = Sys.time()
end - begin

## Analysis
predictive_estimates = apply(trial$predictive, 1, mean)


plot(sort(predictive_estimates), type="l", main="Predictive of the residuals")

eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui Gaussiani\\Results\\
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui Mistura di T\\Results\\
#C:\\Users\\valen\\Desktop\\Simulazioni\\Residui T1\\Results\\
m <- read.table("betas.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")


# Autocorrelation of the estimates of a coefficient
to_study = 1
beta_est = m[,to_study]
acf(beta_est)


plot(beta_est, type="l")
abline(h=beta_true[to_study], col=2)
abline(h=mean(beta_est, na.rm=T), col=3, lty=2)

# ALL the coefficients

d <- data.frame(cbind(beta_true[to_study], mean(beta_est)))
colnames(d) <- c("True Value", "Estimate")
rownames(d) <- paste("Beta", to_study)
betas_estimates <- apply(m, 2, mean)
d <- cbind(betas_estimates, beta_true, coef(lm(Y~X[,-1])))
colnames(d) = c("estimates", "true values", "OLS")
rownames(d) = NULL
library(xtable)
print(xtable(round(d,3)))


# Empirical Credible Intervals
cred_int = credible_intervals(m)
mean(cred_int[,2]-cred_int[,1])


# Tests
library(xtable)
print(xtable(data.frame(t(2*(trial$single_test_avg[,1]-trial$single_test_avg[,2])))))

# test all betas
2*(logSumExp(trial$global_test[,1])-logSumExp(trial$global_test[,2]))

# Comparing OLS and Our Proposal
x11()
ols_model <- lm(Y~X[,-1])
par(mfrow=c(1,1))
plot(beta_true, type="l", main="OLS coefficients vs Proposed coefficients", ylab="coefficients",
     ylim=c(-3,3))
polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(cred_int[ ,2]), cred_int[ ,1]), col = rgb(1,0,0, alpha=0.2), border = NA)
lines(1:(k+1), cred_int[ ,2], lty = 'dashed', col = 'red')
lines(1:(k+1), cred_int[ ,1], lty = 'dashed', col = 'red')
lines(betas_estimates, col=2)
points(betas_estimates, col=2, pch=19)
polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(confint(ols_model)[ ,2]), confint(ols_model)[ ,1]), col = rgb(0,1,0, alpha=0.1), border = NA)
lines(1:(k+1), confint(ols_model)[ ,2], lty = 'dashed', col = 'green')
lines(1:(k+1), confint(ols_model)[ ,1], lty = 'dashed', col = 'green')
lines(coef(ols_model), col=3)
points(coef(ols_model), col=3, pch=19)


legend("bottomright", y = NULL, legend=c("OLS", "Proposed"), fill = NULL, col = c(3,2),
       pch=19, cex=2, border=NULL, bty="n")

# Empirical Density
x11()
plot(density(dat), type="l", lwd=2,
     main = "Empirical density of the residuals", col="grey",  xlab="")


#lines(density(as.numeric(eps[1,])), type="l", col=2, lty=3, lwd=2)
lines(density(as.numeric(eps[dim(eps)[1],])), type="l", col=1, lty=4, lwd=2)
lines(density(fit$residuals), type="l", col=1, lty=3, lwd=1)

legend("topright", legend = c("True Density", "Proposed", "OLS"), col=c("grey",1,1),
       lty=c(1,4,3), cex=1.5, bty="n")


# Cluster Analysis
freq_clusters <- apply(trial$clusters, 2,table) # list
label <- c()

for (i in 1:length(freq_clusters)){
  #print(as.numeric(names(which(freq_clusters[[i]]==max(freq_clusters[[i]])))))
  label[i] = min(as.numeric(names(which(freq_clusters[[i]]==max(freq_clusters[[i]])))))
}

plot(as.numeric(eps[nrow(eps),]), col=ifelse(label==label[1] | label==label[length(label)], 2,1) , pch=19, ylab="Residual",
     xlab="", main="Clustered residuals", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,1))

length(unique(label)) # number of clusters

means_per_obs <- apply(mu, 2, mean)
plot(means_per_obs, pch=19, col=ifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab=expression(mu), main = "Mean", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))

tau1_per_obs <- apply(tau1, 2, mean)
plot(tau1_per_obs, pch=19, col=ifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab="", main = "First Variance", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))

tau2_per_obs <- apply(tau2, 2, mean)
plot(tau2_per_obs, pch=19, col=ifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab="", main = "Second Variance",  cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))



## If we are in two dimensions, we can plot the model and immediately see its fit
beta_true = coef(ols_model)
plot(X[,2], Y, pch=19, cex=.8, xlab="X", main="Application on Real Data")
abline(a=betas_estimates[1], b=betas_estimates[2], col=2)
abline(a=beta_true[1], b=beta_true[2], col=3)
legend("bottomright", y = NULL, legend=c("OLS", "Proposed"), fill = NULL, col = c(3,2),
       pch=19, cex=0.6, border=NULL)


## Coverage
coverage(beta_true, cred_int)
coverage(coef(ols_model), confint(ols_model))

# Mean distance from true values and mean width of the intervals
mean(confint(ols_model)[,2]-confint(ols_model)[,1])
mean(dist(betas_estimates[-1]-beta_true[-1]))

mean(dist(coef(ols_model)[-1]-beta_true[-1]))
plot(trial$sigma0, type="l")

