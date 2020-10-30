setwd("C:\\Users\\valen\\Desktop\\Misture di Misture")
source("code complete/bnp_residuals_blockedGS.R")

# C:\Users\valen\Desktop\Misture di Misture\data
# Data
# Mfunds

library(zoo)
load("data/FamaFrench.RData")
load("data/mfunds.RData")
load("data/GAFA.RData")
load("data/SP500.RData")

GSATX <- mfunds$GSATX

data <- mfunds$GSATX
data <- as.vector(data[4153:length(data)])
GSATX <- rep(0,length(data)-1)
for(i in 2:length(data)){
  GSATX[i-1] <- 100*(data[i]-data[i-1])/data[i-1] #compute returns
}

dat <- window(FamaFrench, start = "2007-06-25") #start of GSATX
dat <- as.data.frame(dat)
dat <- data.frame(dat, GSATX = GSATX[1:2775]) #end of FamaFrench

Y <- dat$GSATX - dat$RF
X <- cbind(1,dat$Mkt.RF,dat$SMB,dat$HML)


### DATA - THE SECOND DATASET
R <- 100 * diff(GAFA) / head(GAFA, -1)

# merge Fama-French data with GAFA's returns
# keep required time window
dat <- merge(FamaFrench, R)
dat <- window(dat, start = "2016-04-30", end = "2018-05-31") 
dat <- as.data.frame(dat)

# create variables for CAPM model for Google's stock 
Y <- dat$GOOGLE - dat$RF # E(Ri) - Rf
X <- dat$Mkt.RF # E(Rm) - Rf
X <- cbind(1, X) # add column of 1s for intercept



b = coef(lm(Y~X))
beta_true = coef(lm(Y~X))


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
Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
complete_matrix = matrix(1, nrow=n, ncol=3)


initializer <- initializer(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix)

## Algorithm
begin = Sys.time()
print(begin)
trial <- sampler(X, Y, initializer, iter=5000, burn_in = 100, thin = 10, log_rate = log_s )
end = Sys.time()
end - begin 


## Analysis
predictive_estimates = apply(trial, 1, mean)
plot(predictive_estimates, type="l", main="Predictive of the residuals")

eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
m <- read.table("betas.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")


# Autocorrelation of the estimates of a coefficient
to_study = 2
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
print(xtable(round(d,3)))


# Empirical Credible Intervals
cred_int = credible_intervals(m)
mean(cred_int[,2]-cred_int[,1])

# Comparing OLS and Our Proposal

ols_model <- lm(Y~X[,-1])

plot(beta_true, type="l", main="OLS coefficients vs Proposed coefficients")
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
       pch=19, cex=0.6, border=NULL)


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