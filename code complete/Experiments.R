# Experiments

# Import bnpResiduals - be sure to have the latest version installed!
library(bnpResiduals)
library(reshape2)
library(ggplot2)
library(xtable)


# Check the documentation - suggestions are welcome!
?bnp.lm
?bnp.ar
?bnp.GP


########### Generate synthetic data ###################

# Some functions to generate weird residuals
ugly_residuals <- function(n){
  mu <- 0
  res <- c()
  s = 1
  for(i in 1:n){

    res[i] <- rnorm(1, s*(mu+i), sd= sqrt(i))
    s = s*(-1)
  }
  return(res)
}

gen_residuals_fourgaussians <- function(n, mu1, mu2, tau1, tau2, tau3, tau4){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(T,F),n, replace = T)
  print(bol)
  res[bol] <- gen_residuals(sum(bol), mu1, tau1, tau2)
  res[!bol] <- gen_residuals(n-sum(bol),mu2, tau3, tau4)
  return(res)
}

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

gen_residuals_trimodal <- function(n, mu, tau1, tau2, tau3){
  #mu: mean (positive value)
  #tau1: variance of positive part
  #tau2: variance of negative part
  res <- rep(0,n)
  bol <- sample(c(1,2,3),n, replace = T)
  res[bol==1] <- rnorm(sum(bol==1), mean=0, sd = tau1)
  res[bol==2] <- rnorm(sum(bol==2), mean=-mu, sd = tau2)
  res[bol==3] <- rnorm(sum(bol==3), mean=mu, sd = tau3)
  return(res)
}

# Generate data

n <- 500 # number of observations
k <- 9 # number of covariates
beta_true = beta = c(rep(-1, 3), rep(0,1), rep(1,6)) # real coefficients of the linear model to simulate from
intercept = T
if (!intercept){ beta_true = beta = beta[-1]}

# Generate residuals

# Set some hyperparameters
muu = 100
tau1 = 1
tau2 = 50


residuals <- ugly_residuals(n)
#dat <- gen_residuals_fourgaussians(n, 1000,50,20,20,5,5)
#dat <- gen_residuals_trimodal(n, 200, 10, 20, 10)
#dat <- rnorm(n)
#dat <- runif(n, -50, 50)
#dat <- gen_residuals(n, muu, tau1, tau2)
#dat <- rt(n, df=1)
#dat <- ifelse(sample(size=n, x=c(0,1), replace=T)==1, gen_residuals(n, 30, 2, 10), gen_residuals(n, 200, 150, 50))
#dat <- gen_residuals_T(n, 2, 3, 2)
# Density of the residuals

# Insert outliers in the residuals outlier
#dat[c(1,length(dat))] <- c(-30, 100)
#dat[sample(1:n, 15)] <- -20:-6
#dat[sample(1:n, 15)] <- 6:20

plot(density(residuals), main = "Density of the synthetic residuals")


# Generate data
genLinear_data <- function(X,beta,residuals){
  y <- X%*%beta+ residuals
  return(y)
}

# Covariates
X <- round(matrix(rnorm(n*k), nrow = n, ncol = k),2) #k
if (intercept){ X <- cbind(1, X)}

# Linear data
Y <- genLinear_data(X,beta_true,residuals)
fit <- lm(Y~.,data = data.frame(X[,-1])) # ols

# Non-linear data

X <- matrix(seq(-100, 100), by=.1, ncol=1)
n = length(X)
Y = sin(X) + rnorm(n, 0, 0.2)
x.star <- seq(-20, 20, by = 0.5) # predictive grid

# Initialization

n <- nrow(X) # just to be sure
k <- ncol(X)

# For linear models, initialize the coefficients of the prior
b0 <-  rep(0, k)
B0 <-  diag(10, k)

s1 <-  2
s2 <-  2
S1 <- 1#0.01
S2 <-  1#30

sigma2=10
theta = 1

b = coef(fit)

if (intercept){
  epsilon = Y-X%*%b # ols residuals
}
if (!intercept) {
  epsilon = rnorm(n) # ols residuals
}

comp <- ifelse(epsilon > 0, 1, -1)

N = 200 # Number of elements for the blocked Gibbs sampler
log_s = T # compute the probabilities in logscale or not

# Initialiaze blocked gibbs sampler
V <- rbeta(N-1,shape1 = 1,theta)
W <- cumprod(1-V)
V <- c(V,1)
W <- c(1,W)
p <- V*W

mu0 = 1
sigma0 = sqrt(25)

Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
complete_matrix = matrix(1, nrow=n, ncol=3)

# Create the initialize object
initializer_ <- bnpResiduals:::initializer(Z, p, comp,  N, mu0, sigma0, complete_matrix, b0, B0)


############################################## Gaussian Process ##############################################################

ss <- bnp.GP(X, Y, initializer_, iter = 1000, Xstar=x.star)

values <- cbind(x = x.star, ss$preds)
prediction = apply(ss$preds, 1, mean)

values <- melt(values,id="x")

df <- data.frame(cbind(x.star = ss$Xstar[,1], f.star.bar=prediction))
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
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

plot(mu[,10], type="l", main="Trace plot", ylab="mu")
plot(tau1[,40], type="l", main="Trace plot", ylab="tau1")
plot(tau2[,40], type="l", main="Trace plot", ylab="tau2")


################################################## Linear model #####################################################

model <- bnp.lm(X, Y, initializer_, cluster = T, iter=5000, burn_in = 500, thin=10, sigma=0.8)
# sigma = .5 non cambia un granch? - ma la densit? stimata ? pi? smooth
# sigma = .8 bene! anche densit?
# sigma = 1 bene! anche densit?, ma un solo cluster ((???))

# ugly residuals
# sigma=0 --- SCHIFO - forse solo la shape della densit? si salva
# sigma=0.5 --- Non bene - forse troppa  variabilit??


## Analysis
trial <- bnpResiduals:::sampler(X, Y, initializer_, iter=1000, burn_in = 100, thin = 10, log_rate = log_s )
predictive_estimates = apply(trial$predictive, 1, mean)
plot(sort(predictive_estimates), type="l", main="Predictive of the residuals")

eps <- read.table("epsilon.csv", header=F, skip=1, sep=";")
m <- read.table("betas.csv", header=F, skip=1, sep=";")
tau1 <-  read.table("tau1.csv", header=F, skip=1, sep=";")
tau2 <-  read.table("tau2.csv", header=F, skip=1, sep=";")
mu <-  read.table("mu.csv", header=F, skip=1, sep=";")

# Real vs empirical densities
par(mfrow=c(1,1))
plot(density(residuals))
lines(density(as.numeric(eps[dim(eps)[1],])), type="l", col=2, lty=4, lwd=2, main = "Density of the Empirical Residuals")

# Check the traceplot, estimate and the autocorrelation of the to_study coefficient
to_study = 2
beta_est = m[,to_study]
acf(beta_est)
plot(beta_est, type="l")
abline(h=beta_true[to_study], col=2)
abline(h=mean(beta_est, na.rm=T), col=3, lty=2)


# Let's see the estimates for all the coefficients
d <- data.frame(cbind(beta_true[to_study], mean(beta_est)))
colnames(d) <- c("True Value", "Estimate")
rownames(d) <- paste("Beta", to_study)

betas_estimates <- apply(m, 2, mean)
d <- cbind(betas_estimates, beta_true, coef(lm(Y~X[,-1])))
colnames(d) = c("estimates", "true values", "OLS")
rownames(d) = NULL
print(xtable(round(d,3)))

# Empirical Credible Intervals
cred_int = bnpResiduals:::credible_intervals(m)
mean(cred_int[,2]-cred_int[,1])
print(xtable(data.frame(t(2*(trial$single_test_avg[,1]-trial$single_test_avg[,2])))))

# Test of fit
2*(logSumExp(trial$global_test[,1])-logSumExp(trial$global_test[,2]))


# Comparing OLS and Our Proposal
x11()
ols_model <- lm(Y~X[,-1])
par(mfrow=c(1,1))
plot(beta_true, type="l", main="OLS coefficients vs Proposed coefficients", ylab="coefficients",
     ylim=c(-3,3))
polygon(c(rev(1:k), 1:k), c(rev(cred_int[ ,2]), cred_int[ ,1]), col = rgb(1,0,0, alpha=0.2), border = NA)
lines(1:k, cred_int[ ,2], lty = 'dashed', col = 'red')
lines(1:k, cred_int[ ,1], lty = 'dashed', col = 'red')
lines(betas_estimates, col=2)
points(betas_estimates, col=2, pch=19)
polygon(c(rev(1:k), 1:k), c(rev(confint(ols_model)[ ,2]), confint(ols_model)[ ,1]), col = rgb(0,1,0, alpha=0.1), border = NA)
lines(1:k, confint(ols_model)[ ,2], lty = 'dashed', col = 'green')
lines(1:k, confint(ols_model)[ ,1], lty = 'dashed', col = 'green')
lines(coef(ols_model), col=3)
points(coef(ols_model), col=3, pch=19)
legend("bottomright", y = NULL, legend=c("OLS", "Proposed"), fill = NULL, col = c(3,2),
       pch=19, cex=2, border=NULL, bty="n")

# Empirical Density
x11()
plot(density(residuals), type="l", lwd=2,
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

plot(as.numeric(eps[nrow(eps),]), col=label,#ifelse(label==label[1] | label==label[length(label)], 2,1) , pch=19, ylab="Residual",
     xlab="", main="Clustered residuals", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,1), pch=19, ylab="Residual")

label[c(1, length(label))] = label[c(1, length(label))]  + 1
means_per_obs <- apply(mu, 2, mean)
plot(means_per_obs, pch=19, col=label,#ifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab=expression(mu), main = "Mean", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))

tau1_per_obs <- apply(tau1, 2, mean)
plot(tau1_per_obs, pch=19, col=label, #ifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab="", main = "First Variance", cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))

tau2_per_obs <- apply(tau2, 2, mean)
plot(tau2_per_obs, pch=19, col=label,#lifelse(label==label[1] | label==label[length(label)], 2, 1),
     xlab="Observation", ylab="", main = "Second Variance",  cex=ifelse(label==label[1] | label==label[length(label)], 1.5,.6))



## If we are in two dimensions, we can plot the model and immediately see its fit
beta_true = coef(ols_model)
plot(X[,2], Y, pch=19, cex=.8, xlab="X", main="Application on Real Data")
abline(a=betas_estimates[1], b=betas_estimates[2], col=2)
abline(a=beta_true[1], b=beta_true[2], col=3)
legend("bottomright", y = NULL, legend=c("OLS", "Proposed"), fill = NULL, col = c(3,2),
       pch=19, cex=0.6, border=NULL)

## Coverage
bnpResiduals:::coverage(beta_true, cred_int)
bnpResiduals:::coverage(coef(ols_model), confint(ols_model))

# Mean distance from true values and mean width of the intervals
mean(confint(ols_model)[,2]-confint(ols_model)[,1])
mean(dist(betas_estimates[-1]-beta_true[-1]))
mean(dist(coef(ols_model)[-1]-beta_true[-1]))
plot(trial$sigma0, type="l", main= "Trace Plot")


######################################## LEFTOVERS ################################################################

# Replicating tables

initializer <- function(Z, p, comp, b0, B0, N, mu0, sigma0, complete_matrix){
  return(list(Z = Z, p = p, comp = comp, b0 = b0, B0 = B0, N=N, mu0=mu0, sigma0=sigma0, complete_matrix=complete_matrix))

}


muu <- c(10,30,100)
dff<- c(2.1, 5,10,50)
total <- matrix(nrow = 4, ncol= 6)
samplesize <- c(15,50,500,1000)

#n = 100
library(attempt)
#beta_true <- c(-1, beta_true)#beta_true[-1]
fig=F
for (s in 1:length(samplesize)){
  #print(muu[s])
  cover <- c()
  cover_ols <- c()
  mean_width <- c()
  mean_width_ols <- c()
  mean_dist <- c()
  mean_dist_ols <- c()
  test_gs <- c()
  test_ols <- c()
  #k = length(beta_true)
  n =samplesize[s]
  X <- round(matrix(rnorm(n*k), nrow = n, ncol = k),2) #k
  X <- cbind(1, X)#cbind(1,X)


  for (i in 1:20){

    print(i)

    Y <- gen_data(X, beta_true, rnorm(n))#gen_residuals_T(n, 2, 3, muu[s]))#rnorm(n))#rt(n, df=dff[s]))# rnorm(n))##gen_residuals_T(n, 2, 3, muu[s]))####gen_residuals_T(n, 2, 3, muu[s]))
    #print(Y)#)#df=2.1
    epsilon = Y-X%*%coef(lm(Y~X-1)) # ols residuals
    comp <- ifelse(epsilon > 0, 1, -1)
    ols_model <- lm(Y~X)
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

    b <- coef(lm(Y~X-1))

    N = 150

    log_s = T
    V <- rbeta(N-1,shape1 = 1,theta)
    W <- cumprod(1-V)
    V <- c(V,1)
    W <- c(1,W)
    p <- V*W

    mu0 = 0
    sigma0 = 1
    K <- sample(1:N, size = n, replace = T, prob = p)
    Z <- cbind(rtruncnorm(N,a = 0, mean = mu0, sd = sigma0), sqrt(rinvgamma(N,shape = s1, rate = S1)), sqrt(rinvgamma(N,shape = s2, rate = S2)))
    temp <- bnpResiduals:::sampler_triplet_blockedGS(Y,epsilon, comp, p, Z, mu0, sigma0)
    K <- temp$K
    Z <- temp$Z
    complete_matrix = Z[K,]

    initializer <- list(Z = Z, p = p, comp = comp, b0 = b0, B0 = B0, N=N, mu0=mu0, sigma0=sigma0, complete_matrix=complete_matrix)


    attempt(trial <- sampler(X, Y, initializer, iter=500, burn_in = 50, thin = 10, log_rate = log_s ), "error")
    m <- read.table("betas.csv", header=F, skip=1, sep=";")#[,-1]
    intervals = credible_intervals(m)#[-1,]
    #print(apply(m, 2, mean))
    betas_estimates <- apply(m, 2, mean)
    #cover[i] <- sum(beta_true>=intervals[,1] & beta_true[-1]<=intervals[,2])/length(beta_true)
    intervals_ols <- confint(lm(Y~X-1))
    #cover_ols[i] <- sum(beta_true>=intervals_ols[,1] & beta_true<=intervals_ols[,2])/length(beta_true)
    mean_width[i] <- mean(intervals[,2]-intervals[,1])
    mean_width_ols[i] <- mean(intervals_ols[,2]-intervals_ols[,1])
    mean_dist[i] <- dist(apply(m, 2, mean)-beta_true)
    mean_dist_ols[i] <- dist(coef(lm(Y~X-1))-beta_true)
    #print(apply(m, 2, mean))

    #
    #print(intervals)
    beta_true_bol <- (beta_true==0)
    intervals_bol <- (intervals[,1]<0 & intervals[,2]>0)

    test_gs[i] <- sum(beta_true_bol == intervals_bol)
    #print(intervals_bol)
    #print(beta_true_bol)
    # print(test_gs)
    #test_gs[i] <- sum((intervals[,1]<0 & intervals[,2]>0 & beta_true==0) |
    #                   (((intervals[,1]<0 & intervals[,2]<0) |(intervals[,1]>0 & intervals[,2]>0))& beta_true!=0))/length(beta_true)


    test_ols[i] <- sum((intervals_ols[,1]<0 & intervals_ols[,2]>0 & beta_true==0) |
                         (((intervals_ols[,1]<0 & intervals_ols[,2]<0) |(intervals_ols[,1]>0 & intervals_ols[,2]>0))& beta_true!=0))/length(beta_true)
    fig=F
    if (fig){

      x11()
      ols_model <- lm(Y~X[,-1])
      par(mfrow=c(1,1))
      plot(beta_true, type="l", main="OLS coefficients vs Proposed coefficients", ylab="coefficients",
           ylim=c(-3,3))
      polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(intervals[ ,2]), intervals[ ,1]), col = rgb(1,0,0, alpha=0.2), border = NA)
      lines(1:(k+1), intervals[ ,2], lty = 'dashed', col = 'red')
      lines(1:(k+1), intervals[ ,1], lty = 'dashed', col = 'red')
      lines(betas_estimates, col=2)
      points(betas_estimates, col=2, pch=19)
      polygon(c(rev(1:(k+1)), 1:(k+1)), c(rev(confint(ols_model)[ ,2]), confint(ols_model)[ ,1]), col = rgb(0,1,0, alpha=0.1), border = NA)
      lines(1:(k+1), confint(ols_model)[ ,2], lty = 'dashed', col = 'green')
      lines(1:(k+1), confint(ols_model)[ ,1], lty = 'dashed', col = 'green')
      lines(coef(ols_model), col=3)
      points(coef(ols_model), col=3, pch=19) }


  }


  total[s,] <- rbind(mean(test_gs/length(beta_true)),
                     mean(test_ols),
                     mean(mean_width),
                     mean(mean_width_ols),

                     mean(mean_dist),
                     mean(mean_dist_ols))

}


#### AUTOREGRESSIVE MODEL - CHECK ############
n <- 100 # length of the time serie
k <- 3 # let's include k covariates

# Simulate the starting data
Y <- rnorm(n)
X <- matrix(rnorm(n*k), ncol=k)



qq = 2 # AR(qq)
# Let's compute the exact data
for (t in qq:(length(Y)-1)){
  Y[t+1] =  .5*Y[t] + .5*Y[t-1] + X[t+1,1] + X[t+1,2]+rnorm(1, 0, 5)#X[t,]%*%delta +Y[(t-count):(t-count+qq-1)]%*%gamma + rnorm(1)#gen_residuals(1, 100, 30, 3) #+
}


X_design <- bnpResiduals:::designmatrix.ar(X, Y, qq)
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




