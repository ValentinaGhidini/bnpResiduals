setwd("C:\\Users\\valen\\Desktop\\Misture di Misture")
source("code complete/bnp_residuals_blockedGS.R")

# C:\Users\valen\Desktop\Misture di Misture\data
# Data
# Mfunds

library(zoo)
load("data/FamaFrench.RData")
load("data/mfunds.RData")
load("data/SP500.RData")
load("data/GAFA.RData")
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


b = coef(lm(Y~X))
beta_true = coef(lm(Y~X))
