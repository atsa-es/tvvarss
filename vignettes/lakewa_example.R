## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  devtools::install_github("eric-ward/tvvarss")

## ----install, results="hide"---------------------------------------------
library(knitr)
library(rstan)
library(tvvarss)
# for optimizing stan on your machine,
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## ------------------------------------------------------------------------
library(MARSS)
data(lakeWAplankton)
lakeWAplankton = lakeWAplanktonRaw[,c("Diatoms", "Greens", "Bluegreens", "Daphnia", "Cyclops", "Diaptomus", "Non.daphnid.cladocerans")]
head(lakeWAplankton)

## ---- echo=FALSE---------------------------------------------------------
m = matrix("", 7, 1)
m[,1] = c("Diatoms","Greens","Bluegreens","Cyclops","Daphnia","Diaptomus","Non-daphnid cladocerans")
colnames(m) = "Species"
kable(m)

## ------------------------------------------------------------------------
for(i in 1:ncol(lakeWAplankton)) {
  zeros = which(lakeWAplankton[,i]==0)
  lakeWAplankton[zeros,i] = runif(n = length(zeros), min = 0, max = 0.5 * min(lakeWAplankton[,i], na.rm=T))
}

# log transform
lakeWAplankton_log = log(lakeWAplankton)

## ----bmatrix-------------------------------------------------------------
B = matrix("zero", ncol(lakeWAplankton_log), ncol(lakeWAplankton_log))
diag(B) = "dd"
B[1,3] = "cf"
B[4,c(5,7)] = "cf"
B[5,c(1,3)] = "cf"
B[6,c(1,5,7)] = "cf"
B[7,c(1,5)] = "cf"

## ---- eval = FALSE-------------------------------------------------------
#  B = B0_lfc
#  
#  stanmod = tvvarss(y = dat_obs, include_trend = FALSE, de_mean = TRUE, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 1000, mcmc_warmup = 500, mcmc_thin = 1, mcmc_chain = 3)

## ---- eval = FALSE-------------------------------------------------------
#  pred = apply(extract(stanmod, c("pred"))$pred, c(3,4), mean)
#  
#  par(mfrow = c(2,2), mgp=c(2,1,0), mai=c(0.3,0.3,0.1,0.1))
#  for(i in 1:4) {
#    plot(pred[,i], type="l", ylim=range(c(dat_obs[,i], pred[,i])))
#    points(dat_obs[,i], col="red")
#  }

