## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE, warning=FALSE, message=FALSE, results="hide"-----------
#  library(devtools)
#  devtools:::install_github("eric-ward/tvvarss")

## ----install, warning=FALSE, message=FALSE, results="hide"---------------
library(rstan)
library(tvvarss)
# for optimizing stan on your machine,
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## ------------------------------------------------------------------------
library(MARSS)
data(ivesDataByWeek)
ivesDataByWeek = as.data.frame(ivesDataByWeek)
# filter months 17 - 38
ivesDataByWeek = ivesDataByWeek[which(ivesDataByWeek$`Year week` %in% seq(17,37)),]
# data from west long lake (low planktivory)
spp = c("Large Phyto", "Small Phyto", "Daphnia", "Non-daphnia")
dat = ivesDataByWeek[,spp]

# log transform
dat = as.matrix(log(dat))

## ------------------------------------------------------------------------
vecY = c(dat)
set.seed(100)

test_ind = matrix(0, nrow(dat), ncol(dat))
# sample observations for each spp separately
for(i in 1:ncol(dat)) {
  pos = as.numeric(which(!is.na(dat[,i])))
  test_ind[sample(pos, size=round(0.1*length(pos)), replace=F),i] = 1
}
test_ind = seq(1,ncol(dat)*nrow(dat))[which(test_ind==1)]

training_data = vecY
training_data[test_ind] = NA
training_data = matrix(training_data, ncol=ncol(dat))

test_data = vecY
test_data[-test_ind] = NA
test_data = matrix(test_data, ncol=ncol(dat))

## ---- eval = TRUE, echo=FALSE--------------------------------------------
## Fitting

B = matrix("zero",ncol(dat),ncol(dat))
diag(B) = "dd"
B[1,2] = "td"
B[2,c(3,4)] = "td"
B[4,2] = "bu"

mod = tvvarss(y=dat, topo=B, dynamicB=FALSE, mcmc_chain=2)

## ----eval=FALSE, echo=FALSE----------------------------------------------
#  if(file.exists("../vignettes/ives_constantB.rds")==FALSE) {
#  stanmod = tvvarss(y = training_data, B = B, include_trend = FALSE, de_mean = TRUE, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 3000, mcmc_warmup = 2000, mcmc_thin = 1, mcmc_chain = 3, dynamicB=FALSE)
#  saveRDS(stanmod, "../vignettes/ives_constantB.rds")
#  }
#
#  if(file.exists("../vignettes/ives.rds")==FALSE) {
#  stanmod = tvvarss(y = training_data, B = B, include_trend = FALSE, de_mean = TRUE, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 3000, mcmc_warmup = 2000, mcmc_thin = 1, mcmc_chain = 3)
#  saveRDS(stanmod, "../vignettes/ives.rds")
#  }

## ---- echo=FALSE, fig.pos="placeHere", fig.cap="Estimated B matrix of Ives et al. 2003. Black lines represent mean estimates, and 95% CIs; purple dashed lines represent the mean estimates from a model with static B."----
stanmod = readRDS(file="../vignettes/ives.rds")
stanmod_constant = readRDS(file="../vignettes/ives_constantB.rds")
Best = apply(extract(stanmod, c("B"))$B, c(2, 3, 4), mean)
Blow = apply(extract(stanmod, c("B"))$B, c(2, 3, 4), quantile,0.025)
Bhigh = apply(extract(stanmod, c("B"))$B, c(2, 3, 4), quantile, 0.975)

Best_constant = apply(extract(stanmod_constant, c("B"))$B, c(2, 3, 4), mean)

par(mfrow = c(ncol(dat),ncol(dat)), mgp=c(2,1,0), mai=c(0.1,0.1,0.1,0.1))
for(i in 1:ncol(dat)) {
  for(j in 1:ncol(dat)) {
  plot(Best[,i,j], type="l", ylim=range(c(Blow[,i,j], Bhigh[,i,j])), lwd=3)
    lines(Blow[,i,j])
    lines(Bhigh[,i,j])
    lines(Best_constant[,i,j], col="purple", lwd=2, lty=3)
  }
}

## ---- echo=FALSE, fig.pos="placeHere", fig.cap = "Model fitted values (line) and training (red) and test (blue) data sets. Grey lines represent the estimates from a model with static B."----
pred = apply(extract(stanmod, c("pred"))$pred, c(3,4), mean)
predlow = apply(extract(stanmod, c("pred"))$pred, c(3,4), quantile, 0.025)
predhi = apply(extract(stanmod, c("pred"))$pred, c(3,4), quantile, 0.975)

pred_constant = apply(extract(stanmod_constant, c("pred"))$pred, c(3,4), mean)
predlow_constant = apply(extract(stanmod_constant, c("pred"))$pred, c(3,4), quantile, 0.025)
predhi_constant = apply(extract(stanmod_constant, c("pred"))$pred, c(3,4), quantile, 0.975)

par(mfrow = c(2,2), mgp=c(2,1,0), mai=c(0.3,0.3,0.2,0.1))
for(i in 1:4) {
  plot(pred[,i], type="l", ylim=range(c(dat[,i], pred[,i], predlow[,i], predhi[,i]), na.rm=T), lwd=2, main = spp[i])
  lines(pred_constant[,i], col="grey", lwd=2)
  lines(predlow_constant[,i], col="grey", lwd=1)
  lines(predhi_constant[,i], col="grey", lwd=1)
  lines(predlow[,i])
  lines(predhi[,i])

  points(training_data[,i], col="red")
  points(test_data[,i], col="blue")
}

## ---- echo=FALSE, fig.pos="placeHere", fig.cap="Observed (training data) vs. Predicted. Red points represent predictions for dynamic B model; purple points represent predictions for static B model."----
par(mfrow = c(2,2), mgp=c(2,1,0), mai=c(0.5,0.5,0.3,0.1))
for(i in 1:4) {
  plot(training_data[,i], pred[,i], type="p", ylim=range(c(predlow[,i], pred[,i], predhi[,i]), na.rm=T), ylab="Predicted", xlab="Observed",col="white", main=spp[i])

  for(j in 1:nrow(training_data)) {
  lines(c(training_data[j,i], training_data[j,i]), c(predlow[j,i], predhi[j,i]), col="grey")
  }

  points(training_data[,i], pred[,i], col="red")
  points(training_data[,i], pred_constant[,i], col="purple")
  abline(0,1)
}

## ---- echo=FALSE, fig.pos="placeHere", fig.cap="Predicted vs Observed (test data). Purple points represent estimates for static B model."----

par(mfrow = c(2,2), mgp=c(2,1,0), mai=c(0.5,0.5,0.3,0.1))
for(i in 1:4) {
  plot(test_data[,i], pred[,i], type="p", ylim=range(c(predlow[,i], pred[,i], predhi[,i]), na.rm=T), ylab="Predicted", xlab="Observed",col="white", main=spp[i])
  for(j in 1:nrow(training_data)) {
  lines(c(test_data[j,i], test_data[j,i]), c(predlow[j,i], predhi[j,i]), col="grey")
  }
  points(test_data[,i], pred[,i], col="blue")
  points(test_data[,i], pred_constant[,i], col="purple")
  abline(0,1)
}

## ---- echo=FALSE, fig.pos="placeHere", fig.cap="Comparison of CIs (test data) from static vs dynamic B model"----

par(mfrow = c(2,2), mgp=c(2,1,0), mai=c(0.5,0.5,0.3,0.1))
for(i in 1:4) {
  plot(predhi[,i]-predlow[,i], predhi_constant[,i]-predlow_constant[,i], type="p", xlab="Dynamic B", ylab="Static B", main=spp[i])

  abline(0,1)
}

