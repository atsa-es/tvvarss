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

## ---- results="hide"-----------------------------------------------------
for(i in 1:ncol(lakeWAplankton)) {
  zeros = which(lakeWAplankton[,i]==0)
  if(length(zeros)>0) lakeWAplankton[zeros,i] = NA
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

## ------------------------------------------------------------------------
vecY = c(lakeWAplankton_log)
set.seed(100)
test_ind = sample(1:length(vecY), 
  size=round(0.1*length(vecY)), replace=F)

training_data = vecY
training_data[test_ind] = NA
training_data = matrix(training_data, ncol=ncol(lakeWAplankton_log))

test_data = vecY
test_data[-test_ind] = NA
test_data = matrix(test_data, ncol=ncol(lakeWAplankton_log))

## ----fitmodel------------------------------------------------------------
# If model hasn't been run
if(file.exists("vignettes/lakewa.rds")) {
stanmod = tvvarss(y = training_data, include_trend = FALSE, de_mean = TRUE, B = B, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 3000, mcmc_warmup = 2000, mcmc_thin = 1, mcmc_chain = 3)

saveRDS(stanmod, "vignettes/lakewa.rds")
}

## ---- eval = FALSE-------------------------------------------------------
#  Best = apply(extract(stanmod, c("B"))$B, c(2, 3, 4), mean)
#  
#  par(mfrow = c(7,7), mgp=c(2,1,0), mai=c(0.1,0.1,0.1,0.1))
#  for(i in 1:7) {
#    for(j in 1:7) {
#    plot(Best[,i,j], type="l")
#    }
#  }

## ---- eval = FALSE-------------------------------------------------------
#  pred = apply(extract(stanmod, c("pred"))$pred, c(3,4), mean)
#  
#  par(mfrow = c(4,2), mgp=c(2,1,0), mai=c(0.3,0.3,0.1,0.1))
#  for(i in 1:7) {
#    plot(pred[,i], type="l", ylim=range(c(lakeWAplankton_log[,i], pred[,i]), na.rm=T))
#    points(lakeWAplankton_log[,i], col="red", cex=0.1)
#  }

