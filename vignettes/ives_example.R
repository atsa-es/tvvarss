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

