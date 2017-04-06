## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE-------------------------------------------------
#  library(rstan)
#  library(devtools)
#  devtools::install_github("eric-ward/tvvarss")
#  library(tvvarss)
#  # for optimizing stan on your machine,
#  rstan_options(auto_write = TRUE)
#  options(mc.cores = parallel::detectCores())

## ------------------------------------------------------------------------
set.seed(123)
#' ## number of time steps
TT <- 30
#' ## number of spp/guilds
nn <- 4
## CASE 1: linear food chain
B0_lfc <- matrix(list(0),nn,nn)
for(i in 1:(nn-1)) {
  B0_lfc[i,i+1] <- "td"
  B0_lfc[i+1,i] <- "bu"
}
for(i in 1:nn) {
  B0_lfc[i,i] = "bu"
}

## simulate & plot states
lfc <- simTVVAR(B0_lfc,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
matplot(lfc$states,type="l")

## ------------------------------------------------------------------------
dat_obs = lfc$states + matrix(rnorm(dim(lfc$states)[1]*dim(lfc$states)[2], 0, 0.1), nrow=dim(lfc$states)[1])

