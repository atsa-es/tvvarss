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
## simulate & plot states
lfc <- simTVVAR(B0_lfc,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
matplot(t(lfc$states),type="l")

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
#' ## CASE 2: 1 consumer & n-1 producers
#' B0_cp <- matrix(list("cf"),nn,nn)
#' B0_cp[1:(nn-1),nn] <- "td"
#' B0_cp[nn,1:(nn-1)] <- "bu"
#' diag(B0_cp) <- 0
#' ## inspect B0
#' B0_cp
#' ## simulate & plot states
#' cp <- simTVVAR(B0_cp,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#' matplot(t(cp$states),type="l")
#'
#' ## simulate a second realization of CASE 2 using same B
#' cp2 <- simTVVAR(cp$B_mat,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#' matplot(t(cp2$states),type="l")

