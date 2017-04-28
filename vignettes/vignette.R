## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE-------------------------------------------------
#  library(rstan)
#  library(tvvarss)
#  ## for optimizing stan on your machine,
#  rstan_options(auto_write = TRUE)
#  options(mc.cores = parallel::detectCores())

## ----data----------------------------------------------------------------
y = array(rnorm(240), dim=c(3, 20, 4))

## ----function, echo=TRUE, eval=FALSE-------------------------------------
#  tvvarss(y = y)

## ----function2, echo=TRUE, eval=FALSE------------------------------------
#  tvvarss(y = y, include_trend = TRUE, de_mean = TRUE, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 1000, mcmc_warmup = 500, mcmc_thin = 1, mcmc_chain = 3)

## ---- echo=FALSE---------------------------------------------------------
matrix(c(rep(0,4), rep(10,4), rep(5,4)), 4, 3)

## ---- echo=FALSE---------------------------------------------------------
matrix(rep(1:7,3), 7, 3)

## ---- echo=FALSE---------------------------------------------------------
matrix(seq(1,7*3), 7, 3)

## ---- eval=FALSE---------------------------------------------------------
#  tvvarss(y, shared_q = matrix(1, 4, 3), shared_r = matrix(1, 4, 3), shared_u = matrix(1, 4, 3))

## ---- eval=FALSE---------------------------------------------------------
#  plot(fitted_model, pars = c("B[3,2,2]", "u[1]"))
#  traceplot(fitted_model, pars = c("B[3,2,2]", "u[1]"))

## ---- eval=FALSE---------------------------------------------------------
#  pars = extract(fitted_model, permuted = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  library(broom)
#  tidy_pars = tidy(fitted_model)

## ---- eval=FALSE---------------------------------------------------------
#  library(loo)
#  loo(extract_log_lik(fitted_model))$looic

