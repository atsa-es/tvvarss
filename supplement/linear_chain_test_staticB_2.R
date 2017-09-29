if(!require("tvvarss")) {
  devtools::install_github("nwfsc-timeseries/tvvarss")
  library("tvvarss")
}
if(!require("broom")) {
  install.packages("broom")
  library("broom")
}
if(!require("rstan")) {
  install.packages("rstan")
  library("rstan")
}

library(broom)
library(rstan)
library(tvvarss)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 4
n_year <- 60
n_site <- 1
n_simulations = 100
## min log-density threshold
dens_min <- -3
## max log-density threshold
dens_max <- 3

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_species,n_species)
diag(B0_lfc) <- "dd"
for(i in 1:(n_species-1)) {
  B0_lfc[i,i+1] <- "td"
  B0_lfc[i+1,i] <- "bu"
}

## initial B
B0_init <- matrix(0,n_species,n_species)
diag(B0_init) <- c(0.9480746 0.7842900 0.7794805 0.8233147)
B0_init[1,2] <- -0.05
B0_init[2,1] <- 0.15
B0_init[2,3] <- -0.1
B0_init[3,2] <- 0.15
B0_init[3,4] <- -0.3
B0_init[4,3] <- 0.15

## empty list for results
saved_output = vector("list", n_simulations)

## function for fitting/extracting
get_coefs <- function(y, topo, shared_r, shared_q, mcmc_chain, mcmc_iter, mcmc_warmup) {
  fitted_model <- tvvarss(y=y, topo=topo, shared_r=shared_r, shared_q=shared_q,
    mcmc_chain=mcmc_chain, mcmc_iter=mcmc_iter, mcmc_warmup=mcmc_warmup, dynamicB = FALSE)
  #coef <- tidy(fitted_model, intervals=TRUE, prob=0.9)
  return(fitted_model)
}

for(ns in 1:n_simulations) {
  ## simulate process
  lfc <- simTVVAR(Bt = B0_init,
                  topo = B0_lfc,
                  TT = n_year,
                  var_QX = 0.02,
                  cov_QX = 0,
                  var_QB = 0.01,
                  cov_QB = 0)
  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    lfc <- simTVVAR(Bt = B0_init,
                    topo = B0_lfc,
                    TT = n_year,
                    var_QX = 0.02,
                    cov_QX = 0,
                    var_QB = 0.01,
                    cov_QB = 0)
  }
  ## add obs error
  Y <- sim2fit(lfc, n_site, sd=0.1)
  ## fit model to only 2nd half of data & save param summaries

  y_thin = array(0, dim=c(1, n_year/2, 4))
  y_thin[1,,] = Y[,-c(1:(n_year/2)),]
  fitted <- get_coefs(y = y_thin, topo = B0_lfc,
                    shared_r = matrix(1,4,1),
                    shared_q = matrix(1,4,1),
                    mcmc_chain = 3, mcmc_iter = 3000, mcmc_warmup = 2000)
  coef = tidy(fitted, intervals=TRUE, prob=0.9)
  ## save data (Y), simulation output (lfc), model coefficients
  saved_output[[ns]] <- list('data' = Y, 'sim_output' = lfc, 'estimate' = coef)
}
## save results
save(saved_output, file = "output_linear_chain_staticB2.Rdata")


# Plot output
# B0_init contains true values

par(mfrow=c(4,4), mgp=c(2,1,0), mai=c(0.1,0.01,0.1,0.01))
for(i in 1:4) {
  for(j in 1:4) {
    idx = which(startsWith(paste0("B[",1,",",i,",",j,"]"), saved_output[[1]]$estimate$term))
    z = 0
    for(k in 1:100) {
      z[k] = saved_output[[k]]$estimate$estimate[idx]
    }
    if(var(z) == 0) {
      plot(0,0,col="white",axes=F,main=paste0("B[",i,",",j,"]"))
    } else {
    hist(z, col="grey70", axes=F, main=paste0("B[",i,",",j,"]"))
    lines(rep(B0_init[i,j], 2), c(0, 10000), lwd=3, col="red")
    }
  }
}


## NOT RUN
## plot states
## lvls <- c("TC","SC","PC","PP")
## clr <- c("purple","darkred","blue","darkgreen")
## par(mai=c(0.9,0.9,0,0), omi=c(0.1,0.1,0.1,1.5))
## matplot(Y[,,], type="l", lty="solid", lwd=2, xpd=NA,
##         col=clr, ylab="Log density")
## legend("right", legend=lvls, lty="solid", lwd=2, bty="n",
##        col=clr, inset=-0.2, xpd=NA, cex=0.9)
