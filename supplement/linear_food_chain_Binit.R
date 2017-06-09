devtools::install_github("nwfsc-timeseries/tvvarss")
library(tvvarss)
library(tseries)
library(broom)
library(rstan)

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
diag(B0_init) <- 0.6
B0_init[1,2] <- -0.03
B0_init[2,1] <- 0.1
B0_init[2,3] <- -0.1
B0_init[3,2] <- 0.1
B0_init[3,4] <- -0.3
B0_init[4,3] <- 0.1

## empty list for results
saved_output = vector("list", n_simulations)

for(ns in 1:n_simulations) {
  ## simulate process
  lfc <- list(states=2*rep(dens_max,2))
  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    lfc <- simTVVAR(Bt = B0_init,
                    topo = B0_lfc,
                    TT = n_year,
                    var_QX = rev(seq(1, 4) / 40),
                    cov_QX = 0,
                    var_QB = 0.01,
                    cov_QB = 0)
  }
  ## add obs error
  dat <- sim2fit(lfc, n_site, sd=0.1)
  ## drop the burn in (1st half of the states/B matrix)
  Bmat = lfc$B_mat[,,-c(1:(n_year/2 + 1))]
  Y = array(0, dim = c(dim(dat)[1], dim(dat)[2]/2, dim(dat)[3]))
  for(j in 1:n_site) {
    Y[j,,] = dat[j,-c(1:(n_year/2)),]
  }
  ## data has been generated -- now fit the model using tvvarss()
  fitted_model = tvvarss(y = Y, topo = B0_lfc, shared_r = matrix(1, n_species, n_site),
                         mcmc_iter = 3000, mcmc_warmup = 2000)
  ## tidy the stanfit object to save space -- don't save every mcmc draw
  coef = tidy(fitted_model, intervals = TRUE, prob = 0.9)
  ## save data (Y), simulation output (lfc), model coefficients
  saved_output[[ns]] = list('data' = Y, 'sim_output' = lfc, 'estimate' = coef)
}
## save results
save(saved_output, file = "output_linear_chain.Rdata")


## NOT RUN
## plot states
## lvls <- c("TC","SC","PC","PP")
## clr <- c("purple","darkred","blue","darkgreen")
## par(mai=c(0.9,0.9,0,0), omi=c(0.1,0.1,0.1,1.5))
## matplot(Y[,,], type="l", lty="solid", lwd=2, xpd=NA,
##         col=clr, ylab="Log density")
## legend("right", legend=lvls, lty="solid", lwd=2, bty="n",
##        col=clr, inset=-0.2, xpd=NA, cex=0.9)
