##-------
## inits
##-------

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
if(!require("doParallel")) {
  install.packages("doParallel")
  library("doParallel")
}
if(!require("foreach")) {
  install.packages("foreach")
  library("foreach")
}

rstan_options(auto_write = TRUE)
n_cores <- parallel::detectCores()
registerDoParallel(n_cores-1)
## options(mc.cores = n_cores)

##----------------
## function defns
##----------------

## define function for simulating/fitting
simfit <- function(sim, fit) {
  ## extract list elements
  for(i in 1:length(sim)) {
    assign(names(sim[i]), sim[[i]], inherits = TRUE)
  }
  for(i in 1:length(fit)) {
    assign(names(fit[i]), fit[[i]], inherits = TRUE)
  }
  ## simulate process
  lfc <- simTVVAR(Bt = B0_init, topo = B0_lfc, TT = n_year,
                  var_QX = var_QX, cov_QX = cov_QX,
                  var_QB = var_QB, cov_QB = cov_QB)
  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    lfc <- simTVVAR(Bt = B0_init, topo = B0_lfc, TT = n_year,
                    var_QX = var_QX, cov_QX = cov_QX,
                    var_QB = var_QB, cov_QB = cov_QB)
  }
  ## add obs error
  Y <- sim2fit(lfc, n_site, sd=0.1, new_real=FALSE)
  ## fit the model to first half of data
  fitted_model <- tvvarss(y = Y[,-c(1:(n_year/2)),],
                          topo = topo,
                          shared_r = shared_r,
                          mcmc_chain = mcmc_chain,
                          mcmc_iter = mcmc_iter,
                          mcmc_warmup = mcmc_warmup)
  coef <- tidy(fitted_model,
               conf.int = intervals,
               conf.level = prob,
               conf.method = "HPDinterval",
               rhat = TRUE)
  ## save data (Y), simulation output (lfc), model coefficients
  return(list('data' = Y, 'sim_output' = lfc, 'estimate' = coef))
}

## define function of summary of B fits
smry <- function(sf) {
  ## estimated params
  ee <- sf$estimate
  ee <- ee[grepl("B",ee$term),]
  ## true params
  bb <- sf$sim_output$B_mat
  ## convert B_i_j_t from array to vec to match ee
  bvec <- NULL
  for(j in 1:n_species) {
    for(i in 1:n_species) {
      bvec <- c(bvec,bb[i,j,-c(1:(n_year/2+1),n_year)])
    }
  }
  ## drop intxns with 0's
  ii <- bvec!=0
  bvec <- bvec[ii]
  ee <- ee[ii,]
  ## chk for CI inclusion
  return(ee$conf.low < bvec & bvec < ee$conf.high)
}

## mean absolute percentage error
MASE <- function(sf) {
  ## estimated params
  ee <- sf$estimate
  ee <- ee[grepl("B",ee$term),]
  ## true params
  bb <- sf$sim_output$B_mat
  ## convert B_i_j_t from array to vec to match ee
  bvec <- NULL
  for(j in 1:n_species) {
    for(i in 1:n_species) {
      bvec <- c(bvec,bb[i,j,-c(1:(n_year/2+1),n_year)])
    }
  }
  ## drop intxns with 0's
  ii <- bvec!=0
  bvec <- bvec[ii]
  ee <- ee[ii,]
  return((sum(abs(bvec-ee$estimate)/ee$estimate))/length(bvec))
}

## function to estimate prop of non-converged params
pRhat <- function(sf, thresh=1.1) {
  ## estimated params
  ee <- sf$estimate
  ## Rhat values
  rh <- ee[!is.nan(ee$rhat),"rhat"]
  ## prop above threshold
  return(round(sum(rh>1.1)/length(ee$rhat), 2))
}

## better rnd
Re2prec <- function(x,map="round",prec=1) {
  ## 'map' can be round, floor, or ceiling
  ## 'prec' is nearest value (eg, 0.1 means to nearest tenth); default 1 gives normal behavior
  if(prec<=0) { stop("\"prec\" cannot be less than or equal to 0") }
  do.call(map,list(x/prec))*prec
}

##-------------
## user inputs
##-------------

## number of species/guilds
n_species <- 4
## number of years to simulate
n_year <- 60
## number of sites
nS <- 1
## number of MC simulations
n_sims <- 100

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

## define list of inputs for simulation
sim_list <- list(B0_init = B0_init,
                 B0_lfc = B0_lfc,
                 n_year = n_year,
                 n_site = nS,
                 var_QX = rev(seq(1, 4) / 40),
                 cov_QX = 0,
                 var_QB = 0.01,
                 cov_QB = 0,
                 dens_max = 3,
                 dens_min = -3)

## define list of inputs for simulations
fit_list <- list(topo = B0_lfc,
                 shared_r = matrix(1, n_species, nS),
                 mcmc_chain = 3,
                 mcmc_iter = 3000,
                 mcmc_warmup = 2000,
                 intervals = TRUE,
                 prob = 0.9)

##-----------
## sim & fit
##-----------

## start timer
timer_start <- proc.time()
## fit and save all expts
saved_output <- foreach(i=1:n_sims,
                        .export=c("sim_list", "fit_list"),
                        .packages=c("tvvarss","rstan","broom"),
                        .inorder=FALSE) %dopar% simfit(sim_list, fit_list)
## shut down workers
stopImplicitCluster()
## stop timer
(run_time_in_hrs <- round(((proc.time()-timer_start)/3600)["elapsed"], 1))
cat(run_time_in_hrs, file="run_time_in_hrs.txt")
## save output
save("saved_output",file="lfc_sim_fit_saved_output.RData")

## proportion of experiments where truth inside CI
props <- apply(sapply(saved_output, smry),1,sum)/n_sims

## plot summary
pdf("lfc_sim_fit_plots_2sites.pdf", height=6.5, width=6.5)

par(mfrow=c(n_species,n_species),
    mai=c(0.3,0.3,0.1,0.1),
    omi=c(0,0.3,0.3,0))
cnt <- 0
idx <- n_year/2-1
for(i in 1:n_species) {
  for(j in 1:n_species) {
    if(B0_init[j,i]!=0) {
      tmp <- props[1:idx+cnt*idx]
      plot.ts(tmp, ylim=c(0,1), xaxt="n", yaxt="n",
              cex.axis=0.9, col=ifelse(i==j,"blue","darkgreen"))
      axis(1,seq(0,Re2prec(idx,prec=10),by=10), las=1)
      axis(2,c(0,0.5,1), las=1)
      text(n_year/2-1, 0.05, round(mean(tmp),2), adj=c(1,0), cex=0.9)
      cnt <- cnt + 1
    } else {
      plot.ts(1:idx, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    }
    if(i==1) { mtext(side=3,j, line=1) }
    if(j==1) { mtext(side=2,i, las=1, line=3) }
  }
}

## proportion of params > Rhat by expt
pbad <- sapply(saved_output, pRhat)
par(mfrow=c(1,1), mai=c(0.9,0.9,0.1,0.1), omi=c(0,0,0,0))
hist(pbad, breaks=seq(0,20)/20, main="", xlab="Prop. of non-converged parameters")

dev.off()
