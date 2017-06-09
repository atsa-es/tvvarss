devtools::install_github("nwfsc-timeseries/tvvarss")
library(tvvarss)
library(tseries)
library(broom)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 4
n_year <- 100
n_site <- 1
n_simulations = 100

#' ## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_species,n_species)
diag(B0_lfc) <- "dd"
for(i in 2:(n_species-1)) {
  B0_lfc[i,i+1] <- "td"
  B0_lfc[i+1,i] <- "bu"
}

saved_output = list()

for(i in 1:n_simlations) {

# simulate states -- nest in loop to assure stationarity
stationary = FALSE

while(stationary == FALSE) {
lfc <-
  simTVVAR(
    Bt = NULL,
    topo = B0_lfc,
    TT = n_year,
    var_QX = rev(seq(1, 4) / 40),
    cov_QX = 0,
    var_QB = 0.01,
    cov_QB = 0
  )

dat <- sim2fit(lfc, n_site, sd=0.1)

# drop the burn in (1st half of the states/B matrix)
Bmat = lfc$B_mat[,,-c(1:(n_year/2 + 1))]
Y = array(0, dim = c(dim(dat)[1], dim(dat)[2]/2, dim(dat)[3]))
for(i in 1:n_site) {
  Y[i,,] = dat[i,-c(1:(n_year/2)),]
}

trend_test = 0
# Apply ADF and PP tests to evaluate stationarity
for (i in 1:n_species) {
  for (j in 1:n_site) {
    trend_test = trend_test + ifelse(adf.test(Y[j,,i])$p.value < 0.05, 0, 1) +
      ifelse(adf.test(Y[j,,i])$p.value < 0.05, 0, 1) +
      ifelse(kpss.test(Y[j,,i])$p.value > 0.05, 0, 1) +
      ifelse(coef(summary(lm(c(NA,log(abs(diff(Y[j,,i])))) ~ seq(1,dim(Y)[2]))))[2,4] > 0.05, 0, 1)
  }
}
if(trend_test == 0) stationary=TRUE
}

# data has been generated -- now fit the model using tvvarss()
fitted_model = tvvarss(y = Y, topo = B0_lfc,
  shared_r = matrix(1, n_species, n_site), mcmc_iter = 3000, mcmc_warmup = 2000)
# tidy the stanfit object to save space -- don't save every mcmc draw
coef = tidy(fitted_model, intervals = TRUE, prob = 0.9)

# save data (Y), simulation output (lfc), model coefficients
saved_output[[i]] = list('data' = Y,
  'sim_output' = lfc, 'estimate' = coef)
}

save(saved_output, file = "output_linear_chain.Rdata")
