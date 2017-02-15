# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
tvvarss <- function(y, include_trend = TRUE, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 1000, mcmc_warmup = 500, mcmc_thin = 1, mcmc_chain = 3) {
  n_year = dim(y)[2]
  n_spp = dim(y)[3]
  n_site = ifelse(is.na(dim(y)[1]), 1, dim(y)[1])

  # vec B matrix, so we need to create a matrix of indices
  n_spp2 = n_spp*n_spp
  row_indices = rep(seq(1,n_spp), n_spp)
  col_indices = sort(row_indices)
  #Bindices = matrix(seq(1,n_spp2),n_spp,n_spp)

  # vector of 0s and 1s indicating whether element of vecB is on diagonal
  b_diag = rep(0, n_spp2)
  b_diag[seq(1,n_spp2,by=(n_spp+1))] = 1
  b_diag = b_diag + 1 # plus 1 because this is index in stan, 1 or 2

  if(is.null(x0)) x0 = y[,1,] # means on initial states,could also be set to 0
  est_trend = ifelse(include_trend,1,0); # estimate the trend

  if(is.null(shared_q)) shared_q = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_q = max(shared_q);

  if(is.null(shared_r)) shared_r = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_r = max(shared_r);

  if(is.null(shared_u)) shared_u = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_u = max(shared_u);

  stan_dir = find.package("tvvarss")
  model = paste0(stan_dir, "/exec/tvvarss.stan")

  model = paste0("exec/tvvarss.stan")
  datalist = list(n_year,
    n_spp,
    n_site,
    y,
    b_diag,
    x0,
    est_trend,
    shared_q,
    n_q,
    shared_r,
    n_r,
    shared_u,
    n_u,
    row_indices,
    col_indices)

  pars = c("sigma")
  mod = stan(data = datalist, pars = pars,
    file = model,
    chains = 1,
    iter = 500)

}
