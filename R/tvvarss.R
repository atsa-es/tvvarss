#' tvvarss is the primary function for fitting tvvarss models to time series data.
#'
#' @param y The data (array, with dimensions = site, year, species)
#' @param include_trend Whether to include time trends. Defaults to TRUE
#' @param x0 The location matrix (mean) of priors on initial states. Defaults to centered on observed data
#' @param shared_q Optional matrix (number of species x number of sites) with integers indicating which process variance parameters are shared. Defaults to shared process variances across sites, but unique to each species
#' @param shared_r Optional matrix (number of species x number of sites) with integers indicating which observation variance parameters are shared. Defaults to shared observation variances across sites, but unique to each species
#' @param shared_u Optional matrix (number of species x number of sites) with integers indicating which trend parameters are shared. Defaults to shared trends across sites, but unique to each species
#' @param mcmc_iter Number of MCMC iterations, defaults to 1000
#' @param mcmc_warmup Warmup / burn in phase, defaults to 500
#' @param mcmc_thin MCMC thin, defaults to 1
#' @param mcmc_chains MCMC chains, defaults to 3
#'

#'
#' @return an object of class 'stanfit'
#' @export
#'
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
  if(length(is.na(x0)) > 0) {
    x0 = matrix(0, n_site, n_spp)
  }

  est_trend = ifelse(include_trend,1,0); # estimate the trend

  if(is.null(shared_q)) shared_q = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_q = max(shared_q);

  if(is.null(shared_r)) shared_r = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_r = max(shared_r);

  if(is.null(shared_u)) shared_u = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  n_u = max(shared_u);

  # we also need to find values that aren't NA to allow for missing data
  spp_indices_pos = 0
  site_indices_pos = 0
  year_indices_pos = 0
  vec_y = 0
  count = 1
  for(i in 1:n_site) {
    for(j in 1:n_year) {
      for(k in 1:n_spp) {
        if(!is.na(y[i,j,k])) {
          spp_indices_pos[count] = k
          site_indices_pos[count] = i
          year_indices_pos[count] = j
          vec_y[count] = y[i,j,k]
          count = count + 1
        }
      }
    }
  }
  n_pos = length(spp_indices_pos)
  y = vec_y

  stan_dir = find.package("tvvarss")
  model = paste0(stan_dir, "/exec/tvvarss.stan")

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
    col_indices,
    spp_indices_pos,
    site_indices_pos,
    year_indices_pos,
    n_pos)

  pars = c("sigma_rw_pars", "resid_process_sd", "obs_sd", "B", "pred", "log_lik")
  if(include_trend) pars = c(pars, "u")

  mod = rstan::stan(data = datalist, pars = pars,
    file = model,
    chains = mcmc_chain,
    iter = mcmc_iter,
    warmup = mcmc_warmup,
    thin = mcmc_thin)
  return(mod)
}
