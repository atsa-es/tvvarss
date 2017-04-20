#' tvvarss is the primary function for fitting tvvarss models to time series data.
#'
#' @param y The data (array, with dimensions = site, year, species)
#' @param include_trend Whether to include time trends. Defaults to TRUE
#' @param de_mean Whether or not to de_mean the process model, defaults to TRUE. For example, X_t+1 = B * (X_t - pred[X_t]) versus X_t+1 = B * (X_t)
#' @param B The list matrix describing optionally whether elements are 'zero', top-down ('td'), bottom up ('bu'), competitive-facilitative ('cf'), or density dependent ('dd)
#' @param x0 The location matrix (mean) of priors on initial states. Defaults to centered on observed data
#' @param shared_q Optional matrix (number of species x number of sites) with integers indicating which process variance parameters are shared. Defaults to shared process variances across sites, but unique to each species
#' @param shared_r Optional matrix (number of species x number of sites) with integers indicating which observation variance parameters are shared. Defaults to shared observation variances across sites, but unique to each species
#' @param shared_u Optional matrix (number of species x number of sites) with integers indicating which trend parameters are shared. Defaults to shared trends across sites, but unique to each species
#' @param mcmc_iter Number of MCMC iterations, defaults to 1000
#' @param mcmc_warmup Warmup / burn in phase, defaults to 500
#' @param mcmc_thin MCMC thin, defaults to 1
#' @param mcmc_chain MCMC chains, defaults to 3
#' @param family observation model, defaults to 'gaussian'
#' @param dynamicB boolean, whether to fit a dynamic B matrix that varies through time (or a static B matrix that doesn't), defaults to TRUE
#'
#' @return an object of class 'stanfit'
#' @export
#' @import rstan
#' @import methods
#' @import Rcpp
#'
tvvarss <- function(y, include_trend = TRUE, de_mean = TRUE, B = NULL, x0 = NULL, shared_q = NULL, shared_r = NULL, shared_u = NULL, mcmc_iter = 1000, mcmc_warmup = 500, mcmc_thin = 1, mcmc_chain = 3, family="gaussian", dynamicB=TRUE) {

  dist = c("gaussian", "binomial", "poisson", "gamma", "lognormal")
  family = which(dist==family)

  if(length(dim(y))==3) {
    # multiple sites in array, site=1st d
    n_year = dim(y)[2]
    n_spp = dim(y)[3]
    n_site = ifelse(is.na(dim(y)[1]), 1, dim(y)[1])
  }
  if(length(dim(y)) == 2) {
    # matrix for 1 site, coerce to 3d
    n_year = dim(y)[1]
    n_site = 1
    n_spp = dim(y)[2]
    y_new = array(0, dim = c(2, dim(y)[1], dim(y)[2]))
    y_new[1,,] = y # make the first array data, rest = NA
    y = y_new
  }

  #
  if(is.null(B)) {
    # matrix constrained 0-1 on diagonal, and no constraints elsewhere
    B = matrix("cf", n_spp, n_spp)
    diag(B) = "dd"
  }

  # convert the character matrix to integers
  # zero = 1, td = 2, bu = 3, cf = 4
  if("zero"%in%B) B[which(B=="zero")] = 1
  if("td"%in%B) B[which(B=="td")] = 2
  if("bu"%in%B) B[which(B%in% c("bu"))] = 3
  if("dd"%in%B) B[which(B%in% c("dd"))] = 3
  if("cf"%in%B) B[which(B=="cf")] = 4
  class(B) = "numeric"

  b_limits = matrix(0, 4, 2)
  b_limits[1,] = c(0,0)
  b_limits[2,] = c(-1, 0)
  b_limits[3,] = c(0, 1)
  b_limits[4,] = c(-1, 1)
  # vec B matrix, so we need to create a matrix of indices
  n_spp2 = n_spp*n_spp
  row_indices = rep(seq(1,n_spp), n_spp)
  col_indices = sort(row_indices)
  b_indices = c(B)
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
  demean = ifelse(de_mean, 1, 0);

  if(is.null(shared_q)) shared_q = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  shared_q = cbind(shared_q, 0);
  n_q = max(shared_q);

  if(is.null(shared_r)) shared_r = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  shared_r = cbind(shared_r, 0);
  n_r = max(shared_r);

  if(is.null(shared_u)) shared_u = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  shared_u = cbind(shared_u, 0);
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
  y_int = round(y)

  stan_dir = find.package("tvvarss")
  model = paste0(stan_dir, "/exec/tvvarss.stan")

  fit_dynamicB = as.integer(dynamicB) # convert 0 or 1

  datalist = list(n_year,
    n_spp,
    n_site,
    y,
    b_diag,
    x0,
    est_trend,
    demean,
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
    n_pos,
    y_int,
    family,
    b_indices,
    b_limits,
    fit_dynamicB)

  pars = c("sigma_rw_pars", "resid_process_sd", "obs_sd", "B", "pred")
  if(include_trend) pars = c(pars, "u")

  mod = rstan::stan(data = datalist, pars = pars,
    file = model,
    chains = mcmc_chain,
    iter = mcmc_iter,
    warmup = mcmc_warmup,
    thin = mcmc_thin)
  return(mod)
}
