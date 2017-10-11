#' Fit a TVVARSS model to multivariate time series data
#'
#' \code{tvvarss} is the primary function for fitting TVVARSS models data.
#'
#' @param y The data (array, with dimensions = site, year, species)
#' @param de_mean Whether or not to de_mean the process model; defaults to TRUE.
#'   For example, X_t+1 = B_t * (X_t - pred[X_t]) versus X_t+1 = B_t * (X_t).
#' @param topo Optional list matrix describing the presumed topology of the
#'   community. Pairwise interactions are specified as density-dependent ("dd"),
#'   top-down ("td"), bottom-up ("bu"), competitive/facilitative ("cf"), or
#'   absent ("zero").
#' @param dynamicB Logical indicator of whether to fit a dynamic B matrix that
#'   varies through time (or a static B matrix that does not); defaults to TRUE.
#' @param family Statistical distribution for the observation model, defaults to
#'  "gaussian".
#' @param x0 The location matrix (mean) of priors on initial states; defaults to
#'   centered on observed data.
#' @param shared_q Optional matrix (number of species x number of sites) with
#'   integers indicating which process variance parameters are shared; defaults
#'   to unique process variances for each species that are shared across sites.
#' @param shared_r Optional matrix (number of species x number of sites) with
#'   integers indicating which observation variance parameters are shared;
#'   defaults to unique observation variances for each species that are shared
#'   across sites.
#' @param process Vector that optionally maps sites to states. Defaults to each site as its own state
#' @param mcmc_iter Number of MCMC iterations, defaults to 1000
#' @param mcmc_warmup Warmup / burn in phase, defaults to 500
#' @param mcmc_thin MCMC thin, defaults to 1
#' @param mcmc_chain MCMC chains, defaults to 3
#'
#' @return an object of class 'stanfit'
#'
#' @import rstan
#' @import methods
#' @import Rcpp
#'
#' @export
tvvarss <- function(y, de_mean = TRUE, topo = NULL, dynamicB=TRUE, family="gaussian",
  x0 = NULL, shared_q = NULL, shared_r = NULL, process = NULL,
  mcmc_iter = 1000, mcmc_warmup = 500, mcmc_thin = 1, mcmc_chain = 3) {
  #@useDynLib tvvarss, .registration = TRUE
  include_trend = FALSE # not used as argument, but passed to STAN
  shared_u = NULL # not used as argument
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

  if(is.null(process)) process = seq(1, dim(y)[1])
  n_process = max(process)

  if(is.null(topo)) {
    # matrix constrained 0-1 on diagonal, and no constraints elsewhere
    topo = matrix("cf", n_spp, n_spp)
    diag(topo) = "dd"
  }

  # convert the character matrix to integers
  # zero = 1, td = 2, bu = 3, cf = 4
  B = matrix(1, n_spp, n_spp)
  if("zero"%in%topo) B[which(topo=="zero")] = 1
  if("td"%in%topo) B[which(topo=="td")] = 2
  if("bu"%in%topo) B[which(topo%in% c("bu"))] = 3
  if("dd"%in%topo) B[which(topo%in% c("dd"))] = 3
  if("cf"%in%topo) B[which(topo=="cf")] = 4

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
    x0 = matrix(0, n_process, n_spp)
  }

  est_trend = ifelse(include_trend,1,0); # estimate the trend
  demean = ifelse(de_mean, 1, 0);

  if(is.null(shared_q)) shared_q = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  shared_q = cbind(shared_q, 0);
  n_q = max(shared_q);

  if(is.null(shared_r)) shared_r = matrix(rep(1:n_spp,n_site), n_spp, n_site) # default to shared acros site,  unique by spp
  shared_r = cbind(shared_r, 0);
  n_r = max(shared_r);

  if(is.null(shared_u)) shared_u = matrix(rep(1:n_spp,n_process), n_spp, n_process) # default to shared acros site,  unique by spp
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
          site_indices_pos[count] = process[i]#i
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
    fit_dynamicB,
    n_process = n_process,
    process = process)

  pars = c("sigma_rw_pars", "resid_process_sd", "obs_sd", "B", "pred")
  if(include_trend) pars = c(pars, "u")

  mod = rstan::stan(data = datalist, pars = pars,
    file = model,
    chains = mcmc_chain,
    iter = mcmc_iter,
    warmup = mcmc_warmup,
    thin = mcmc_thin,
    control=list(adapt_delta=0.99, max_treedepth=20))
  return(mod)
}
