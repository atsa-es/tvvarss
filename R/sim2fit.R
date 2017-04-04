#' Simulate TVVAR model and add observation error
#'
#' \code{sim2fit} adds observation error to a simulated TVVAR process and
#'   converts it to a form suitable for fitting with \code{tvvarss}.
#'
#' This is a helper function that takes a fitted \code{simTVVAR} object and
#' simulates multiple realizations of the process before adding Gaussian
#' obsveration errors.
#'
#' @param obj A fitted \code{simTVVAR} object.
#' @param n_sim The number of realizations of the TVVAR process.
#' @param sd The standard deviation of the Gaussian observation errors.
#'   Can be set to 0 for no observation error.
#'
#' @return An array with dimensions \code{c(n_sim, TT, n_spp)}.
#'
#' @examples
#' set.seed(123)
#' ## number of time steps
#' TT <- 30
#' ## number of spp/guilds
#' nn <- 4
#' ## CASE 1: linear food chain
#' B0_lfc <- matrix(list(0),nn,nn)
#' for(i in 1:(nn-1)) {
#'   B0_lfc[i,i+1] <- "td"
#'   B0_lfc[i+1,i] <- "bu"
#' }
#' ## simulate process
#' lfc <- simTVVAR(B0_lfc, TT, var_QX=rev(seq(1,4)/40), cov_QX=0, var_QB=0.05, cov_QB=0)
#' ## create data array with 3 realizations of the process
#' dat <- sim2fit(lfc, 3)
#'
#' @export
sim2fit <- function(obj, n_sim, sd=0.1) {
  if(n_sim < 1 | round(n_sims) - n_sims != 0) {
    stop("'n_sims' must be a positive integer.")
  }
  if(sd < 0) {
    stop("'sd' must be non-negative.")
  }
  nC <- ncol(obj$states)-1
  nR <- nrow(obj$states)
  yy <- array(NA,c(n_sim,nC,nR))
  for(i in 1:n_sim) {
    y[i,,] <- t(eval(as.expression(obj$call))$states[,-1]) + matrix(rnorm(nC*nR,0,sd),nC,nR)
  }
  return(yy)
}
