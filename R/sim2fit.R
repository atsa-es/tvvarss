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
#' @param n_sims The number of realizations of the TVVAR process.
#' @param sd The standard deviation of the Gaussian observation errors.
#'   Can be set to 0 for no observation error.
#' @param new_real If n_sims > 1, logical indicator of whether to base the new
#'   observations on a new realization of the TVVAR process.
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
#' topo <- matrix(list(0),nn,nn)
#' for(i in 1:(nn-1)) {
#'   topo[i,i+1] <- "td"
#'   topo[i+1,i] <- "bu"
#' }
#' ## simulate process
#' lfc <- simTVVAR(topo, TT, var_QX=rev(seq(1,4)/40), cov_QX=0, var_QB=0.05, cov_QB=0)
#' ## create data array with 3 realizations of the process
#' dat <- sim2fit(lfc, 3)
#'
#' @importFrom stats rnorm
#'
#' @export
sim2fit <- function(obj, n_sims, sd=0.1, new_real=TRUE) {
  if(n_sims < 1 | round(n_sims) - n_sims != 0) {
    stop("'n_sims' must be a positive integer.")
  }
  if(sd < 0) {
    stop("'sd' must be non-negative.")
  }
  nC <- ncol(obj$states)-1
  nR <- nrow(obj$states)
  yy <- array(NA,c(n_sims,nC,nR))
  for(i in 1:n_sims) {
    if(new_real) {
      yy[i,,] <- t(eval(as.expression(obj$call))$states[,-1]) + matrix(rnorm(nC*nR,0,sd),nC,nR)
    } else {
      yy[i,,] <- t(obj$states[,-1]) + matrix(rnorm(nC*nR,0,sd),nC,nR)
    }
  }
  return(yy)
}
