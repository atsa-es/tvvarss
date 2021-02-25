#' Simulate the process component of a TVVARSS model
#'
#' \code{simTVVAR} simulates the process (state) component of a TVVARSS model.
#'
#' \code{Bt} can be used in one of two ways when simulating a TVVAR model:
#' \enumerate{
#'   \item An \eqn{n x n} \code{matrix} with initial numeric values of B (i.e., B0).
#'     If \code{QQ_BB = matrix(0, n, n)} then, a time-invariant (MARSS) model is
#'     simulated based on these values.
#'   \item An \eqn{n x n x (T+1)} \code{array} with actual values of B for each time
#'     step, including B0. This is useful for simulating multiple realizations
#'     of the same process.
#' }
#' \code{topo} can be used to specify the food web topology by passing an
#' \eqn{n x n} \code{matrix} with a combination of \code{character} and
#' \code{numeric} values in the off-diagonal elements; the diagonal should
#' always contain \code{"dd"} as density-dependence is implicit in this
#' model. Use 0 or "zero" to indicate no interaction and the following
#' \code{character} codes for ecological interactions:
#'     \itemize{
#'       \item \code{"td"} to indicate a top-down interaction
#'       \item \code{"bu"} to indicate a bottom-up interaction
#'       \item \code{"cf"} to indicate a competitive/facilitative
#'         interaction
#'     }

#' See 'Examples' for details on formatting \code{B0}.
#'
#' @param Bt A matrix describing the topology of the food web (see 'Details').
#'   If \code{Bt == NULL}, then the food web topology must be specified and
#'   passed as \code{topo}. See 'Details'.
#' @param topo Optional list matrix describing the presumed topology of the
#'   community. Pairwise interactions are specified as density-dependent ("dd"),
#'   top-down ("td"), bottom-up ("bu"), competitive/facilitative ("cf"), or
#'   absent ("zero"). If specified, pairwise interactions will be constrained
#'   in an approporiate manner (e.g., top-down effects are between -1 and 0).
#' @param TT Number of time steps to simulate.
#' @param var_QX Scalar or vector of variances for process errors of states.
#' @param cov_QX Covariance, if any, of the process errors of the states; if \code{cov_QX} > 0, then \code{var_QX} must be a scalar.
#' @param var_QB Scalar or vector of variances for process errors of \strong{B}.
#' @param cov_QB Covariance, if any, of process errors of \strong{B}; if \code{cov_QB} > 0, then \code{var_QB} must be a scalar.
#' @param QQ_XX Optionally specify the explicit form for the var-cov matrix \strong{Q} of the process errors of the states.
#' @param QQ_BB Optionally specify the explicit form for the var-cov matrix \strong{Q} of the process errors of \strong{B}.
#' @param X0 Optionally specify vector of initial states; \code{nrow(X0)} must equal \code{nrow(Bt)}.
#' @param CC Optionally specify matrix of covariate effects on states.
#' @param cc Optionally specify matrix of covariates.
#'
#' @return A list with the following components:
#' \describe{
#' \item{\code{B_mat}}{An array of the \strong{B} matrix over time; \code{dim(B_mat) = c(n,n,T+1)}.}
#' \item{\code{WW_BB}}{The process errors for \strong{B}; \code{dim(WW_BB) = c(n^2,T)}.}
#' \item{\code{QQ_BB}}{Variance-covariance matrix of the process errors for \strong{B}; \code{dim(QQ_BB) = c(n^2,n^2)}.}
#' \item{\code{states}}{A matrix of the states over time; \code{dim(states) = c(n,T+1)}.}
#' \item{\code{WW_XX}}{The process errors (innovations) for the states; \code{dim(WW_XX) = c(n,T)}.}
#' \item{\code{QQ_XX}}{Variance-covariance matrix of the process errors for the states; \code{dim(QQ_XX) = c(n,n)}.}
#' \item{\code{call}}{The function call as returned by \code{match.call()}.}
#' }
#'
#' @examples
# set.seed(123)
# ## number of time steps
# TT <- 30
# ## number of spp/guilds
# nn <- 4
# ## CASE 1: linear food chain; starting values are random
# B0_lfc <- matrix(list(0),nn,nn)
# diag(B0_lfc) <- "dd"
# for(i in 1:(nn-1)) {
#   B0_lfc[i,i+1] <- "td"
#   B0_lfc[i+1,i] <- "bu"
# }
# ## inspect B0
# B0_lfc
# ## simulate & plot states
# lfc <- simTVVAR(Bt=NULL,topo=B0_lfc,TT=TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
# matplot(t(lfc$states),type="l")
#
# ## CASE 2: 1 consumer & n-1 producers; starting values are random
# B0_cp <- matrix(list("cf"),nn,nn)
# B0_cp[1:(nn-1),nn] <- "td"
# B0_cp[nn,1:(nn-1)] <- "bu"
# diag(B0_cp) <- "dd"
# ## inspect B0
# B0_cp
# ## simulate & plot states
# cp <- simTVVAR(Bt=NULL,topo=B0_lfc,TT=TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
# matplot(t(cp$states),type="l")
#
# ## simulate a second realization of CASE 2 using same B
# cp2 <- simTVVAR(Bt=cp$B_mat,topo=B0_lfc,TT=TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#'
#' @importFrom stats plogis qlogis rnorm runif
#' @importFrom MASS mvrnorm
#'
#' @export
#'
simTVVAR <- function(Bt=NULL, topo=NULL, TT, var_QX, cov_QX, var_QB, cov_QB = 0,
                     QQ_XX = NULL, QQ_BB = NULL, X0 = NULL,
                     CC = NULL, cc = NULL) {
  if(!is.null(Bt)) {
    if(class(Bt)[1] %in%c("matrix","array")==FALSE) {
      stop("'Bt' must be an [n x n] matrix or [n x n x T] array of interaction strengths. Otherwise, it must be set to NULL with 'topo' passed as well.\n\n")
    }
    if(length(dim(Bt)) < 2 | length(dim(Bt)) > 3 | dim(Bt)[1] != dim(Bt)[2]) {
        stop("'Bt' must be an [n x n] matrix or [n x n x T] array of interaction strengths. Otherwise, it must be set to NULL with 'topo' passed as well.\n\n")
      }
    ## number of spp/guilds
    nn <- dim(Bt)[1]
  } else {
    nn <- dim(topo)[1]
  }
  ## if no var-cov matrix for proc errors of states was passed, create one
  if (is.null(QQ_XX)) {
    ## var-cov matrix for proc errors of states
    if (length(var_QX) > 1 & cov_QX != 0) {
      stop("If 'var_QX is a vector, then 'cov_QX' must be 0.\n\n", call. = FALSE)
    }
    QQ_XX <- matrix(cov_QX, nn, nn)
    diag(QQ_XX) <- var_QX
  } else {
    if (!all(dim(QQ_XX) == nn)) {
      stop("'QQ_XX' must be an [n x n] matrix.\n\n")
    }
  }
  ## if no var-cov matrix for proc errors of BB was passed, create one
  if (is.null(QQ_BB)) {
    ## var-cov matrix for proc errors of BB
    if (length(var_QB) > 1 & cov_QB != 0) {
      stop("If 'var_QB is a vector, then 'cov_QB' must be 0.\n\n", call. = FALSE)
    }
    QQ_BB <- matrix(cov_QB, nn * nn, nn * nn)
    diag(QQ_BB) <- var_QB
  } else {
    if (!all(dim(QQ_BB) == nn * nn)) {
      stop("'QQ_BB' must be an [nn x nn] matrix.\n\n")
    }
  }
  ## STATES
  XX <- matrix(NA, nn, TT + 1)
  ## initial states
  if (is.null(X0)) {
    XX[,1] <- matrix(runif(nn, -1, 1), nn, 1)
  } else {
    XX[,1] <- X0
  }
  ## proc errors for states
  WW_XX <- t(MASS::mvrnorm(TT, matrix(0, nn, 1), QQ_XX))
  ## BB
  BB <- array(0, c(nn, nn, TT + 1))
  if (length(dim(Bt)) == 2) { ## then Bt is an [n x n] init matrix
    ## initial BB
    BB[,,1] <- Bt
  }
  if(is.null(topo)) { ## matrix constrained 0:1 on diagonal, and -1:1 elsewhere
    topo <- matrix("cf", nn, nn)
    diag(topo) <- "dd"
  }
  ## find top-down interactions
  i_td <- sapply(topo, function(x) { x == "td" })
  ## find bottom-up interactions
  i_bu <- sapply(topo, function(x) { x == "bu" })
  ## find competitive/facilitative interactions
  i_cf <- sapply(topo, function(x) { x == "cf" })
  ## if no Bt, make one
  if(is.null(Bt)) {
    diag(BB[, , 1]) <- plogis(rnorm(nn, 0, 1))
    BB[, , 1][i_td] <- -plogis(rnorm(sum(i_td), 0, 1))
    BB[, , 1][i_bu] <- plogis(rnorm(sum(i_bu), 0, 1))
    BB[, , 1][i_cf] <- plogis(rnorm(sum(i_cf), 0, 1)) * 2 - 1
  }
  ## proc errors for BB
  WW_BB <- t(MASS::mvrnorm(TT, matrix(0, nn*nn, 1), QQ_BB))
  WW_BB[which(BB[,,1] == 0), ] <- 0
  if (length(dim(Bt)) == 3) { ## then Bt is [n x n x T] array of interactions
    BB <- Bt
    WW_BB <- NULL
  }
  ## covariates, if missing
  if (is.null(CC)) {
    CC <- matrix(0, nn, 1)
    cc <- matrix(0, 1, TT + 1)
  }
  ## evolutions
  for (t in 1:TT + 1) {
    ## BB
    if (length(dim(Bt)) != 3) {
      ## constrain diagonals to [0,1]
      diag(BB[,,t]) <- plogis(qlogis(diag(BB[,,t-1])) + diag(matrix(WW_BB[,t-1],nn,nn)))
      ## constrain top-down effects to [-1,0]
      BB[,,t][i_td] <- -plogis(qlogis(-BB[,,t-1][i_td]) + WW_BB[i_td,t-1])
      ## constrain bottom-up effects to [0,1]
      BB[,,t][i_bu] <- plogis(qlogis(BB[,,t-1][i_bu]) + WW_BB[i_bu,t-1])
      ## constrain comp-facil effects to [-1,1]
      BB[,,t][i_cf] <- plogis(qlogis((BB[,,t-1][i_cf] + 1) / 2) + WW_BB[i_cf,t-1]) * 2 - 1
    }
    ## state
    XX[,t] <- BB[,,t] %*% (XX[,t-1,drop = FALSE] - CC %*% cc[,t-1]) + CC %*% cc[,t] + WW_XX[,t-1]
  }
  return(list(
    B_mat = BB,
    WW_BB = WW_BB,
    QQ_BB = QQ_BB,
    states = XX,
    WW_XX = WW_XX,
    QQ_XX = QQ_XX,
    call = match.call()
  ))
} ## end function
