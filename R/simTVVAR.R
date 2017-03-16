#' Simulate the process component of a TVVARSS model
#' 
#' \code{simTVVAR} simulates the process (state) component of a TVVARSS model.
#'
#' \code{B0} can be used in one of two ways when simulating a TVVAR model:
#' \enumerate{
#'   \item an [n x n] \code{'matrix'} with a combination of \code{character} and
#'     \code{numeric} values in the off-diagonal elements; the diagonal should
#'     always contain 0's as density-dependence is implicit in this model. Use
#'     0 to indicate no interaction and the following \code{character} codes for
#'     interactions:
#'     \itemize{
#'       \item use \code{'td'} to indicate a "top-down" interaction
#'       \item use \code{'bu'} to indicate a "bottom-up" interaction
#'       \item use \code{'cf'} to indicate a "competitive/facilitative"
#'         interaction
#'     }
#'   \item an [n x n x (T+1)] \code{array} with actual values of B for each time
#'     step, including B0. This is useful for simulating multiple realizations
#'     of the same process.
#' }
#' See 'Examples' for details on formatting \code{B0}.
#'
#' @param \code{B0} A matrix describing the topology of the food web (see 'Details').
#' @param \code{TT} Number of time steps to simulate.
#' @param \code{var_QX} Scalar or vector of variances for process errors of states.
#' @param \code{cov_QX} Covariance, if any, of the process errors of the states; if \code{cov_QX} > 0, then \code{var_QX} must be a scalar. 
#' @param \code{var_QB} Scalar or vector of variances for process errors of \strong{B}.
#' @param \code{cov_QB} Covariance, if any, of process errors of \strong{B}; if \code{cov_QB} > 0, then \code{var_QB} must be a scalar.
#' @param \code{QQ_XX} [optional] Specify the explicit form for the var-cov matrix \strong{Q} of the process errors. 
#' @param \code{QQ_BB} [optional] Specify the explicit form for the var-cov matrix \strong{Q} of \strong{B}. 
#' @param \code{X0} [optional] Specify vector of initial states; \code{nrow(X0)} must equal \code{nrow(B0)}. 
#' @param \code{CC} [optional] Specify matrix of covariate effects on states. 
#' @param \code{cc} [optional] Specify matrix of covariates. 
#'
#' @return A list with the following components:
#' \describe{
#' \item{\code{B_mat}}{An array of the \strong{B} matrix over time; \code{dim(B_mat) = c(n,n,T+1)}.}
#' \item{\code{WW_BB}}{The process errors for \strong{B}.}
#' \item{\code{QQ_BB}}{Variance-covariance matrix of the process errors for \strong{B}.}
#' \item{\code{states}}{A matrix of the states over time; \code{dim(states) = c(n,T+1)}.}
#' \item{\code{WW_XX}}{The process errors (innovations) for the states.}
#' \item{\code{QQ_XX}}{Variance-covariance matrix of the process errors for the states.}
#' }
#'
#' @examples
#' set.seed(123)
#' ## number of time steps
#' TT <- 30
#' ## number of spp/guilds
#' nn <- 4
#' 
#' ## CASE 1: linear food chain
#' B0_lfc <- matrix(list(0),nn,nn)
#' for(i in 1:(nn-1)) {
#'   B0_lfc[i,i+1] <- "td"
#'   B0_lfc[i+1,i] <- "bu"
#' }
#' ## inspect B0
#' B0_lfc
#' ## simulate & plot states
#' lfc <- simTVVAR(B0_lfc,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#' matplot(t(lfc$states),type="l")
#' 
#' ## CASE 2: 1 consumer & n-1 producers
#' B0_cp <- matrix(list("cf"),nn,nn)
#' B0_cp[1:(nn-1),nn] <- "td"
#' B0_cp[nn,1:(nn-1)] <- "bu"
#' diag(B0_cp) <- 0
#' ## inspect B0
#' B0_cp
#' ## simulate & plot states
#' cp <- simTVVAR(B0_cp,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#' matplot(t(cp$states),type="l")
#'
#' ## simulate a second realization of CASE 2 using same B
#' cp2 <- simTVVAR(cp$B_mat,TT,var_QX=rev(seq(1,4)/40),cov_QX=0,var_QB=0.05,cov_QB=0)
#' matplot(t(cp2$states),type="l")
#' 
#' @export
simTVVAR <- function(B0,TT,var_QX,cov_QX,var_QB,cov_QB=0,QQ_XX=NULL,QQ_BB=NULL,X0=NULL,CC=NULL,cc=NULL) {
  if(class(B0)!="matrix" & class(B0)!="array") {
    stop("'B0' must be an [n x n] matrix or [n x n x T] array.\n\n")
  } else {
    if(length(dim(B0))<2 | length(dim(B0))>3 | dim(B0)[1]!=dim(B0)[2]) {
      stop("'B0' must be an [n x n] matrix or [n x n x T] array.\n\n")
    }
  }
  ## number of spp/guilds
  nn <- dim(B0)[1]
  ## if no var-cov matrix for proc errors of states was passed, create one
  if(is.null(QQ_XX)) {
    ## var-cov matrix for proc errors of states
    if(length(var_QX)>1 & cov_QX != 0) {
      stop("If 'var_QX is a vector, then 'cov_QX' must be 0.\n\n", call.=FALSE)
    }
    QQ_XX <- matrix(cov_QX,nn,nn)
    diag(QQ_XX) <- var_QX	
  } else {
    if(!all(dim(QQ_XX)==nn)) {
      stop("'QQ_XX' must be an [n x n] matrix.\n\n")
    }
  }
  ## if no var-cov matrix for proc errors of BB was passed, create one
  if(is.null(QQ_BB)) {
    ## var-cov matrix for proc errors of BB
    if(length(var_QB)>1 & cov_QB != 0) {
      stop("If 'var_QB is a vector, then 'cov_QB' must be 0.\n\n", call.=FALSE)
    }
    QQ_BB <- matrix(cov_QB,nn*nn,nn*nn)
    diag(QQ_BB) <- var_QB
  } else {
    if(!all(dim(QQ_BB)==nn)) {
      stop("'QQ_BB' must be an [nn x nn] matrix.\n\n")
    }
  }
  ## STATES
  XX <- matrix(NA,nn,TT+1)
  ## initial states
  if(is.null(X0)) {
    XX[,1] <- matrix(runif(nn,-5,5),nn,1)
  } else {
    XX[,1] <- X0
  }
  ## proc errors for states
  WW_XX <- t(mvrnorm(TT,matrix(0,nn,1),QQ_XX))
  ## BB
  if(length(dim(B0))==2) { ## then B0 is food web topology
    ## initial BB
    BB <- array(0,c(nn,nn,TT+1))
    ## build B0 based on topology
    ## fill in diagonal
    diag(BB[,,1]) <- plogis(rnorm(nn,0,1))
    ## fill in top-down interactions
    i_td <- sapply(B0,function(x) {x=="td"})
    BB[,,1][i_td] <- -plogis(rnorm(sum(i_td),0,1))
    ## fill in bottom-up interactions
    i_bu <- sapply(B0,function(x) {x=="bu"})
    BB[,,1][i_bu] <- plogis(rnorm(sum(i_bu),0,1))
    ## fill in competitive/facilitative interactions
    i_cf <- sapply(B0,function(x) {x=="cf"})
    BB[,,1][i_cf] <- plogis(rnorm(sum(i_cf),0,1))*2-1
    ## proc errors for BB
    WW_BB <- t(mvrnorm(TT,matrix(0,nn*nn,1),QQ_BB))
    WW_BB[which(BB[,,1]==0),] <- 0
  } else { ## B0 is a passed [n x n x T] array of interaction strengths
    BB <- B0
    WW_BB <- NULL
  }
  ## covariates, if missing
  if(is.null(CC)) {
    CC <- matrix(0,nn,1)
    cc <- matrix(0,1,TT+1)
  }
  ## evolutions
  for(t in 1:TT+1) {
    ## evolution of BB
    if(length(dim(B0))==2) {
      ## constrain diagonals to [0,1]
      diag(BB[,,t]) <- plogis(qlogis(diag(BB[,,t-1])) + diag(matrix(WW_BB[,t-1],nn,nn)))
      ## constrain off-diag effects to [-1,1]
      BB[,,t][i_td] <- -plogis(qlogis(-BB[,,t-1][i_td]) + WW_BB[i_td,t-1])
      BB[,,t][i_bu] <- plogis(qlogis(BB[,,t-1][i_bu]) + WW_BB[i_bu,t-1])
      BB[,,t][i_cf] <- plogis(qlogis((BB[,,t-1][i_cf]+1)/2) + WW_BB[i_cf,t-1])*2-1
    }
    ## state
    XX[,t] <- BB[,,t]%*%(XX[,t-1,drop=FALSE] - CC%*%cc[,t-1]) + CC%*%cc[,t] + WW_XX[,t-1]
  }
  return(list(B_mat=BB,WW_BB=WW_BB,QQ_BB=QQ_BB,states=XX,WW_XX=WW_XX,QQ_XX=QQ_XX))
} ## end function
