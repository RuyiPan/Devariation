#' Devariation-adjusted data via low-rank extraction
#'
#' Constructs the devariation-adjusted data by separating the observed matrix
#' into (i) a dominant structured component and (ii) a residual
#' component. The residual is interpreted as the
#' devariation-adjusted data, where large-scale structured variation is removed.
#'
#' The structured component \code{S} is obtained by extracting dominant
#' low-rank variation from \code{X}. If \code{lambda} is not provided, it is
#' set using an estimated noise scale \code{tauX} (via \code{estim_sigma()})
#' and the high-dimensional scaling factor \eqn{\sqrt{n} + \sqrt{p}} for an
#' \eqn{n \times p} matrix \code{X}.
#'
#' @param X A numeric matrix (e.g., observed data containing structured signal and noise).
#' @param lambda Optional nonnegative tuning parameter controlling the strength of
#'   devariation (how much dominant structure is removed). If \code{NULL}, the
#'   function estimates \code{tauX} and sets
#'   \eqn{\lambda = \tau_X(\sqrt{n} + \sqrt{p})}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{dX}{The devariation-adjusted data matrix, defined as \code{dX = X - S}.}
#'   \item{lambda}{The tuning parameter used for devariation.}
#'   \item{tauX}{Estimated noise scale. If \code{lambda} is supplied, \code{tauX = lambda / (\sqrt{n} + \sqrt{p})}.}
#' }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 80), 100, 80)
#' out <- dev(X)
#' dim(out$dX)
#'
#' # User-specified devariation strength
#' out2 <- dev(X, lambda = 5)
#'
#' @export


dev <- function(X, lambda=NULL) {

  svdx <- svd(X)
  if (is.null(lambda)) {
    tauX <- denoiseR::estim_sigma(X, method="MAD")
    lambda <- tauX*(sqrt(dim(X)[1]) + sqrt(dim(X)[2]))
  }else{
    tauX <- lambda/(sqrt(dim(X)[1]) + sqrt(dim(X)[2]))
  }
  S=svdx$u%*%diag(pmax(svdx$d - lambda,0))%*%t(svdx$v)
  # rank=sum(svdx$d >= lambda)
  dX = X - S
  return(list(dX = dX, lambda=lambda, tauX=tauX))
}






#' Devariation RV (dRV) coefficient
#'
#' Computes the devariation RV (dRV) coefficient between two data matrices.
#' The procedure first removes dominant structured variation from each matrix
#' using the \code{dev()} function, producing devariation-adjusted data
#' \eqn{dX} and \eqn{dY}. The RV coefficient is then computed on the
#' devariation-adjusted matrices to quantify residual cross-structure association.
#'
#' This approach aims to mitigate the influence of strong individual
#' low-rank components and assess association after devariation.
#'
#' @param X A numeric matrix (n × p).
#' @param Y A numeric matrix (n × q). Must have the same number of rows as \code{X}.
#' @param lambdaX Optional tuning parameter controlling the strength of
#'   devariation applied to \code{X}. If \code{NULL}, it is automatically
#'   determined within \code{dev()}.
#' @param lambdaY Optional tuning parameter controlling the strength of
#'   devariation applied to \code{Y}. If \code{NULL}, it is automatically
#'   determined within \code{dev()}.
#'
#' @return A list containing:
#' \describe{
#'   \item{drvRes}{The devariation RV coefficient computed from
#'   the devaried matrices \code{dX} and \code{dY}.}
#' }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 80), 100, 80)
#' Y <- matrix(rnorm(100 * 60), 100, 60)
#' out <- drv(X, Y)
#' out$drvRes
#' out$drvRes$p.value
#'
#' @seealso \code{\link{dev}}, \code{\link{coeffRV}}
#'
#' @export
drv <- function(X, Y, lambdaX=NULL, lambdaY=NULL) {

  dX = dev(X, lambdaX)[[1]]
  dY = dev(Y, lambdaY)[[1]]

  return(list(drvRes = FactoMineR::coeffRV(dX,dY)))
}
