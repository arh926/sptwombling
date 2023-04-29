#' Spatiotemporal Gaussian Covariance Kernel
#'
#' Computes an adapted non-separable Gaussian covariance kernel
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param phis spatial range
#' @param phit temporal range
#' @param sig2 variance parameter
#' @keywords st_cov_gaussian
#' @examples
st_cov_gaussian <- function(delta = NULL, Delta = NULL, phis = NULL, phit = NULL, sig2 = NULL){
  Nt = nrow(delta)
  Ns = nrow(Delta)
  N = Ns * Nt

  A = ((phit * delta)^2 + 1)
  Sig = matrix(NA, nrow = N, ncol = N)
  for(i in 1 : Ns){
    for(j in i : Ns){
      Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
        Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = sig2/A *
        exp(- phis * Delta[i, j]^2 / A)
    }
  }
  Sig + 1e-10 * diag(N)
}
