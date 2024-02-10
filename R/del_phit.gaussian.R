#' Partial derivative with respect to $\eqn{\phi_t}$ of the spatiotemporal Gaussian Covariance Kernel
#'
#' Computes an adapted non-separable Gaussian covariance kernel
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param phis spatial range
#' @param phit temporal range
#' @param sig2 variance parameter
#' @keywords st_cov_matern1
#' @examples
#'
del.phit.st_cov_gaussian <- function(delta = NULL,
                                     Delta = NULL,
                                     phis = NULL,
                                     phit = NULL,
                                     sig2 = NULL){
  Nt <- nrow(delta)
  Ns <- nrow(Delta)
  N <- Ns * Nt

  A <- ((phit * delta)^2 + 1)
  Sig <- matrix(NA, nr = N, nc = N)
  for(i in 1:Ns){
    for(j in i:Ns){
      Sig[(i - 1) * Nt + 1:Nt, (j - 1) * Nt + 1:Nt] =
        Sig[(j - 1) * Nt + 1:Nt, (i - 1) * Nt + 1:Nt] =
        - 2 * sig2/A^2 * phit * delta^2 * exp(-phis * Delta[i,j]^2/A) * (1 - phis * Delta[i, j]/A)
    }
  }
  Sig
}
