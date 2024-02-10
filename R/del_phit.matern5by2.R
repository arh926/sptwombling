#' Partial derivative with respect to $\eqn{\phi_t}$ of the spatiotemporal Mat\'ern Covariance Kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' Computes an adapted non-separable Mat\'ern covariance kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param phis spatial range
#' @param phit temporal range
#' @param sig2 variance parameter
#' @keywords st_cov_matern1
#' @examples
#'
del.phit.st_cov_matern2 <- function(delta = NULL,
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
        - sig2/A^(3/2) * phit * delta^2 * exp(-sqrt(5) * phis * Delta[i,j]/sqrt(A)) *
        (1 + sqrt(5) * phis * Delta[i,j]/sqrt(A) - 5 * sqrt(5)/3 * phis^3 * Delta[i,j]^3/A^(3/2))
    }
  }
  Sig
}
