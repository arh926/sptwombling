#' Spatiotemporal Mat\'ern Covariance Kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' Computes an adapted non-separable Mat\'ern covariance kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param lphis log of spatial range
#' @param lphit log of temporal range
#' @param lsig2 log of spatiotemporal variance parameter
#' @keywords st_cov_matern2
st_cov_matern2 <- function(delta = NULL,
                           Delta = NULL,
                           lphis = NULL,
                           lphit = NULL,
                           lsig2 = NULL){
  Nt = nrow(delta)
  Ns = nrow(Delta)
  N = Ns * Nt

  l5by2 = log(5)/2
  l5 = 2 * l5by2
  l3 = log(3)

  A = (exp(2 * (lphit + log(delta))) + 1)
  lA = log(A)

  Sig = matrix(NA, nrow = N, ncol = N)
  for(i in 1 : Ns){
    for(j in i : Ns){
      if(Delta[i, j] > 0){
        lD = log(Delta[i, j])

        term.1 = lsig2 - lA
        term.3 = exp(l5by2 + lphis + lD - lA/2)
        term.2 = log(1 + term.3 + exp(l5 + 2 * lphis + 2 * lD - l3 - lA))

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1 + term.2 - term.3)
      }else{
        term.1 = lsig2 - lA

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1)
      }
      # sig2/A *
      # (1 + sqrt(5) * phis * Delta[i, j]/sqrt(A) + 5 * phis^2 * Delta[i,j]^2/(3 * A)) *
      # exp( - sqrt(5) * phis * Delta[i, j]/sqrt(A))
    }
  }
  Sig
}
