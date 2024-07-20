#' Spatiotemporal Mat\'ern Covariance Kernel with fractal parameter, \eqn{nu = 3/2}
#'
#' Computes an adapted non-separable Mat\'ern covariance kernel with fractal parameter, \eqn{nu = 3/2}
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param lphis log of spatial range
#' @param lphit log of temporal range
#' @param lsig2 log of spatiotemporal variance parameter
#' @keywords st_cov_matern1
#' @examples

st_cov_matern1 <- function(delta = NULL,
                           Delta = NULL,
                           lphis = NULL, # log(phis)
                           lphit = NULL, # log(phit)
                           lsig2 = NULL){ # log(sigma2)
  Nt = nrow(delta)
  Ns = nrow(Delta)
  N = Ns * Nt

  l3by2 = log(3)/2

  A = (exp(2 * (lphit + log(delta))) + 1)
  lA = log(A)

  Sig = matrix(NA, nrow = N, ncol = N)
  for(i in 1 : Ns){
    for(j in i : Ns){
      if(Delta[i, j] > 0){
        lD = log(Delta[i, j])

        term.1 = lsig2 - lA
        term.3 = exp(l3by2 + lphis + lD - lA/2)
        term.2 = log(1 + term.3)

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] = Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1 + term.2 - term.3)
      }else{

        term.1 = lsig2 - lA

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1)
      }
      # exp(log(sig2) - log(A) +
      # log(1 + sqrt(3) * phis * Delta[i, j]/sqrt(A)) - sqrt(3) * phis * Delta[i, j]/sqrt(A))
    }
  }
  Sig
}
