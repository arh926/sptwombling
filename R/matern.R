#' Spatiotemporal Mat\'ern Covariance Kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' Computes an adapted non-separable Mat\'ern covariance kernel with fractal parameter, \eqn{nu = 5/2}
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param lphis log spatial range
#' @param lphit log temporal range
#' @param lsig2 log spatiotemporal variance parameter
#' @param nu fractal parameter
#' @keywords st_cov_matern
st_cov_matern <- function(delta = NULL,
                          Delta = NULL,
                          lphis = NULL,
                          lphit = NULL,
                          lsig2 = NULL,
                          nu = NULL){
  Nt = nrow(delta)
  Ns = nrow(Delta)
  N = Ns * Nt

  A = (exp(2 * (lphit + log(delta))) + 1)
  lA = log(A)

  lnuby2 = log(nu)/2
  l2 = log(2)

  Sig = matrix(NA, nrow = N, ncol = N)
  for(i in 1 : Ns){
    for(j in i : Ns){
      if(Delta[i, j] > 0){

        lD = log(Delta[i, j])

        term.1 = lsig2 - lA
        term.31 = lnuby2 + l2/2 + lphis + lD - lA/2
        term.2 = (1 - nu) * l2 - lfactorial(nu - 1) + nu * (term.31)
        term.3 = log(besselK(as.matrix(term.31), nu = nu, expon.scaled = TRUE)) - as.matrix(term.31)

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1 + term.2 + term.3)

        # sig2 * 2^(1 - nu) / gamma(nu) / A *
        #   (sqrt(2 * nu) * phis * Delta[i, j]/sqrt(A))^nu * besselK(as.matrix(sqrt(2 * nu) * phis * Delta[i, j]/sqrt(A)),
        #                                                            nu = nu, expon.scaled = TRUE)/exp(as.matrix(sqrt(2 * nu) * phis * Delta[i, j]/sqrt(A)))
      }else{

        term.1 = lsig2 - lA

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1)
      }
    }
  }
  Sig
}
