#' Spatiotemporal Gaussian Covariance Kernel
#'
#' Computes an adapted non-separable Gaussian covariance kernel
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param lphis log of spatial range
#' @param lphit log of temporal range
#' @param lsig2 log of spatiotemporal variance parameter
#' @keywords st_cov_gaussian
#' @examples
st_cov_gaussian <- function(delta = NULL,
                            Delta = NULL,
                            lphis = NULL,
                            lphit = NULL,
                            lsig2 = NULL){
  Nt = nrow(delta)
  Ns = nrow(Delta)
  N = Ns * Nt

  A = (exp(2 * (lphit + log(delta))) + 1)
  lA = log(A)

  Sig = matrix(NA, nrow = N, ncol = N)
  for(i in 1 : Ns){
    for(j in i : Ns){
      if(Delta[i, j] > 0){
        lD = log(Delta[i, j])

        term.1 = lsig2 - lA
        term.2 = exp(2 * (lphis +  lD) - lA)

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1 - term.2)
      }else{

        term.1 = lsig2 - lA

        Sig[(i - 1) * Nt + 1 : Nt, (j - 1) * Nt + 1 : Nt] =
          Sig[(j - 1) * Nt + 1 : Nt, (i - 1) * Nt + 1 : Nt] = exp(term.1)
      }
    }
  }
  Sig
}
