#' Spatiotemporal Exponential Covariance Kernel
#'
#' Computes a separable exponential covariance kernel
#'
#' @param delta temporal distance matrix
#' @param Delta spatial distance matrix
#' @param phis spatial range
#' @param phit temporal range
#' @param sig2 variance parameter
#' @param cov.type.s spatial covariance type (four available choices: Exponential, Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @param cov.type.t spatial covariance type (five available choices: Exponential, Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2}, Inverse/polynomial)
#' @keywords st_cor_sepbl
#' @examples
st_cor_sepbl <- function(delta = NULL,
                         Delta = NULL,
                         phis = NULL,
                         phit = NULL,
                         cov.type.s = c("exponential", "gaussian", "matern1", "matern2"),
                         cov.type.t = c("exponential", "gaussian", "matern1", "matern2", "poly")){

  if(cov.type.t == "exponential"){
    R.t = exp(-phit * delta)
  }else if(cov.type.t == "matern1"){
    R.t = (1 + sqrt(3) * phit * delta) * exp(-sqrt(3) * phit * delta)
  }else if(cov.type.t == "matern2"){
    R.t = (1 + sqrt(5) * phit * delta + 5/3 * phit^2 * delta^2) * exp(-sqrt(5) * phit * delta)
  }else if(cov.type.t == "gaussian"){
    R.t = exp(-phit^2 * delta^2)
  }else if(cov.type.t == "poly"){
    R.t = 1/(phit^2 * delta^2 + 1)
  }

  chol.r.t = chol(R.t + 1e-4 * diag(nrow(delta)))
  inv.r.t = chol2inv(chol.r.t)

  if(cov.type.s == "exponential"){
    R.s = exp(-phis * Delta)
  }else if(cov.type.s == "matern1"){
    R.s = (1 + sqrt(3) * phis * Delta) * exp(-sqrt(3) * phis * Delta)
  }else if(cov.type.s == "matern2"){
    R.s = (1 + sqrt(5) * phis * Delta + 5/3 * phis^2 * Delta^2) * exp(-sqrt(5) * phis * Delta)
  }else if(cov.type.s == "gaussian"){
    R.s = exp(-phis^2 * Delta^2)
  }

  chol.r.s = chol(R.s + 1e-4 * diag(nrow(Delta)))
  inv.r.s = chol2inv(chol.r.s)

  R = kronecker(R.s, R.t)
  R.chol  = kronecker(chol.r.s, chol.r.t)
  R.inv = kronecker(inv.r.s, inv.r.t)

  return(list(R.t = R.t,
              R.s = R.s,
              R = R,
              R.chol = R.chol,
              R.inv = R.inv))
}
