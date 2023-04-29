#' A spatiotemporal hierarchical Bayes Markov Chain Monte Carlo sampler
#'
#' Fits a univariate Gaussian spatiotemporal regression model: \eqn{Y(s,t)=x(s)^T\beta+Z(s)+\epsilon}. Parameters not listed are optional.
#'
#' @param t temporal coordinates for observed process (order \eqn{N_t} x \eqn{1})
#' @param coords coordinates for observed process (order \eqn{N_s} x \eqn{2})
#' @param y observed response (order \eqn{N} x  \eqn{1})
#' @param X a matrix of covariates (order \eqn{N} x  \eqn{p})
#' @param z_init starting values of spatial effects (order \eqn{N} x  \eqn{1})
#' @param D (use if replication at co-ordinate level) index for observations
#' @param phis_init starting value for \eqn{\phi_s}
#' @param phit_init starting value for \eqn{\phi_t}
#' @param lower_phis lower bound for uniform prior on \eqn{\phi_s}
#' @param upper_phis upper bound for uniform prior on \eqn{\phi_s}
#' @param lower_phit lower bound for uniform prior on \eqn{\phi_t}
#' @param upper_phit upper bound for uniform prior on \eqn{\phi_t}
#' @param sigma2_init starting value for \eqn{\phi}
#' @param shape_sigma lower bound for uniform prior on \eqn{\phi}
#' @param scale_sigma upper bound for uniform prior on \eqn{\phi}
#' @param tau2_init starting value for \eqn{\phi}
#' @param shape_tau lower bound for uniform prior on \eqn{\phi}
#' @param scale_tau upper bound for uniform prior on \eqn{\phi}
#' @param beta_init starting value for \eqn{\phi}
#' @param mean_beta lower bound for uniform prior on \eqn{\phi}
#' @param prec_beta upper bound for uniform prior on \eqn{\phi}
#' @param verbose if true prints output for batches
#' @param steps_init tuning parameter for \eqn{\phi_s}
#' @param stept_init tuning parameter for \eqn{\phi_t}
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param report batch length
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @keywords
#' @import stats coda
#' @export
#' @examples
##################################################################
### Hierarchical Bayesian spatial model: point referenced data ###
##################################################################
hlmBayes_spt <- function(coords = NULL, t = NULL,
                         y = NULL, X = NULL,
                         D = NULL,
                         z_init = NULL,
                         phis_init = NULL, lower_phis = NULL, upper_phis = NULL,
                         phit_init = NULL, lower_phit = NULL, upper_phit = NULL,
                         sigma2_init = NULL, shape_sigma = NULL, scale_sigma = NULL,
                         tau2_init = NULL, shape_tau = NULL, scale_tau = NULL,
                         beta_init = NULL, mean_beta = NULL, prec_beta = NULL,
                         steps_init = NULL, stept_init = NULL,
                         niter = NULL, nburn = NULL, report = NULL,
                         verbose = TRUE,
                         cov.type = c("exponential", "gaussian", "matern1", "matern2")){
  ##################
  # Chain Defaults #
  ##################
  if(is.null(niter)){
    warning(" Warning: chain length not specified, setting defaults to length = 1e4, burnin = 5e3, report = 1e2. ")
    niter = 1e4
    nburn = niter/2
    report = 1e2
  }
  if(is.null(nburn) & !is.null(niter)){
    warning(" Warning: burn-in not specified, setting default to niter/2. ")
    nburn = niter/2
  }
  if(is.null(report)){
    warning(" Warning: batch length not specified, setting default to 100. ")
    report = 1e2
  }

  ###################
  # Dimensions Etc. #
  ###################
  N = length(y)
  p = ncol(X)
  delta = as.matrix(dist(t))
  Delta = as.matrix(dist(coords))
  XtX = crossprod(X, X)
  if(is.null(D)){
    D = 1:N
    DtD = diag(N)
  }else{
    DtD = diag(as.vector(table(D)))
  }

  ######################
  # Initialize Storage #
  ######################
  res_phis  =  res_phit  =  res_sigma2  =  res_tau2  =  rep(NA, niter)
  res_beta  =  matrix(NA, nrow = niter, ncol = p)
  res_z  =  matrix(NA, nrow = niter, ncol = N)

  ##############################
  # Starting Values & Defaults #
  ##############################
  phis  = ifelse(is.null(phis_init), 25, phis_init)
  if(is.null(lower_phis)) lower_phis = 0
  if(is.null(upper_phis)) upper_phis = 300

  phit = ifelse(is.null(phit_init), 10, phit_init)
  if(is.null(lower_phit)) lower_phit = 0
  if(is.null(upper_phit)) upper_phit = 300

  sigma2 = ifelse(is.null(sigma2_init), 1, sigma2_init)
  if(is.null(shape_sigma)) shape_sigma = 2
  if(is.null(scale_sigma)) scale_sigma = 1

  tau2 = ifelse(is.null(tau2_init), 1, tau2_init)
  if(is.null(shape_tau)) shape_tau = 2
  if(is.null(scale_tau)) scale_tau = 0.1

  if(p == 1){
    beta  =  ifelse(is.null(beta_init), 0, beta_init)
  }else{
    if(is.null(beta_init)){
      beta = matrix(0, ncol = 1, nrow = p)
    }else{
      beta = matrix(beta_init, ncol = 1, nrow = p)
    }
  }
  if(is.null(mean_beta)) mean_beta = rep(0, ncol(X))
  if(is.null(prec_beta)) prec_beta = 1e-6 * diag(ncol(X))

  if(is.null(z_init)){
    z = matrix(0, ncol = 1, nrow = N)
  }else{
    z = matrix(z_init, ncol = 1, nrow = N)
  }

  accepts_vec  =  acceptt_vec  =  c()
  accepts  =  acceptt   =  0
  steps  =  steps_init; stept  =  stept_init

  if(cov.type == "exponential"){
    R  =  st_cov_exponential(delta = delta, Delta = Delta,
                             phis = phis, phit = phit,
                             sig2 = 1)
    chol.R  =  chol(R)
    R.inv  =  chol2inv(chol.R)
    E  =  diag(chol.R)^2
  }
  if(cov.type == "gaussian"){
    R  =  st_cov_gaussian(delta = delta, Delta = Delta,
                          phis = phis, phit = phit,
                          sig2 = 1)
    chol.R  =  chol(R)
    R.inv  =  chol2inv(chol.R)
    E  =  diag(chol.R)^2
  }
  if(cov.type == "matern1"){
    R  =  st_cov_matern1(delta = delta, Delta = Delta,
                         phis = phis, phit = phit,
                         sig2 = 1)
    chol.R  =  chol(R)
    R.inv  =  chol2inv(chol.R)
    E  =  diag(chol.R)^2
  }
  if(cov.type == "matern2"){
    R  =  st_cov_matern2(delta = delta, Delta = Delta,
                         phis = phis, phit = phit,
                         sig2 = 1)
    chol.R  =  chol(R)
    R.inv  =  chol2inv(chol.R)
    E  =  diag(chol.R)^2
  }

  for(i in 1:niter){
    ###############
    # Update Beta #
    ###############
    post_sd_beta  =  chol2inv(chol(prec_beta + XtX/tau2))
    post_mean_beta  =  post_sd_beta %*% (prec_beta %*% mean_beta + crossprod(X, y - z[D])/tau2)
    res_beta[i,]  =  beta  =  crossprod(chol(post_sd_beta), rnorm(p)) + post_mean_beta
    # as.vector(mvrnorm(1, post_mean_beta, post_sd_beta))

    ############
    # Update Z #
    ############
    Sig.Z.in  =  R.inv/sigma2 + DtD/tau2
    Sig.Z  =  chol2inv(chol(Sig.Z.in))
    mu.Z  =  crossprod(Sig.Z, as.vector(y - crossprod(t(X), beta)))/tau2
    # mu.Z  =  crossprod(Sig.Z,aggregate(as.vector(y - crossprod(t(X),beta)),list(ct$zip),sum)[,2])/tau2
    Z  =  crossprod(chol(Sig.Z), rnorm(N)) + mu.Z
    res_z[i,]  =  z  =  as.vector(t(Z))

    ################
    # Update sigma #
    ################
    post_shape_sigma  =  shape_sigma + N/2
    post_rate_sigma  =  1/scale_sigma + crossprod(t(crossprod(z, R.inv)), z)/2
    res_sigma2[i]  =  sigma2  =  1/rgamma(1, shape  =  post_shape_sigma, rate  =  post_rate_sigma)

    ##############
    # Update tau #
    ##############
    post_shape_tau  =  shape_tau + N/2
    post_rate_tau  =  1/scale_tau + sum((y - crossprod(t(X),beta) - z[D])^2)/2
    res_tau2[i]  =  tau2  =  1/rgamma(1, shape  =  post_shape_tau, rate  =  post_rate_tau)

    ###########################
    # Metropolis update: phis #
    ###########################
    phis.draw  =  runif(1, - 1, 1) * steps + phis
    if(cov.type == "exponential"){
      R.draw_s  =  st_cov_exponential(delta = delta, Delta = Delta,
                                      phis = phis.draw, phit = phit,
                                      sig2 = 1)
    }
    if(cov.type == "gaussian"){
      R.draw_s  =  st_cov_gaussian(delta = delta, Delta = Delta,
                                   phis = phis.draw, phit = phit,
                                   sig2 = 1)
    }
    if(cov.type == "matern1"){
      R.draw_s  =  st_cov_matern1(delta = delta, Delta = Delta,
                                  phis = phis.draw, phit = phit,
                                  sig2 = 1)
    }else{
      R.draw_s  =  st_cov_matern2(delta = delta, Delta = Delta,
                                  phis = phis.draw, phit = phit,
                                  sig2 = 1)
    }

    chol.R.draw_s  =  chol(R.draw_s)
    R.draw.inv_s  =  chol2inv(chol.R.draw_s)
    E.draw_s  =  diag(chol.R.draw_s)^2

    ra  =  sum(log(sqrt(E/E.draw_s)))
    rb  =   - crossprod(t(crossprod(z, (R.draw.inv_s - R.inv))), z)/(2 * sigma2)
    accept.prob  =  min((ra + rb) + ifelse(test  =  (phis.draw > lower_phis & phis.draw < upper_phis), 0, - Inf), 0)

    if(log(runif(1))  <  accept.prob){
      res_phis[i]  =  phis  =  phis.draw
      R  =  R.draw_s; R.inv  =  R.draw.inv_s; E  =  E.draw_s
      accepts  =  accepts + 1
    }else{
      res_phis[i]  =  phis
    }

    ###########################
    # Metropolis update: phit #
    ###########################
    phit.draw  =  runif(1, - 1, 1) * stept + phit
    if(cov.type == "exponential"){
      R.draw_t  =  st_cov_exponential(delta = delta, Delta = Delta,
                                      phis = phis, phit = phit.draw,
                                      sig2 = 1)
    }
    if(cov.type == "gaussian"){
      R.draw_t  =  st_cov_gaussian(delta = delta, Delta = Delta,
                                   phis = phis, phit = phit.draw,
                                   sig2 = 1)
    }
    if(cov.type == "matern1"){
      R.draw_t  =  st_cov_matern1(delta = delta, Delta = Delta,
                                  phis = phis, phit = phit.draw,
                                  sig2 = 1)
    }else{
      R.draw_t  =  st_cov_matern2(delta = delta, Delta = Delta,
                                  phis = phis, phit = phit.draw,
                                  sig2 = 1)
    }

    chol.R.draw_t  =  chol(R.draw_t)
    R.draw.inv_t  =  chol2inv(chol.R.draw_t)
    E.draw_t  =  diag(chol.R.draw_t)^2

    ra  =  sum(log(sqrt(E/E.draw_t)))
    rb  =   - crossprod(t(crossprod(z, (R.draw.inv_t - R.inv))), z)/(2 * sigma2)
    accept.prob  =  min(ra + rb + ifelse(test  =  (phit.draw > lower_phit & phit.draw < upper_phit), 0, - Inf), 0)

    if(log(runif(1))  <  accept.prob){
      res_phit[i]  =  phit  =  phit.draw
      R  =  R.draw_t; R.inv  =  R.draw.inv_t; E  =  E.draw_t
      acceptt  =  acceptt + 1
    }else{
      res_phit[i]  =  phit
    }

    if(i %% report == 0){
      #########################################
      # Adaptive scaling of proposal variance #
      #########################################
      accepts  =  accepts/report
      accepts  =  max(0.1667, min(accepts, 0.75))
      accepts_vec  =  c(accepts_vec, accepts)
      if(accepts > 0.5) steps  =  steps * accepts/0.5
      else if(accepts < 0.25) steps  =  steps * accepts/0.25


      acceptt  =  acceptt/report
      acceptt  =  max(0.1667, min(acceptt, 0.75))
      acceptt_vec  =  c(acceptt_vec, acceptt)
      if(acceptt > 0.5) stept  =  stept * acceptt/0.5
      else if(acceptt < 0.25) stept  =  stept * acceptt/0.25


      #####################
      # Printing Progress #
      #####################
      if(verbose){
        cat("Iteration ", i, "\n")
        if(i <= nburn){
          cat("-------------------------","\n",
              "phi.s:", "\t", round(median(res_phis[1:i]), 3), "\t",
              "phi.t:", "\t", round(median(res_phit[1:i]), 3), "\n",
              "sigma.2: ", "\t", round(median(res_sigma2[1:i]), 3), "\t",
              "tau.2:", "\t", round(median(res_tau2[1:i]), 3), "\n",
              "Acceptance Rate (phi.s):","\t", accepts * 100, "%", "\t",
              ", (phi.t):", "\t", acceptt * 100, "%", "\n",
              "Overall Acceptance Rate:", median(accepts_vec[1:(i/report)]) * 100, "%", "\t", median(acceptt_vec[1:(i/report)]) * 100,"%", "\n",
              "-------------------------", "\n")
        }else{
          cat("-------------------------","\n",
              "phi.s:", "\t", round(median(res_phis[nburn:i]), 3), "\t",
              "phi.t:", "\t", round(median(res_phit[nburn:i]), 3), "\n",
              "sigma.2:", "\t", round(median(res_sigma2[nburn:i]), 3), "\t",
              "tau.2:", "\t", round(median(res_tau2[nburn:i]), 3), "\n",
              "Acceptance Rate (phi.s):","\t", accepts * 100, "%", "\t",
              ", (phi.t):","\t", acceptt * 100, "%", "\n",
              "Overall Acceptance Rate:", median(accepts_vec[(nburn/report):(i/report)]) * 100, "%", "\t", median(acceptt_vec[(nburn/report):(i/report)]) * 100, "%", "\n",
              "-------------------------", "\n")
        }
        accepts  =  acceptt  =  0
      }
    }
  }
  return(list(post_phis = as.mcmc(res_phis),
              post_phit = as.mcmc(res_phit),
              post_sigma2 = as.mcmc(res_sigma2),
              post_tau2 = as.mcmc(res_tau2),
              post_beta = as.mcmc(res_beta),
              post_z = as.mcmc(res_z),
              ar_s = accepts_vec,
              ar_t = acceptt_vec))
}
