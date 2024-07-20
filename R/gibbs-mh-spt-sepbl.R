#' A spatiotemporal hierarchical Bayes Markov Chain Monte Carlo sampler
#'
#' Fits a univariate Gaussian spatiotemporal regression model: \eqn{Y(s,t)=x(s,t)^T\beta+Z(s,t)+\epsilon}. Parameters not listed are optional.
#'
#'
#'
#' The posterior being,
#'
#' \eqn{U(\phi_s|a_{\phi_s},b_{\phi_s})\times U(\phi_t|a_{\phi_t},b_{\phi_t}) \times IG(\sigma^2|a_{\sigma},b_{\sigma})\times IG(\tau^2|a_{\tau},b_{\tau})\times N(\beta|\mu_{\beta},\sigma_{\beta}^2)\times N_{N_s}(Z|0,\Sigma)\times\prod_{i_s=1}^{N_s} \prod_{i_t=1}^{N_t} N(y(s_{i_s},t_{i_t})|x(s_{i_s},t_{i_t})^T\beta+Z(s_{i_s},t_{i_t}),\tau^2)}
#'
#' @param t temporal coordinates for observed process (order \eqn{N_t} x \eqn{1})
#' @param coords coordinates for observed process (order \eqn{N_s} x \eqn{2})
#' @param y observed response (order \eqn{N} x  \eqn{1}), \eqn{N= N_s\times N_t}
#' @param X a matrix of covariates (order \eqn{N} x  \eqn{p})
#' @param z_init starting values of spatial effects (order \eqn{N} x  \eqn{1})
#' @param phis_init starting value for \eqn{\phi_s}
#' @param phit_init starting value for \eqn{\phi_t}
#' @param lower_phis lower bound for uniform prior on \eqn{\phi_s}
#' @param upper_phis upper bound for uniform prior on \eqn{\phi_s}
#' @param lower_phit lower bound for uniform prior on \eqn{\phi_t}
#' @param upper_phit upper bound for uniform prior on \eqn{\phi_t}
#' @param sigma2_init starting value for \eqn{\sigma^2}
#' @param shape_sigma shape parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param scale_sigma scale parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param tau2_init starting value for \eqn{\tau^2}
#' @param shape_tau shape parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param scale_tau scale parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param beta_init starting value for \eqn{\beta}
#' @param mean_beta mean parameter for normal prior on \eqn{\beta}
#' @param prec_beta precision (1/variance) parameter for normal prior on \eqn{\beta}
#' @param verbose if true prints output for batches
#' @param steps_init tuning parameter for \eqn{\phi_s}
#' @param stept_init tuning parameter for \eqn{\phi_t}
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param digits number of significant digits
#' @param report batch length
#' @param cov.type.s covariance type (three available choices: Exponential, Gaussian, Mat\'ern(\eqn{\nu=1/2}), Mat\'ern(\eqn{\nu=3/2}), Mat\'ern(\eqn{\nu=5/2})
#' @param cov.type.t covariance type (three available choices: Exponential, Gaussian, Mat\'ern(\eqn{\nu=1/2}), Mat\'ern(\eqn{\nu=3/2}), Mat\'ern(\eqn{\nu=5/2}), Polynomial
#' @keywords hlmBayes_spt
#' @import stats coda
#' @export
##################################################################
### Hierarchical Bayesian spatial model: point referenced data ###
##################################################################
hlmBayes_spt <- function(coords = NULL, t = NULL,
                         y = NULL, X = NULL,
                         z_init = NULL,
                         phis_init = NULL, lower_phis = NULL, upper_phis = NULL,
                         phit_init = NULL, lower_phit = NULL, upper_phit = NULL,
                         sigma2_init = NULL, shape_sigma = NULL, scale_sigma = NULL,
                         tau2_init = NULL, shape_tau = NULL, scale_tau = NULL,
                         beta_init = NULL, mean_beta = NULL, prec_beta = NULL,
                         steps_init = NULL, stept_init = NULL,
                         niter = NULL, nburn = NULL, report = NULL,
                         verbose = FALSE, digits = NULL,
                         cov.type.s = c("exponential", "gaussian", "matern1", "matern2"),
                         cov.type.t = c("exponential", "gaussian", "matern1", "matern2")){
  ##################
  # Chain Defaults #
  ##################
  if(is.null(niter)){
    warning(" Warning: chain length not specified, setting defaults to length = 5e3, burnin = 2.5e3, report = 1e2. ")
    niter = 5e3
    nburn = niter/2
    report = 1e2
  }
  if(is.null(digits)) digits = 3
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
  delta = as.matrix(dist(t))
  Delta = as.matrix(dist(coords))
  if(is.null(X)) p = 1 else p = ncol(X)


  ######################
  # Initialize Storage #
  ######################
  res_phis  =  res_phit  =  res_sig2  =  res_tau2  =  rep(NA, niter)
  if(p == 1){
    res_beta  =  rep(NA, niter)
  }else{
    res_beta  =  matrix(NA, nrow = niter, ncol = p)
  }
  res_z  =  matrix(NA, nrow = niter, ncol = N)

  ##############################
  # Starting Values & Defaults #
  ##############################
  phis  = ifelse(is.null(phis_init), 1, phis_init); lphis = log(phis)
  if(is.null(lower_phis)) lower_phis = 0
  if(is.null(upper_phis)) upper_phis = 30

  phit = ifelse(is.null(phit_init), 1, phit_init); lphit = log(phit)
  if(is.null(lower_phit)) lower_phit = 0
  if(is.null(upper_phit)) upper_phit = 30

  sig2 = ifelse(is.null(sigma2_init), 1, sigma2_init)
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
  if(is.null(mean_beta)) mean_beta = rep(0, p)
  if(is.null(prec_beta)) prec_beta = 1e-6 * diag(p)

  if(is.null(z_init)){
    z = matrix(0, ncol = 1, nrow = N)
  }else{
    z = matrix(z_init, ncol = 1, nrow = N)
  }

  steps_init = ifelse(is.null(steps_init), 1, steps_init)
  stept_init = ifelse(is.null(stept_init), 1, stept_init)

  accepts_vec  =  acceptt_vec  =  c()
  accepts  =  acceptt   =  0
  steps  =  steps_init; stept  =  stept_init


  corr.mat = st_cor_sepbl(delta = delta, Delta = Delta,
                          phis = phis, phit = phit,
                          cov.type.s = cov.type.s,
                          cov.type.t = cov.type.t)

  R = corr.mat$R
  chol.R = corr.mat$R.chol
  R.inv  = corr.mat$R.inv
  E = diag(chol.R)^2

  for(i in 1:niter){

    ############
    # Update Z #
    ############
    Sig.Z.in  =  R.inv/sig2 + diag(N)/tau2
    Sig.Z  =  chol2inv(chol(Sig.Z.in))
    mu.Z  =  crossprod(Sig.Z, as.vector(y))/tau2
    Z  =  crossprod(chol(Sig.Z), rnorm(N)) + mu.Z
    res_beta[i] = beta = mean(Z)
    res_z[i,]  =  z  =  as.vector(t(Z))

    ################
    # Update sigma #
    ################
    post_shape_sigma  =  shape_sigma + N/2
    post_rate_sigma  =  1/scale_sigma + crossprod(t(crossprod(z, R.inv)), z)/2
    res_sig2[i]  =  sig2  =  1/rgamma(1, shape  =  post_shape_sigma, rate  =  post_rate_sigma)

    ##############
    # Update tau #
    ##############
    post_shape_tau  =  shape_tau + N/2
    post_rate_tau  =  1/scale_tau + sum((y - z)^2)/2
    res_tau2[i]  =  tau2  =  1/rgamma(1, shape  =  post_shape_tau, rate  =  post_rate_tau)

    ###########################
    # Metropolis update: phis #
    ###########################
    lphis.draw  =  rnorm(1) * steps + lphis
    if(exp(lphis.draw) < lower_phis | exp(lphis.draw) > upper_phis) accept.prob = -Inf
    else{
      corr.mat.draw_s = st_cor_sepbl(delta = delta, Delta = Delta,
                                     phis = exp(lphis.draw), phit = phit,
                                     cov.type.s = cov.type.s,
                                     cov.type.t = cov.type.t)

      R.draw_s = corr.mat.draw_s$R
      chol.R.draw_s = corr.mat.draw_s$R.chol
      R.inv.draw_s  = corr.mat.draw_s$R.inv
      E.draw_s = diag(chol.R.draw_s)^2


      ra  =  sum(log(sqrt(E/E.draw_s)))
      rb  =   - (lphis.draw - lphis) - crossprod(t(crossprod(z, (R.inv.draw_s - R.inv))), z)/(2 * sig2)
      accept.prob  =  min((ra + rb), 0)
    }


    if(log(runif(1))  <  accept.prob){
      lphis = lphis.draw
      res_phis[i]  =  phis  =  exp(lphis.draw)
      R  =  R.draw_s; R.inv  =  R.inv.draw_s; E  =  E.draw_s
      accepts  =  accepts + 1
    }else{
      res_phis[i]  =  phis
    }

    ###########################
    # Metropolis update: phit #
    ###########################
    lphit.draw  =  rnorm(1) * stept + lphit
    if(exp(lphit.draw) < lower_phit | exp(lphit.draw) > upper_phit) accept.prob = -Inf
    else{
      corr.mat.draw_t = st_cor_sepbl(delta = delta, Delta = Delta,
                                     phis = phis, phit = exp(lphit.draw),
                                     cov.type.s = cov.type.s,
                                     cov.type.t = cov.type.t)

      R.draw_t = corr.mat.draw_t$R
      chol.R.draw_t = corr.mat.draw_t$R.chol
      R.inv.draw_t  = corr.mat.draw_t$R.inv
      E.draw_t = diag(chol.R.draw_t)^2

      ra  =  sum(log(sqrt(E/E.draw_t)))
      rb  =  - (lphit.draw - lphit) - crossprod(t(crossprod(z, (R.inv.draw_t - R.inv))), z)/(2 * sig2)
      accept.prob  =  min((ra + rb), 0)
    }


    if(log(runif(1))  <  accept.prob){
      lphit = lphit.draw
      res_phit[i]  =  phit  =  exp(lphit.draw)
      R  =  R.draw_t; R.inv  =  R.inv.draw_t; E  =  E.draw_t
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
      if(accepts > 0.45) steps  =  steps * accepts/0.45
      else if(accepts < 0.25) steps  =  steps * accepts/0.25


      acceptt  =  acceptt/report
      acceptt  =  max(0.1667, min(acceptt, 0.75))
      acceptt_vec  =  c(acceptt_vec, acceptt)
      if(acceptt > 0.45) stept  =  stept * acceptt/0.45
      else if(acceptt < 0.25) stept  =  stept * acceptt/0.25


      #####################
      # Printing Progress #
      #####################
      if(verbose){
        if(i <= nburn){
          cat("Iteration::", i, "Acceptance:", accepts, acceptt, "\t", "Tuning:", round(c(steps, stept), digits = digits), "\t", "Overall Acceptance Rate:", median(accepts_vec[1:(i/report)]) * 100, "%", median(acceptt_vec[1:(i/report)]) * 100, "%", "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
              "Beta::", round(median(res_beta[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_beta[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_beta[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phis::", round(median(res_phis[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phis[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phit::", round(median(res_phit[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phit[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")"," \n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
        }else{
          cat("Iteration::", i, "Acceptance:", accepts, acceptt, "\t", "Tuning:", round(c(steps, stept), digits),"\t", "Overall Acceptance Rate:", median(accepts_vec[1:(i/report)]) * 100, "%", median(acceptt_vec[1:(i/report)]) * 100, "%", "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
              "Beta::", round(median(res_beta[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_beta[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_beta[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phis::", round(median(res_phis[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phis[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phit::", round(median(res_phit[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phit[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[nburn:i], probs = 0.975), digits = (digits - 1)), ")"," \n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
        }
        accepts  =  acceptt  =  0
      }
    }
  }
  samples.id = (nburn + 1):niter

  return(list(cov.type.s = cov.type.s,
              cov.type.t = cov.type.t,
              y = y, X = X,
              coords = coords, t = t,
              phis = as.mcmc(res_phis[samples.id]),
              phit = as.mcmc(res_phit[samples.id]),
              sig2 = as.mcmc(res_sig2[samples.id]),
              tau2 = as.mcmc(res_tau2[samples.id]),
              beta = as.mcmc(res_beta[samples.id]),
              z = as.mcmc(res_z[samples.id,]),
              acpt_rt = cbind(accepts_vec, acceptt_vec)))
}
