#' A collapsed spatiotemporal hierarchical Bayes Markov Chain Monte Carlo (MCMC) sampler
#'
#' Fits the univariate Gaussian spatiotemporal regression model: \eqn{Y(s,t)=x(s,t)^T\beta+Z(s,t)+\epsilon}. Parameters not listed are optional.
#'
#'
#' The collapsed posterior being,
#'
#' \eqn{U(\phi_s|a_{\phi_s},b_{\phi_s}) \times U(\phi_t|a_{\phi_t},b_{\phi_t}) \times IG(\sigma^2|a_{\sigma},b_{\sigma})\times IG(\tau^2|a_{\tau},b_{\tau}) \times N(\beta|\mu_{\beta},\sigma_{\beta}^2) \times  \mathcal{N}(y(s,t)|x(s_{i_s},t_{i_t})^T\beta,\sigma^2 R(\phi_s,\phi_t)+\tau^2I)}
#'
#' @param y observed response (order \eqn{N} x  \eqn{1}), \eqn{N= N_s\times N_t}
#' @param coords spatial coordinates for observed process (order \eqn{N_s} x \eqn{2})
#' @param t temporal coordinates for observed process (order \eqn{N_t} x \eqn{1})
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param report batch length
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=1/2}), Mat\'ern(\eqn{\nu=3/2}), Mat\'ern(\eqn{\nu=5/2})
#' @param lower.phis lower bound for uniform prior on \eqn{\phi_s}
#' @param upper.phis upper bound for uniform prior on \eqn{\phi_s}
#' @param lower.phit lower bound for uniform prior on \eqn{\phi_t}
#' @param upper.phit upper bound for uniform prior on \eqn{\phi_t}
#' @param a.sigma shape parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param b.sigma scale parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param a.tau shape parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param b.tau scale parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param verbose logical prints output for batches
#' @keywords
#' @import stats coda
#' @export
#' @examples
hlmBayes_mh.spt <- function(y = NULL,
                            coords = NULL,
                            t = NULL,
                            niter = NULL, nburn = NULL, report = NULL,
                            cov.type = NULL,
                            mu.beta = NULL, prec.beta = NULL,
                            a.sigma = NULL, b.sigma = NULL,
                            a.tau = NULL, b.tau = NULL,
                            lower.nu = NULL, upper.nu = NULL,
                            lower.phis = NULL, upper.phis = NULL,
                            lower.phit = NULL, upper.phit = NULL,
                            verbose = TRUE,
                            trgt_fn.compute = FALSE,
                            digits = 3){
  ################################
  # Collapsed MH Sampler for SPT #
  ################################
  delta = as.matrix(dist(t))
  Delta = as.matrix(dist(coords))

  N = length(y)

  if(is.null(niter)){
    niter = 5e3
    nburn = niter/2
    report = 1e2
  }

  # learning rates
  e.sig2 = e.tau2 = e.phis = e.phit = 1e-1
  if(is.null(cov.type)) e.nu = 1e-2

  # Priors
  if(is.null(a.sigma)) a.sigma = 2
  if(is.null(a.tau)) a.tau = 2
  if(is.null(b.sigma)) b.sigma = 1
  if(is.null(b.tau)) b.tau = 0.1
  if(is.null(lower.phis)) lower.phis = 1e-3
  if(is.null(lower.phit)) lower.phit = 1e-3
  if(is.null(upper.phis)) upper.phis = 30
  if(is.null(upper.phit)) upper.phit = 30
  if(is.null(cov.type)){
    if(is.null(lower.nu)) lower.nu = 0.1
    if(is.null(upper.nu)) upper.nu = 1e1
  }

  # initial values
  lphis = lphit = lsig2 = ltau2 =  0
  if(is.null(cov.type)) nu = 1
  else{
    if(cov.type == "exponential") nu = 0.5
    if(cov.type == "matern1") nu = 1.5
    if(cov.type == "matern2") nu = 2.5
    if(cov.type == "gaussian") nu = Inf
  }

  # batch acceptance probabilities
  accept.p = rep(0, 4)
  if(is.null(cov.type)) accept.p = rep(0, 5)

  # storage
  res_sig2 = res_tau2 = res_phis = res_phit =  rep(0, niter)
  accept_m = c()
  if(is.null(cov.type)) res_nu = rep(0, niter)

  if(trgt_fn.compute) trgt_fn  = rep(0, niter)

  # covariance matrices for initial values
  if(is.null(cov.type)){
    R = st_cov_matern(delta = delta, Delta = Delta,
                      lphis = lphis, lphit = lphit,
                      nu = nu,
                      lsig2 = 0)
  }else{
    if(cov.type == "exponential"){
      R = st_cov_exponential(delta = delta, Delta = Delta,
                             lphis = lphis, lphit = lphit,
                             lsig2 = 0)
    }

    if(cov.type == "matern1"){
      R = st_cov_matern1(delta = delta, Delta = Delta,
                         lphis = lphis, lphit = lphit,
                         lsig2 = 0)
    }

    if(cov.type == "matern2"){
      R = st_cov_matern2(delta = delta, Delta = Delta,
                         lphis = lphis, lphit = lphit,
                         lsig2 = 0)
    }

    if(cov.type == "gaussian"){
      R = st_cov_gaussian(delta = delta, Delta = Delta,
                          lphis = lphis, lphit = lphit,
                          lsig2 = 0)
    }
  }

  inv.M.c = chol(exp(lsig2) * R + exp(ltau2) * diag(N))
  inv.M = chol2inv(inv.M.c)

  if(trgt_fn.compute){
    trgt_fn.init = -(a.sigma + 1) * lsig2 - b.sigma/exp(lsig2) -
      (a.tau + 1) * ltau2 - b.tau/exp(ltau2) -
      sum(log(diag(inv.M.c)))/4 - t(y) %*% inv.M %*% (y)/2
  }

  for(i in 1:niter){
    #################
    # Update sigma2 #
    #################
    lsig2.draw = lsig2 + e.sig2 * rnorm(1)

    inv.M.c.draw = chol(exp(lsig2.draw) * R + exp(ltau2) * diag(N))
    inv.M.draw = chol2inv(inv.M.c.draw)

    # with the Jacobian for the log-transform
    ra = - b.sigma * (exp(-lsig2.draw) - exp(-lsig2)) - a.sigma * (lsig2.draw - lsig2) -
      crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
      sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

    # old:: wrong ar
    # ra = - b.sigma * (exp(-lsig2.draw) - exp(-lsig2)) - (a.sigma + 1) * (lsig2.draw - lsig2) -
    #   crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
    #   sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

    accept.prob = min(ra, 0)

    if(log(runif(1)) < accept.prob){
      lsig2 = lsig2.draw
      res_sig2[i] = exp(lsig2)
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[1] = accept.p[1] + 1
    }else{
      res_sig2[i] = exp(lsig2)
    }

    ###############
    # Update tau2 #
    ###############
    ltau2.draw = ltau2 + e.tau2 * rnorm(1)

    inv.M.c.draw = chol(exp(lsig2) * R + exp(ltau2.draw) * diag(N))
    inv.M.draw = chol2inv(inv.M.c.draw)

    ra = - b.tau * (exp(-ltau2.draw) - exp(-ltau2)) - a.tau * (ltau2.draw - ltau2) -
      crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
      sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

    accept.prob = min(ra, 0)

    if(log(runif(1)) < accept.prob){
      ltau2 = ltau2.draw
      res_tau2[i] = exp(ltau2)
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[2] = accept.p[2] + 1
    }else{
      res_tau2[i] = exp(ltau2)
    }

    ###############
    # Update phis #
    ###############
    lphis.draw = lphis + e.phis * rnorm(1)

    if(((exp(lphis.draw) < lower.phis) | (exp(lphis.draw) > upper.phis))){
      accept.prob = -Inf
    }else{
      if(is.null(cov.type)){
        R.draw = st_cov_matern(delta = delta, Delta = Delta,
                               lphis = lphis.draw, lphit = lphit,
                               lsig2 = 0,
                               nu = nu)
      }else{
        if(cov.type == "exponential"){
          R.draw = st_cov_exponential(delta = delta, Delta = Delta,
                                      lphis = lphis.draw, lphit = lphit,
                                      lsig2 = 0)
        }

        if(cov.type == "matern1"){
          R.draw = st_cov_matern1(delta = delta, Delta = Delta,
                                  lphis = lphis.draw, lphit = lphit,
                                  lsig2 = 0)
        }

        if(cov.type == "matern2"){
          R.draw = st_cov_matern2(delta = delta, Delta = Delta,
                                  lphis = lphis.draw, lphit = lphit,
                                  lsig2 = 0)
        }

        if(cov.type == "gaussian"){
          R.draw = st_cov_gaussian(delta = delta, Delta = Delta,
                                   lphis = lphis.draw, lphit = lphit,
                                   lsig2 = 0)
        }
      }

      inv.M.c.draw = chol(exp(lsig2) * R.draw + exp(ltau2) * diag(N))
      inv.M.draw = chol2inv(inv.M.c.draw)

      ra =  - (lphis.draw - lphis) - crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

      # old:: wrong ar
      # ra =  - crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
      #   sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

      accept.prob = min(ra, 0)
    }

    if(log(runif(1)) < accept.prob){
      lphis = lphis.draw
      res_phis[i] = exp(lphis)
      R = R.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[3] = accept.p[3] + 1
    }else{
      res_phis[i] = exp(lphis)
    }

    ###############
    # Update phit #
    ###############
    lphit.draw = lphit + e.phit * rnorm(1)
    if(((exp(lphit.draw) < lower.phit) | (exp(lphit.draw) > upper.phit))){
      accept.prob = -Inf
    }else{
      if(is.null(cov.type)){
        R.draw = st_cov_matern(delta = delta, Delta = Delta,
                               lphis = lphis, lphit = lphit.draw,
                               lsig2 = 0,
                               nu = nu)
      }else{
        if(cov.type == "exponential"){
          R.draw = st_cov_exponential(delta = delta, Delta = Delta,
                                      lphis = lphis, lphit = lphit.draw,
                                      lsig2 = 0)
        }

        if(cov.type == "matern1"){
          R.draw = st_cov_matern1(delta = delta, Delta = Delta,
                                  lphis = lphis, lphit = lphit.draw,
                                  lsig2 = 0)
        }

        if(cov.type == "matern2"){
          R.draw = st_cov_matern2(delta = delta, Delta = Delta,
                                  lphis = lphis, lphit = lphit.draw,
                                  lsig2 = 0)
        }

        if(cov.type == "gaussian"){
          R.draw = st_cov_gaussian(delta = delta, Delta = Delta,
                                   lphis = lphis, lphit = lphit.draw,
                                   lsig2 = 0)
        }
      }

      inv.M.c.draw = chol(exp(lsig2) * R.draw + exp(ltau2) * diag(N))
      inv.M.draw = chol2inv(inv.M.c.draw)

      ra = - (lphit.draw - lphit) - crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

      accept.prob = min(ra, 0)
    }

    if(log(runif(1)) < accept.prob){
      lphit = lphit.draw
      res_phit[i] = exp(lphit)
      R = R.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[4] = accept.p[4] + 1
    }else{
      res_phit[i] = exp(lphit)
    }

    if(is.null(cov.type)){
      #############
      # Update nu #
      #############
      nu.draw = nu + e.nu * rnorm(1)

      R.draw = st_cov_matern(delta = delta, Delta = Delta,
                             lphis = lphis, lphit = lphit,
                             lsig2 = 0,
                             nu = nu.draw)

      inv.M.c.draw = chol(exp(lsig2) * R.draw + exp(ltau2) * diag(N))
      inv.M.draw = chol2inv(inv.M.c.draw)

      ra =  - crossprod(t(crossprod(y, (inv.M.draw - inv.M))), y)/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))

      accept.prob = min(ra + ifelse(test = (nu.draw > lower.nu & nu.draw < upper.nu), 0, -Inf), 0)

      if(log(runif(1)) < accept.prob){
        res_nu[i] = nu = nu.draw
        R = R.draw
        inv.M.c = inv.M.c.draw
        inv.M = inv.M.draw
        accept.p[5] = accept.p[5] + 1
      }else{
        res_nu[i] = nu
      }
    }

    ####################################
    # Target Function Evaluation (Log) #
    ####################################
    if(trgt_fn.compute){
      trgt_fn[i] = - (a.sigma + 1) * lsig2 - b.sigma/exp(lsig2) -
        (a.tau + 1) * ltau2 - b.tau/exp(ltau2) -
        sum(log(diag(inv.M.c)))/4 - t(y) %*% inv.M %*% (y)/2
    }

    ############@@@@@@@@@@@@@@@##############
    # Adaptive Scaling of Proposal Variance #
    #     MH: optimal scales (33%)          #
    ############@@@@@@@@@@@@@@@##############
    if(i %% report == 0){
      accept.p = accept.p/report
      accept_m = rbind(accept_m, accept.p)
      if(verbose){
        if(is.null(cov.type)){
          if(i > nburn){
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.sig2, e.tau2, e.phis, e.phit, e.nu), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phis::", round(median(res_phis[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phis[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phit::", round(median(res_phit[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phit[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[nburn:i], probs = 0.975), digits = (digits - 1)), ")"," \n",
                "Nu::", round(median(res_nu[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_nu[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_nu[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.sig2, e.tau2, e.phis, e.phit, e.nu), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phis::", round(median(res_phis[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phis[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phit::", round(median(res_phit[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phit[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")"," \n",
                "Nu::", round(median(res_nu[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_nu[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_nu[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }
        }else{
          if(i > nburn){
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.sig2, e.tau2, e.phis, e.phit), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phis::", round(median(res_phis[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phis[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phit::", round(median(res_phit[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phit[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[nburn:i], probs = 0.975), digits = (digits - 1)), ")"," \n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.sig2, e.tau2, e.phis, e.phit), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phis::", round(median(res_phis[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phis[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phit::", round(median(res_phit[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phit[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")"," \n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }
        }
      }
      ####################
      # Proposal Scaling #
      ####################
      step = sapply(accept.p, function(x){
        out = 1
        x = max(0.17, min(x, 0.75))
        if(x > 0.50) out = x/0.50
        else if(x < 0.25) out = x/0.25
        out
      })
      e.sig2 = e.sig2 * step[1]
      e.tau2 = e.tau2 * step[2]
      e.phis = e.phis * step[3]
      e.phit = e.phit * step[4]
      if(is.null(cov.type)) e.nu = e.nu * step[5]

      accept.p = rep(0, 4)
      if(is.null(cov.type)) accept.p = rep(0, 5)
    }
  }

  if(is.null(cov.type)){
    if(trgt_fn.compute){
      return(list(y = y,
                  coords = coords,
                  trgt_fn = c(trgt_fn.init, trgt_fn),
                  t = t,
                  nu = res_nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phis = res_phis[(nburn + 1):niter],
                  phit = res_phit[(nburn + 1):niter]))
    }else{
      return(list(y = y,
                  coords = coords,
                  t = t,
                  nu = res_nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phis = res_phis[(nburn + 1):niter],
                  phit = res_phit[(nburn + 1):niter]))
    }
  }else{
    if(trgt_fn.compute){
      return(list(y = y,
                  coords = coords,
                  trgt_fn = c(trgt_fn.init, trgt_fn),
                  t = t,
                  nu = nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phis = res_phis[(nburn + 1):niter],
                  phit = res_phit[(nburn + 1):niter]))
    }else{
      return(list(y = y,
                  coords = coords,
                  t = t,
                  nu = nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phis = res_phis[(nburn + 1):niter],
                  phit = res_phit[(nburn + 1):niter]))
    }
  }
}
