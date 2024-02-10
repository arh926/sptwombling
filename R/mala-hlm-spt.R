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
#' @param lower.phis lower bound for uniform prior on \eqn{\phi_s}
#' @param upper.phis upper bound for uniform prior on \eqn{\phi_s}
#' @param lower.phit lower bound for uniform prior on \eqn{\phi_t}
#' @param upper.phit upper bound for uniform prior on \eqn{\phi_t}
#' @param a.sigma shape parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param b.sigma scale parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param a.tau shape parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param b.tau scale parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param mu.beta mean parameter for normal prior on \eqn{\beta}
#' @param prec.beta precision (1/variance) parameter for normal prior on \eqn{\beta}
#' @param verbose if true prints output for batches
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param report batch length
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=1/2}), Mat\'ern(\eqn{\nu=3/2}), Mat\'ern(\eqn{\nu=5/2})
#' @keywords
#' @import stats coda MASS
#' @export
#' @examples
hlmBayes_mala.spt <- function(y = NULL,
                              coords = NULL,
                              t = NULL,
                              X = NULL,
                              niter = NULL, nburn = NULL, report = NULL,
                              cov.type = c("exponential", "matern1", "matern2" ,"gaussian"),
                              digits = 3,
                              verbose = TRUE,
                              track = TRUE,
                              mu.beta = NULL, prec.beta = NULL,
                              a.sigma = NULL, b.sigma = NULL,
                              a.tau = NULL, b.tau = NULL,
                              lower.phis = NULL, upper.phis = NULL,
                              lower.phit = NULL, upper.phit = NULL){

  ##################################
  # Collapsed MALA Sampler for SPT #
  ##################################
  delta = as.matrix(dist(t))
  Delta = as.matrix(dist(coords))

  N = length(y)
  if(is.null(X)) X = matrix(1, nrow = N, ncol = 1)
  p = ncol(X)

  if(is.null(niter)){
    niter = 5e4
    nburn = niter/2
    report = 1e2
  }

  # learning rates
  e.beta = e.sig2 = e.phis = e.phit = e.tau2 = 0.1

  # Priors
  if(is.null(a.sigma)) a.sigma = 2
  if(is.null(a.tau)) a.tau = 2
  if(is.null(b.sigma)) b.sigma = 1
  if(is.null(b.tau)) b.tau = 1
  if(is.null(mu.beta)) mu.beta = rep(0, p)
  if(is.null(prec.beta)) prec.beta = chol2inv(chol(1e6 * diag(p)))
  if(is.null(lower.phis)) lower.phis = 0
  if(is.null(lower.phit)) lower.phit = 0
  if(is.null(upper.phis)) upper.phis = 300
  if(is.null(upper.phit)) upper.phit = 300


  phis.init = phis = 1
  phit.init = phit = 1
  sig2.init = sig2 = 1
  tau2.init = tau2 = 1
  beta.init = beta = rep(0, p)

  accept.p = rep(0, 4) # batch acceptance probabilities

  # storage
  res_beta = matrix(0, nrow = niter, ncol = p)
  res_sig2 = res_tau2 = res_phis = res_phit = trgt_fn  = rep(0, niter)
  accept_m = c()

  if(cov.type == "exponential"){
    R = st_cov_exponential(delta = delta, Delta = Delta,
                           phis = phis, phit = phit,
                           sig2 = 1)
  }

  if(cov.type == "matern1"){
    R = st_cov_matern1(delta = delta, Delta = Delta,
                           phis = phis, phit = phit,
                           sig2 = 1)
  }

  if(cov.type == "matern2"){
    R = st_cov_matern2(delta = delta, Delta = Delta,
                       phis = phis, phit = phit,
                       sig2 = 1)
  }

  if(cov.type == "gaussian"){
    R = st_cov_gaussian(delta = delta, Delta = Delta,
                       phis = phis, phit = phit,
                       sig2 = 1)
  }
  inv.M.c = chol(sig2 * R + tau2 * diag(N))
  inv.M = chol2inv(inv.M.c)


  trgt_fn.init = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 - -(a.tau + 1) * log(tau2) - b.tau/tau2 - log(upper.phis - lower.phis) - log(upper.phit - lower.phit) -
    t(beta - mu.beta) %*% prec.beta %*% (beta - mu.beta)/2 + sum(log(diag(inv.M.c))) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta)/2

  for(i in 1:niter){
    ###############
    # Update beta #
    ###############
    Sigma.beta = prec.beta + t(X) %*% inv.M %*% X
    Sigma.beta.inv = chol2inv(chol(Sigma.beta))
    Mu.beta = Sigma.beta.inv %*% (prec.beta %*% mu.beta + t(X) %*% inv.M %*% y)
    res_beta[i,] = beta = as.vector(mvrnorm(1, Mu.beta, Sigma.beta))

    #################
    # Update sigma2 #
    #################
    nabla.sig2 =  - (a.sigma + 1)/sig2 + b.sigma/sig2^2 - sum(diag(inv.M %*% R))/2 + t(y - X %*% beta) %*% inv.M %*% R %*% inv.M %*% (y - X %*% beta)/2

    sig2.draw = as.numeric(sig2 + e.sig2^2 * nabla.sig2/2 + e.sig2 * rnorm(1))
    if(sig2.draw < 0){
      accept.prob = -Inf
    }else{
      inv.M.c.draw = chol(sig2.draw * R + tau2 * diag(N))
      inv.M.draw = chol2inv(inv.M.c.draw)
      nabla.sig2.draw = - (a.sigma + 1)/sig2.draw + b.sigma/sig2.draw^2 - sum(diag(inv.M.draw %*% R))/2 + t(y - X %*% beta) %*% inv.M.draw %*% R %*% inv.M.draw %*% (y - X %*% beta)/2

      ra = - b.sigma * (1/sig2.draw - 1/sig2) - (a.sigma + 1) * log(sig2.draw/sig2) -
        (t(y - X %*% beta) %*% inv.M.draw %*% (y - X %*% beta) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta))/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))
      rb = (sig2.draw - sig2 - e.sig2^2 * nabla.sig2/2)^2/2/e.sig2^2 - (sig2 - sig2.draw - e.sig2^2 * nabla.sig2.draw/2)^2/2/e.sig2^2

      accept.prob = min(ra + rb, 0)
    }

    if(log(runif(1)) < accept.prob){
      res_sig2[i] = sig2 = sig2.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[1] = accept.p[1] + 1
    }else{
      res_sig2[i] = sig2
    }

    ###############
    # Update tau2 #
    ###############
    nabla.tau2 =  - (a.tau + 1)/tau2 + b.tau/tau2^2 - sum(diag(inv.M))/2 + t(y - X %*% beta) %*% inv.M %*% inv.M %*% (y - X %*% beta)/2

    tau2.draw = as.numeric(tau2 + e.tau2^2 * nabla.tau2/2 + e.tau2 * rnorm(1))
    if(tau2.draw < 0){
      accept.prob = -Inf
    }else{
      inv.M.c.draw = chol(sig2 * R + tau2.draw * diag(N))
      inv.M.draw = chol2inv(inv.M.c.draw)
      nabla.tau2.draw = - (a.tau + 1)/tau2.draw + b.tau/tau2.draw^2 - sum(diag(inv.M.draw))/2 + t(y - X %*% beta) %*% inv.M.draw %*% inv.M.draw %*% (y - X %*% beta)/2

      ra = - b.tau * (1/tau2.draw - 1/tau2) - (a.tau + 1) * log(tau2.draw/tau2) -
        (t(y - X %*% beta) %*% inv.M.draw %*% (y - X %*% beta) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta))/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))
      rb = (tau2.draw - tau2 - e.tau2^2 * nabla.tau2/2)^2/2/e.tau2^2 - (tau2 - tau2.draw - e.tau2^2 * nabla.tau2.draw/2)^2/2/e.tau2^2

      accept.prob = min(ra + rb, 0)
    }

    if(log(runif(1)) < accept.prob){
      res_tau2[i] = tau2 = tau2.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[2] = accept.p[2] + 1
    }else{
      res_tau2[i] = tau2
    }

    ###############
    # Update phis #
    ###############
    if(cov.type == "exponential"){
      del.R.del.phis = del.phis.st_cov_exponential(delta = delta, Delta = Delta,
                                                   phis = phis, phit = phit,
                                                   sig2 = sig2)
    }
    if(cov.type == "matern1"){
      del.R.del.phis = del.phis.st_cov_matern1(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)
    }
    if(cov.type == "matern2"){
      del.R.del.phis = del.phis.st_cov_matern2(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)
    }
    if(cov.type == "gaussian"){
      del.R.del.phis = del.phis.st_cov_gaussian(delta = delta, Delta = Delta,
                                                phis = phis, phit = phit,
                                                sig2 = sig2)
    }

    nabla.phis = -sum(diag(inv.M %*% del.R.del.phis))/2 + t(y - X %*% beta) %*% inv.M %*% del.R.del.phis %*% inv.M %*% (y - X %*% beta)/2

    phis.draw = as.numeric(phis + e.phis^2 * nabla.phis/2 + e.phis * rnorm(1))

    if(cov.type == "exponential"){
      del.R.del.phis.draw = del.phis.st_cov_exponential(delta = delta, Delta = Delta,
                                                        phis = phis.draw, phit = phit,
                                                        sig2 = sig2)
    }
    if(cov.type == "matern1"){
      del.R.del.phis.draw = del.phis.st_cov_matern1(delta = delta, Delta = Delta,
                                                    phis = phis.draw, phit = phit,
                                                    sig2 = sig2)
    }
    if(cov.type == "matern2"){
      del.R.del.phis.draw = del.phis.st_cov_matern2(delta = delta, Delta = Delta,
                                                    phis = phis.draw, phit = phit,
                                                    sig2 = sig2)
    }
    if(cov.type == "gaussian"){
      del.R.del.phis.draw = del.phis.st_cov_gaussian(delta = delta, Delta = Delta,
                                                     phis = phis.draw, phit = phit,
                                                     sig2 = sig2)
    }

    if(cov.type == "exponential"){
      R.draw = st_cov_exponential(delta = delta, Delta = Delta,
                                  phis = phis.draw, phit = phit,
                                  sig2 = 1)
    }
    if(cov.type == "matern1"){
      R.draw = st_cov_matern1(delta = delta, Delta = Delta,
                              phis = phis.draw, phit = phit,
                              sig2 = 1)
    }

    if(cov.type == "matern2"){
      R.draw = st_cov_matern2(delta = delta, Delta = Delta,
                              phis = phis.draw, phit = phit,
                              sig2 = 1)
    }

    if(cov.type == "gaussian"){
      R.draw = st_cov_gaussian(delta = delta, Delta = Delta,
                               phis = phis.draw, phit = phit,
                               sig2 = 1)
    }
    inv.M.c.draw = chol(sig2 * R.draw + tau2 * diag(N))
    inv.M.draw = chol2inv(inv.M.c.draw)

    nabla.phis.draw = -sum(diag(inv.M.draw %*% del.R.del.phis.draw))/2 + t(y - X %*% beta) %*% inv.M.draw %*% del.R.del.phis.draw %*% inv.M.draw %*% (y - X %*% beta)/2

    ra =  - (t(y - X %*% beta) %*% inv.M.draw %*% (y - X %*% beta) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta))/2 -
      sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))
    rb = (phis.draw - phis - e.phis^2 * nabla.phis/2)^2/2/e.phis^2 - (phis - phis.draw - e.phis^2 * nabla.phis.draw/2)^2/2/e.phis^2

    accept.prob = min(ra + rb + ifelse(test = (phis.draw > lower.phis & phis.draw < upper.phis), 0, -Inf), 0)
    if(log(runif(1)) < accept.prob){
      res_phis[i] = phis = phis.draw
      R = R.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[3] = accept.p[3] + 1
    }else{
      res_phis[i] = phis
    }

    ###############
    # Update phit #
    ###############
    if(cov.type == "exponential"){
      del.R.del.phit = del.phit.st_cov_exponential(delta = delta, Delta = Delta,
                                                   phis = phis, phit = phit,
                                                   sig2 = sig2)
    }
    if(cov.type == "matern1"){
      del.R.del.phit = del.phit.st_cov_matern1(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)
    }
    if(cov.type == "matern2"){
      del.R.del.phit = del.phit.st_cov_matern2(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)
    }
    if(cov.type == "gaussian"){
      del.R.del.phit = del.phit.st_cov_gaussian(delta = delta, Delta = Delta,
                                                phis = phis, phit = phit,
                                                sig2 = sig2)
    }
    nabla.phit = - sum(diag(inv.M %*% del.R.del.phit))/2 + t(y - X %*% beta) %*% inv.M %*% del.R.del.phit %*% inv.M %*% (y - X %*% beta)/2

    phit.draw = as.numeric(phit + e.phit^2 * nabla.phit/2 + e.phit * rnorm(1))
    if(cov.type == "exponential"){
      del.R.del.phit.draw = del.phit.st_cov_exponential(delta = delta, Delta = Delta,
                                                        phis = phis, phit = phit.draw,
                                                        sig2 = sig2)
    }
    if(cov.type == "matern1"){
      del.R.del.phit.draw = del.phit.st_cov_matern1(delta = delta, Delta = Delta,
                                                    phis = phis, phit = phit.draw,
                                                    sig2 = sig2)
    }
    if(cov.type == "matern2"){
      del.R.del.phit.draw = del.phit.st_cov_matern2(delta = delta, Delta = Delta,
                                                    phis = phis, phit = phit.draw,
                                                    sig2 = sig2)
    }
    if(cov.type == "gaussian"){
      del.R.del.phit.draw = del.phit.st_cov_gaussian(delta = delta, Delta = Delta,
                                                     phis = phis, phit = phit.draw,
                                                     sig2 = sig2)
    }

    if(cov.type == "exponential"){
      R.draw = st_cov_exponential(delta = delta, Delta = Delta,
                                  phis = phis, phit = phit.draw,
                                  sig2 = 1)
    }
    if(cov.type == "matern1"){
      R.draw = st_cov_matern1(delta = delta, Delta = Delta,
                              phis = phis, phit = phit.draw,
                              sig2 = 1)
    }

    if(cov.type == "matern2"){
      R.draw = st_cov_matern2(delta = delta, Delta = Delta,
                              phis = phis, phit = phit.draw,
                              sig2 = 1)
    }

    if(cov.type == "gaussian"){
      R.draw = st_cov_gaussian(delta = delta, Delta = Delta,
                               phis = phis, phit = phit.draw,
                               sig2 = 1)
    }
    inv.M.c.draw = chol(sig2 * R.draw + tau2 * diag(N))
    inv.M.draw = chol2inv(inv.M.c.draw)
    nabla.phit.draw = -sum(diag(inv.M.draw %*% del.R.del.phit.draw))/2 + t(y - X %*% beta) %*% inv.M.draw %*% del.R.del.phit.draw %*% inv.M.draw %*% (y - X %*% beta)/2

    ra =  - (t(y - X %*% beta) %*% inv.M.draw %*% (y - X %*% beta) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta))/2 -
      sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))
    rb = (phit.draw - phit - e.phit^2 * nabla.phit/2)^2/2/e.phit^2 - (phit - phit.draw - e.phit^2 * nabla.phit.draw/2)^2/2/e.phit^2

    accept.prob = min(ra + rb + ifelse(test = (phit.draw > lower.phit & phit.draw < upper.phit), 0, -Inf), 0)
    if(log(runif(1)) < accept.prob){
      res_phit[i] = phit = phit.draw
      R = R.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p[4] = accept.p[4] + 1
    }else{
      res_phit[i] = phit
    }

    ####################################
    # Target Function Evaluation (Log) #
    ####################################
    trgt_fn[i] = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 - -(a.tau + 1) * log(tau2) - b.tau/tau2 - log(upper.phis - lower.phis) - log(upper.phit - lower.phit) -
      t(beta - mu.beta) %*% prec.beta %*% (beta - mu.beta)/2 + sum(log(diag(inv.M.c))) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta)/2

    ############@@@@@@@@@@@@@@@##############
    # Adaptive Scaling of Proposal Variance #
    #     MALA: optimal scales (58%)        #
    ############@@@@@@@@@@@@@@@##############
    if(i %% report == 0){
      accept.p = accept.p/report
      accept_m = rbind(accept_m, accept.p)
      if(verbose){
        if(i >= nburn){
          cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.beta, e.sig2, e.tau2, e.phis, e.phit), digits), "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
              # "Beta::", round(apply(res_beta[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(e.beta, (report %% 10)),"\n",
              "Beta::", round(median(res_beta[nburn:i,]), digits = (digits - 1)), "(", round(quantile(res_beta[nburn:i,], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_beta[nburn:i,], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phis::", round(median(res_phis[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phis[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Phit::", round(median(res_phit[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phit[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[nburn:i], probs = 0.975), digits = (digits - 1)), ")"," \n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
        }else{
          if(track){
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.beta, e.sig2, e.tau2, e.phis, e.phit), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                # "Beta::", round(apply(res_beta[(i - report + 1):i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(e.beta, (report %% 10)),"\n",
                "Beta::", round(median(res_beta[(i - report + 1):i,]), digits = (digits - 1)), "(", round(quantile(res_beta[(i - report + 1):i,], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_beta[(i - report + 1):i,], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phis::", round(median(res_phis[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phis[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phis[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phit::", round(median(res_phit[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phit[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phit[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")"," \n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "\n")
          }
        }
      }
      ####################
      # Proposal Scaling #
      ####################
      step = sapply(accept.p, function(x){
        out = 1
        x = max(0.25, min(x, 0.75))
        if(x > 0.65) out = x/0.65
        else if(x < 0.45) out = x/0.45
        out
      })
      e.sig2 = e.sig2 * step[1]
      e.tau2 = e.tau2 * step[2]
      e.phis = e.phis * step[3]
      e.phit = e.phit * step[4]

      accept.p = rep(0, 4)
    }
  }

  return(list(y = y,
              coords = coords,
              t = t,
              trgt_fn = c(trgt_fn.init,trgt_fn),
              beta = res_beta,
              sig2 = res_sig2,
              tau2 = res_tau2,
              phis = res_phis,
              phit = res_phit))

}


# Test Examples::
# spt_model0 = hlmBayes_mala.spt(y = y, coords = coords, t = t, cov.type = "exponential")
# spt_model1 = hlmBayes_mala.spt(y = y, coords = coords, t = t, cov.type = "matern1")
# spt_model2 = hlmBayes_mala.spt(y = y, coords = coords, t = t, cov.type = "matern2")
# spt_model3 = hlmBayes_mala.spt(y = y, coords = coords, t = t, cov.type = "gaussian")


# save(mcmc.list, file = "mala-sqnc-mcmc.RData")
#
#
# plot(c(trgt_fn.init,trgt_fn), type="l", col = "darkred", ylab = "Target Probability")
# grid()
# par(mfrow=c(5,3))
# plot_mcmc(samples = res_beta[nburn:niter,], cnames = "beta0")
# plot_mcmc(samples = res_sig2[nburn:niter], cnames = "sigma2")
# plot_mcmc(samples = res_tau2[nburn:niter], cnames = "tau2")
# plot_mcmc(samples = res_phis[nburn:niter], cnames = "phis")
# plot_mcmc(samples = res_phit[nburn:niter], cnames = "phit")
