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
hlmBayes_mala_pc.spt <- function(y = NULL,
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
  ################################## ##############
  # Collapsed Preconditioned MALA Sampler for SPT #
  #################################################
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
  e.theta = 0.1

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

  accept.p = 0 # batch acceptance probabilities

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
    Sigma.beta = prec.beta + t(X) %*% inv.M %*% X                                  # posterior variance
    Sigma.beta.inv = chol2inv(chol(Sigma.beta))                                    # inverse of posterior variance
    Mu.beta = Sigma.beta.inv %*% (prec.beta %*% mu.beta + t(X) %*% inv.M %*% y)    # posterior mean
    res_beta[i,] = beta = as.vector(mvrnorm(1, Mu.beta, Sigma.beta))               # Gibbs step

    #####################@@@@@@@@@@@@@#############################
    # Update all parameters in one step: (sig2, phis, phit, tau2) #
    #####################@@@@@@@@@@@@@#############################
    theta = matrix(c(sig2, tau2, phis, phit), nc = 1)                              # current parameter value
    if(cov.type == "exponential"){
      del.R.del.phis = del.phis.st_cov_exponential(delta = delta, Delta = Delta,
                                                   phis = phis, phit = phit,
                                                   sig2 = sig2)                          # gradient of R w.r.t phis
      del.R.del.phit = del.phit.st_cov_exponential(delta = delta, Delta = Delta,
                                                   phis = phis, phit = phit,
                                                   sig2 = sig2)                          # gradient of R w.r.t phit
    }
    if(cov.type == "matern1"){

      del.R.del.phis = del.phis.st_cov_matern1(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)                          # gradient of R w.r.t phis
      del.R.del.phit = del.phit.st_cov_matern1(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)                          # gradient of R w.r.t phit
    }
    if(cov.type == "matern2"){

      del.R.del.phis = del.phis.st_cov_matern2(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)                          # gradient of R w.r.t phis
      del.R.del.phit = del.phit.st_cov_matern2(delta = delta, Delta = Delta,
                                               phis = phis, phit = phit,
                                               sig2 = sig2)                          # gradient of R w.r.t phit
    }
    if(cov.type == "gaussian"){

      del.R.del.phis = del.phis.st_cov_gaussian(delta = delta, Delta = Delta,
                                                phis = phis, phit = phit,
                                                sig2 = sig2)                          # gradient of R w.r.t phis
      del.R.del.phit = del.phit.st_cov_gaussian(delta = delta, Delta = Delta,
                                                phis = phis, phit = phit,
                                                sig2 = sig2)                          # gradient of R w.r.t phit
    }
    # terms that are required multiple times
    t1 = inv.M %*% R
    t2 = inv.M %*% inv.M
    t3 = inv.M %*% del.R.del.phis
    t4 = inv.M %*% del.R.del.phit
    t5 = inv.M %*% (y - X %*% beta)

    ###################
    # Gradient vector #
    ###################
    g.sig2 = - (a.sigma + 1)/sig2 + b.sigma/sig2^2 - sum(diag(t1))/2 + t(y - X %*% beta) %*% t1 %*% t5/2 # sigma2
    g.tau2 = - (a.tau + 1)/tau2 + b.tau/tau2^2 - sum(diag(inv.M))/2 + t(y - X %*% beta) %*% t2 %*% (y - X %*% beta)/2 # tau2
    g.phis = - sum(diag(t3))/2 + t(y - X %*% beta) %*% t3 %*% t5/2 # phis
    g.phit = - sum(diag(t4))/2 + t(y - X %*% beta) %*% t4 %*% t5/2 # phit
    g.theta = matrix(c(g.sig2, g.tau2, g.phis, g.phit), nc = 1) # gradient vector

    ########################
    # Fisher's Information #
    ########################
    i.11 =  -(a.sigma + 1)/sig2^2 + 2 * b.sigma/sig2^3 + sum(diag(t1 %*% t1))/2    # sigma2, sigma2
    i.12 = sum(diag(t2 %*% R))/2                                                   # sigma2, tau2
    i.13 = sum(diag(t3 %*% t1))/2                                                  # sigma2, phis
    i.14 = sum(diag(t4 %*% t1))/2                                                  # sigma2, phit
    i.22 = -(a.tau + 1)/tau2^2 + 2 * b.tau/tau2^3 + sum(diag(t2))/2                # tau2, tau2
    i.23 = sum(diag(inv.M %*% t3))/2                                               # tau2, phis
    i.24 = sum(diag(inv.M %*% t4))/2                                               # tau2, phit
    i.33 = sum(diag(t3 %*% t3))/2                                                  # phis, phis
    i.34 = sum(diag(t4 %*% t3))/2                                                  # phis, phit
    i.44 = sum(diag(t4 %*% t4))/2                                                  # phit, phit
    i.theta = matrix(c(i.11, i.12, i.13, i.14,
                       i.12, i.22, i.23, i.24,
                       i.13, i.23, i.33, i.34,
                       i.14, i.24, i.34, i.44),
                     ncol = 4, nrow = 4, byrow = T)                                # Fisher's Information
    i.theta.inv = try(chol2inv(chol(i.theta)), silent = TRUE)                      # inverse of Fisher's Information
    if("try-error" %in% class(i.theta.inv)) i.theta.inv = chol2inv(chol(i.theta + diag(4)))
    ###############################
    # Pre-conditioned MALA Update #
    ###############################
    theta.draw = theta  + e.theta^2 * i.theta.inv %*% g.theta/2 + e.theta * chol(i.theta.inv) %*% matrix(rnorm(4), nc = 1)
    sig2.draw = theta.draw[1]
    tau2.draw = theta.draw[2]
    phis.draw = theta.draw[3]
    phit.draw = theta.draw[4]

    if(sig2.draw < 0 | tau2.draw < 0 | phis.draw < lower.phis | phis.draw > upper.phis | phit.draw < lower.phit | phit.draw > upper.phit){
      accept.prob = -Inf # error-handling
    }else{
      # updated correlation matrix
      if(cov.type == "exponential"){
        R.draw = st_cov_exponential(delta = delta, Delta = Delta,
                                    phis = phis.draw, phit = phit.draw,
                                    sig2 = 1)
      }
      if(cov.type == "matern1"){
        R.draw = st_cov_matern1(delta = delta, Delta = Delta,
                                phis = phis.draw, phit = phit.draw,
                                sig2 = 1)
      }

      if(cov.type == "matern2"){
        R.draw = st_cov_matern2(delta = delta, Delta = Delta,
                                phis = phis.draw, phit = phit.draw,
                                sig2 = 1)
      }

      if(cov.type == "gaussian"){
        R.draw = st_cov_gaussian(delta = delta, Delta = Delta,
                                 phis = phis.draw, phit = phit.draw,
                                 sig2 = 1)
      }
      inv.M.c.draw = chol(sig2.draw * R.draw + tau2.draw * diag(N)) # updated Cholesky for determinant
      inv.M.draw = chol2inv(inv.M.c.draw)

      if(cov.type == "exponential"){
        del.R.del.phis.draw = del.phis.st_cov_exponential(delta = delta, Delta = Delta,
                                                     phis = phis.draw, phit = phit.draw,
                                                     sig2 = sig2)                          # gradient of R w.r.t phis
        del.R.del.phit.draw = del.phit.st_cov_exponential(delta = delta, Delta = Delta,
                                                     phis = phis.draw, phit = phit.draw,
                                                     sig2 = sig2)                          # gradient of R w.r.t phit
      }
      if(cov.type == "matern1"){

        del.R.del.phis.draw = del.phis.st_cov_matern1(delta = delta, Delta = Delta,
                                                 phis = phis.draw, phit = phit.draw,
                                                 sig2 = sig2)                          # gradient of R w.r.t phis
        del.R.del.phit.draw = del.phit.st_cov_matern1(delta = delta, Delta = Delta,
                                                 phis = phis.draw, phit = phit.draw,
                                                 sig2 = sig2)                          # gradient of R w.r.t phit
      }
      if(cov.type == "matern2"){

        del.R.del.phis.draw = del.phis.st_cov_matern2(delta = delta, Delta = Delta,
                                                 phis = phis.draw, phit = phit.draw,
                                                 sig2 = sig2)                          # gradient of R w.r.t phis
        del.R.del.phit.draw = del.phit.st_cov_matern2(delta = delta, Delta = Delta,
                                                 phis = phis.draw, phit = phit.draw,
                                                 sig2 = sig2)                          # gradient of R w.r.t phit
      }
      if(cov.type == "gaussian"){

        del.R.del.phis.draw = del.phis.st_cov_gaussian(delta = delta, Delta = Delta,
                                                  phis = phis.draw, phit = phit.draw,
                                                  sig2 = sig2)                          # gradient of R w.r.t phis
        del.R.del.phit.draw = del.phit.st_cov_gaussian(delta = delta, Delta = Delta,
                                                  phis = phis.draw, phit = phit.draw,
                                                  sig2 = sig2)                          # gradient of R w.r.t phit
      }

      t1 = inv.M.draw %*% R.draw
      t2 = inv.M.draw %*% inv.M.draw
      t3 = inv.M.draw %*% del.R.del.phis.draw
      t4 = inv.M.draw %*% del.R.del.phit.draw
      t5 = inv.M.draw %*% (y - X %*% beta)
      g.sig2.draw = - (a.sigma + 1)/sig2.draw + b.sigma/sig2.draw^2 - sum(diag(t1))/2 + t(y - X %*% beta) %*% t1 %*% t5/2
      g.tau2.draw = - (a.tau + 1)/tau2.draw + b.tau/tau2.draw^2 - sum(diag(inv.M.draw))/2 + t(y - X %*% beta) %*% t2 %*% (y - X %*% beta)/2
      g.phis.draw = - sum(diag(t3))/2 + t(y - X %*% beta) %*% t3 %*% t5/2
      g.phit.draw = - sum(diag(t4))/2 + t(y - X %*% beta) %*% t4 %*% t5/2
      g.theta.draw = matrix(c(g.sig2.draw, g.tau2.draw, g.phis.draw, g.phit.draw), nc = 1)

      i.11.draw =  -(a.sigma + 1)/sig2.draw^2 + 2 * b.sigma/sig2.draw^3 + sum(diag(t1 %*% t1))/2
      i.12.draw = sum(diag(t2 %*% R.draw))/2
      i.13.draw = sum(diag(t3 %*% t1))/2
      i.14.draw = sum(diag(t4 %*% t1))/2
      i.22.draw = -(a.tau + 1)/tau2.draw^2 + 2 * b.tau/tau2.draw^3 + sum(diag(t2))/2
      i.23.draw = sum(diag(inv.M.draw %*% t3))/2
      i.24.draw = sum(diag(inv.M.draw %*% t4))/2
      i.33.draw = sum(diag(t3 %*% t3))/2
      i.34.draw = sum(diag(t4 %*% t3))/2
      i.44.draw = sum(diag(t4 %*% t4))/2
      i.theta.draw = matrix(c(i.11.draw, i.12.draw, i.13.draw, i.14.draw,
                              i.12.draw, i.22.draw, i.23.draw, i.24.draw,
                              i.13.draw, i.23.draw, i.33.draw, i.34.draw,
                              i.14.draw, i.24.draw, i.34.draw, i.44.draw),
                            ncol = 4, nrow = 4, byrow = T)

      i.theta.inv.draw = try(chol2inv(chol(i.theta.draw)), silent = TRUE)
      if("try-error" %in% class(i.theta.inv.draw)) i.theta.inv.draw = chol2inv(chol(i.theta.draw + diag(4)))

      ra = - b.sigma * (1/sig2.draw - 1/sig2) - (a.sigma + 1) * log(sig2.draw/sig2) -
        b.tau * (1/tau2.draw - 1/tau2) - (a.tau + 1) * log(tau2.draw/tau2) -
        (t(y - X %*% beta) %*% inv.M.draw %*% (y - X %*% beta) - t(y - X %*% beta) %*% inv.M %*% (y - X %*% beta))/2 -
        sum(log(diag(inv.M.c.draw)/diag(inv.M.c)))
      rb = t(theta.draw - theta - e.theta^2 * i.theta.inv %*% g.theta/2) %*% chol2inv(chol(e.theta^2 * i.theta.inv)) %*% (theta.draw - theta - e.theta^2 * i.theta.inv %*% g.theta/2)/2 -
        t(theta - theta.draw - e.theta^2 * i.theta.inv.draw %*% g.theta.draw/2) %*% chol2inv(chol(e.theta^2 * i.theta.inv.draw)) %*% (theta - theta.draw - e.theta^2 * i.theta.inv.draw %*% g.theta.draw/2)/2

      accept.prob = min(ra + rb, 0)
    }

    if(log(runif(1)) < accept.prob){
      res_sig2[i] = sig2 = sig2.draw
      res_tau2[i] = tau2 = tau2.draw
      res_phis[i] = phis = phis.draw
      res_phit[i] = phit = phit.draw
      R = R.draw
      inv.M.c = inv.M.c.draw
      inv.M = inv.M.draw
      accept.p = accept.p + 1
    }else{
      res_sig2[i] = sig2
      res_tau2[i] = tau2
      res_phis[i] = phis
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
      accept_m = c(accept_m, accept.p)
      if(verbose){
        if(i > nburn){
          cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(e.theta, digits), "\n",
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
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(e.theta, digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                # "Beta::", round(apply(res_beta[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(e.beta, (report %% 10)),"\n",
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
      step = 1
      accept.p = max(0.25, min(accept.p, 0.75))
      if(accept.p > 0.65) step = accept.p/0.65
      else if(accept.p < 0.45) step = accept.p/0.45

      e.theta = e.theta * step

      accept.p = 0
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
# spt_model0 = hlmBayes_mala_pc.spt(y = y, coords = coords, t = t, cov.type = "exponential")
# spt_model1 = hlmBayes_mala_pc.spt(y = y, coords = coords, t = t, cov.type = "matern1")
# spt_model2 = hlmBayes_mala_pc.spt(y = y, coords = coords, t = t, cov.type = "matern2")
# spt_model3 = hlmBayes_mala_pc.spt(y = y, coords = coords, t = t, cov.type = "gaussian")
# save(mcmc.list, file = "mala-precon-mcmc.RData")
# plot(c(trgt_fn.init,trgt_fn), type="l", col = "orange", ylab = "Target Probability")
# grid()
# par(mfrow=c(5,3))
# plot_mcmc(samples = res_beta[nburn:niter,], cnames = "beta0")
# plot_mcmc(samples = res_sig2[nburn:niter], cnames = "sigma2")
# plot_mcmc(samples = res_tau2[nburn:niter], cnames = "tau2")
# plot_mcmc(samples = res_phis[nburn:niter], cnames = "phis")
# plot_mcmc(samples = res_phit[nburn:niter], cnames = "phit")

