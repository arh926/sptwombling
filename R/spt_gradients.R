#' Posterior sampling of spatiotemporal differential processes
#'
#' Performs 1-for-1 sampling of spatiotemporal differential processes
#'
#' @param model a list containing posterior samples from a collapsed or Gibbs sampler that fits a spatiotemporal Bayesian hierarchical model
#' @param grid.points can be an equi-spaced grid or in-fill or extrapolated points
#' @param cov.type three choices of covariance: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @param nbatch batching for posterior samples. Defaults at 300.
#' @param plots logical for plotting. Defaults to FALSE.
#' @param only.grad.no.curv logical for only gradients and no curvature estimation
#' @param true true values, if available. Should be supplied as a list.
#' @keywords spt_gradients
#' @import parallel magic latex2exp coda
#' @export
spt_gradients <- function(model = NULL, # should be a list of MCMC results post burn-in
                          grid.points = NULL,
                          cov.type = c("matern1", "matern2", "gaussian"),
                          nbatch = NULL,
                          plots = FALSE,
                          only.grad.no.curv = FALSE,
                          true = NULL){ # supply as a list

  if(is.null(grid.points)){
    grid.points = expand.grid(x = seq(0, 1, by = 0.1),
                              y = seq(0, 1, by = 0.1))
  }
  if(is.null(nbatch)) nbatch = 300

  coords <- model$coords
  d <- ncol(coords)
  t <- model$t

  delta <- as.matrix(dist(t))
  Delta <- as.matrix(dist(coords))

  Ns <- nrow(coords)
  Nt <- length(t)
  N <- Ns * Nt


  id.mcmc <- 1:length(model$sig2)

  ncores <- detectCores() - 1
  samp.list <- split(id.mcmc, ceiling(seq_along(id.mcmc)/(length(id.mcmc)/ncores)))
  parallel.index <- 1:ncores

  dist.s0 <- sapply(1:nrow(grid.points),function(y) apply(coords, 1, function(x) sqrt(sum((x - grid.points[y,])^2)) ))
  delta.s0 <- sapply(1:nrow(grid.points),function(y) do.call(rbind, t(apply(coords, 1, function(x) x - grid.points[y,]))), simplify = FALSE)

  ################################
  # Configurations and constants #
  ################################
  I.d = diag(d)
  I.tilde = matrix(c(3, 0, 1,
                     0, 1, 0,
                     1, 0, 3),
                   nrow = 3, ncol = 3, byrow = TRUE)
  I.star = c(1, 0, 1)

  if(cov.type == "matern1"){
    ngrad = 2 * (1 + d) - 1
    results.grad.st <- sapply(t, function(t0){
      delta.t0 <- t - t0

      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]

        sig2 <- model$sig2[samp.x]
        phis <- model$phis[samp.x]
        phit <- model$phit[samp.x]
        Z <- model$z[samp.x,]

        n.mcmc <- length(sig2)

        mcmc.grad <- list()

        for(i.mcmc in 1:n.mcmc){

          R.Z <- st_cov_matern1(delta = delta,
                                Delta = Delta,
                                lphis = log(phis[i.mcmc]),
                                lphit = log(phit[i.mcmc]),
                                lsig2 = 0)

          #############################################
          # Error Handling:: in case R.Z is singular #
          #############################################
          R.Z.in <- try(chol2inv(chol(R.Z)), silent = TRUE)
          if("try-error" %in% class(R.Z.in)) R.Z.in <- chol2inv(chol(R.Z + 1e-4 * diag(N)))

          grad.est <- matrix(NA, nrow = nrow(grid.points), ncol = ngrad)

          A.t <- (phit[i.mcmc]^2 * delta.t0^2 + 1)
          for(i in 1:nrow(grid.points)){

            ###########################
            # gradients and curvature #
            ###########################
            # \nabla_s
            nabKs <- as.vector(unlist(sapply(dist.s0[,i], function(x) -3 * phis[i.mcmc]^2 * exp(-sqrt(3) * phis[i.mcmc] * x/sqrt(A.t))/A.t^2, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKs) <- NULL
            # \nabla_t
            nabKt <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phit[i.mcmc]^2 * delta.t0 * exp(-sqrt(3) * phis[i.mcmc] * x/sqrt(A.t)) * (1/A.t^2 + sqrt(3) * phis[i.mcmc] * x/A.t^(5/2) - 3 * phis[i.mcmc]^2 * x^2/A.t^3/2), simplify = FALSE)))
            # \nabla_t\nabla_s
            nabKst <- as.vector(unlist(sapply(dist.s0[,i], function(x) 6 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(3) * phis[i.mcmc] * x/sqrt(A.t)) * (2/A.t^3 - sqrt(3) * phis[i.mcmc] * x/A.t^(7/2)/2) * delta.t0, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKst) <- NULL

            nabK0 <- as.matrix(cbind(nabKs, nabKt, nabKst))
            nabK0.t <- as.matrix(cbind(-nabKs, -nabKt, nabKst))

            HK <- adiag(3 * phis[i.mcmc]^2 * I.d,
                        2 * phit[i.mcmc]^2,
                        12 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.d)

            m.tmp <- t(crossprod(nabK0.t, R.Z.in))
            muZ0 <- crossprod(m.tmp, Z[i.mcmc,])
            varZ0 <- HK - crossprod(m.tmp, nabK0)

            chol.Z0 = try(chol(varZ0))
            if("try-error" %in% class(chol.Z0)){
              grad.est[i,] = muZ0
            }else{
              grad.est[i,] <- sqrt(sig2[i.mcmc]) * crossprod(chol(varZ0), rnorm(ngrad)) + muZ0
            }
          }

          mcmc.grad[[i.mcmc]] <- grad.est

          if(i.mcmc/100 == floor(i.mcmc/100)){
            cat('Iteration: ',i.mcmc,'\n')
          }
        }
        return(mcmc.grad)
      }, mc.cores = ncores)
      return(results.grad)
    }, simplify = FALSE)
  }else if(cov.type == "matern2"){
    ngrad = 3 * (1 + d + d * (d + 1)/2) - 1

    results.grad.st <- sapply(t, function(t0){
      delta.t0 <- t - t0

      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]

        sig2 <- model$sig2[samp.x]
        phis <- model$phis[samp.x]
        phit <- model$phit[samp.x]
        Z <- model$z[samp.x,] # z + beta

        n.thin.mcmc <- length(sig2)

        mcmc.grad <- list()

        for(i.mcmc in 1:n.thin.mcmc){

          R.Z <- st_cov_matern2(delta = delta,
                                Delta = Delta,
                                lphis = log(phis[i.mcmc]),
                                lphit = log(phit[i.mcmc]),
                                lsig2 = 0)
          #############################################
          # Error Handling:: in case R.Z is singular #
          #############################################
          R.Z.in <- try(chol2inv(chol(R.Z)), silent = TRUE)
          if("try-error" %in% class(R.Z.in)) R.Z.in <- chol2inv(chol(R.Z + 1e-4 * diag(N)))

          grad.est <- matrix(NA, nrow = nrow(grid.points), ncol = ngrad)

          A.t <- (phit[i.mcmc]^2 * delta.t0^2 + 1)
          for(i in 1:nrow(grid.points)){
            # Data Object:: Matrix of: ||\Delta||, \Delta_1^2, \Delta_1\Delta_2, \Delta_2^2
            dist.delta.mat <- cbind(dist.s0[,i], delta.s0[[i]][,1]^2, delta.s0[[i]][,1] * delta.s0[[i]][,2], delta.s0[[i]][,2]^2)
            ###########################
            # gradients and curvature #
            ###########################
            # \nabla_s
            nabKs <- as.vector(unlist(sapply(dist.s0[,i], function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x/sqrt(A.t))/A.t^2, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKs) <- NULL
            # \nabla_s^2
            nabK2s.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[2]/A.t)/A.t^2, simplify = FALSE)))
            nabK2s.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 25/3 * phis[i.mcmc]^4 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * x[3]/A.t^3, simplify = FALSE)))
            nabK2s.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[4]/A.t)/A.t^2, simplify = FALSE)))

            # \nabla_t
            nabKt <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) + 5/6 * phis[i.mcmc]^2 * x^2/A.t - 5 * sqrt(5)/6 * phis[i.mcmc]^3 * x^3/A.t^(3/2)) * delta.t0/A.t^2, simplify = FALSE)))
            # \nabla_t\nabla_s
            nabKst <- as.vector(unlist(sapply(dist.s0[,i], function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 5 * phis[i.mcmc]^2 * x^2/A.t)/A.t^3 * delta.t0, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKst) <- NULL
            # \nabla_t\nabla_s^2
            nabK2st.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[2]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[2] * x[1]/A.t^(3/2)) * delta.t0/A.t^3, simplify = FALSE)))
            nabK2st.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -25/3 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (6 - sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * x[3] * delta.t0/A.t^4, simplify = FALSE)))
            nabK2st.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[4]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[4] * x[1]/A.t^(3/2)) * delta.t0/A.t^3, simplify = FALSE)))

            # \nabla_t^2
            nabK2t <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * ((1 + sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) + 5/6 * phis[i.mcmc]^2 * x^2/A.t - 5 * sqrt(5)/6 * phis[i.mcmc]^3 * x^3/A.t^(3/2)) -
                                                                                                                                             (4 + 4 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 20 * sqrt(5)/3 * phis[i.mcmc]^3 * x^3/A.t^(3/2) + 25/6 * phis[i.mcmc]^4 * x^4/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^2, simplify = FALSE)))
            # \nabla_t^2\nabla_s
            nabKs2t <- as.vector(unlist(sapply(dist.s0[,i], function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * ((4 + 4 * sqrt(5) * phis[i.mcmc] * x/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x^2/A.t) -
                                                                                                                                                                (24 + 24 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 60 * phis[i.mcmc]^2 * x^2/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x^3/A.t^(3/2)) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKs2t) <- NULL
            # \nabla_t^2\nabla_s^2
            nabK2s2t.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * ((4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[2]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[2] * x[1]/A.t^(3/2)) -
                                                                                                                                                                            (24 + 24 * sqrt(5) * phis[i.mcmc] * x[1]/A.t^(1/2) - 60 * phis[i.mcmc]^2 * x[1]^2/A.t - 240 * phis[i.mcmc]^2 * x[2]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[1]^3/A.t^(3/2) + 75 * sqrt(5) * phis[i.mcmc]^3 * x[2] * x[1]/A.t^(3/2) - 25 * phis[i.mcmc]^4 * x[2] * x[1]^2/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3,
                                                  simplify = FALSE)))
            nabK2s2t.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -25/3 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * ((6 - sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) -
                                                                                                                                                                              (48 - 15 * sqrt(5) * phis[i.mcmc] * x[1]/A.t^(1/2) + 5 * phis[i.mcmc]^2 * x[1]^2/A.t) * phit[i.mcmc]^2 * delta.t0^2/A.t) * x[3]/A.t^4,
                                                  simplify = FALSE)))
            nabK2s2t.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * ((4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[4]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[4] * x[1]/A.t^(3/2)) -
                                                                                                                                                                            (24 + 24 * sqrt(5) * phis[i.mcmc] * x[1]/A.t^(1/2) - 60 * phis[i.mcmc]^2 * x[1]^2/A.t - 240 * phis[i.mcmc]^2 * x[4]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[1]^3/A.t^(3/2) + 75 * sqrt(5) * phis[i.mcmc]^3 * x[4] * x[1]/A.t^(3/2) - 25 * phis[i.mcmc]^4 * x[4] * x[1]^2/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3,
                                                  simplify = FALSE)))

            ####################
            # HK:: Calculation #
            ####################
            # Diagonal Entries
            nabK0 <- as.matrix(cbind(nabKs,  nabK2s.11,  nabK2s.12,  nabK2s.22,
                                     nabKt, nabKst, nabK2st.11, nabK2st.12, nabK2st.22,
                                     nabK2t,  nabKs2t, nabK2s2t.11, nabK2s2t.12, nabK2s2t.22))
            nabK0.t <- as.matrix(cbind(-nabKs,  nabK2s.11,  nabK2s.12,  nabK2s.22,
                                       -nabKt, nabKst, -nabK2st.11, -nabK2st.12, -nabK2st.22,
                                       nabK2t,  -nabKs2t, nabK2s2t.11, nabK2s2t.12, nabK2s2t.22))
            colnames(nabK0) <- colnames(nabK0.t) <- c("nabKs.1", "nabKs.2",  "nabK2s.11",  "nabK2s.12",  "nabK2s.22",
                                                      "nabKt", "nabKst.1", "nabKst.2", "nabK2st.11", "nabK2st.12", "nabK2st.22",
                                                      "nabK2t",  "nabKs2t.1", "nabKs2t.2", "nabK2s2t.11", "nabK2s2t.12", "nabK2s2t.22")

            # HK is not diagonal
            # Comment:: iijj terms for \nabK2s entries are non-zero
            HK <- adiag(5/3 * phis[i.mcmc]^2 * I.d,
                        25/3 * phis[i.mcmc]^4 * I.tilde, # changed here 05/12/24 added/3
                        2 * phit[i.mcmc]^2,
                        20/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.d,
                        50 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * I.tilde,
                        24 * phit[i.mcmc]^4,
                        120 * phis[i.mcmc]^2 * phit[i.mcmc]^4 * I.d,
                        1200 * phis[i.mcmc]^4 * phit[i.mcmc]^4 * I.tilde)


            # Non-zero off-diagonal entries for HK
            HK[1:2, 13:14] <- HK[13:14, 1:2] <- -20/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.d
            HK[3:5, 12] <- HK[12, 3:5] <- 20/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.star
            HK[3:5, 15:17] <-  HK[15:17, 3:5] <- -50 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * I.tilde
            HK[6, 9:11] <- HK[9:11, 6] <- -20/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.star
            HK[12, 15:17] <- HK[15:17, 12] <- -120 * phis[i.mcmc]^2 * phit[i.mcmc]^4 * I.star

            ##################
            # posterior mean #
            ##################
            m.tmp <- t(crossprod(nabK0.t, R.Z.in))
            muZ0 <- crossprod(m.tmp, Z[i.mcmc,])
            ######################
            # posterior variance #
            ######################
            varZ0 <- HK - crossprod(m.tmp, nabK0)
            ############
            # sampling #
            ############
            chol.varZ0 = try(chol(varZ0), silent = TRUE)
            if("try-error" %in% class(chol.varZ0)){
              # Alternative Approach: takes longer
              # MASS::mvrnorm(1, mu = muZ0, Sigma = sig2[i.mcmc] * varZ0, tol = 1e-1)
              grad.est[i,] <- muZ0
            }else{
              grad.est[i,] <- sqrt(sig2[i.mcmc]) * crossprod(chol.varZ0, rnorm(ngrad)) + muZ0
            }
          }

          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }, mc.cores = ncores)
      return(results.grad)
    }, simplify = FALSE)
  }else if(cov.type == "gaussian"){
    ngrad = 3 * (1 + d + d * (d + 1)/2) - 1

    results.grad.st <- sapply(t, function(t0){
      delta.t0 <- t - t0

      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]

        sig2 <- model$sig2[samp.x]
        phis <- model$phis[samp.x]
        phit <- model$phit[samp.x]
        Z <- model$z[samp.x,] # z + beta

        n.thin.mcmc <- length(sig2)

        mcmc.grad <- list()

        for(i.mcmc in 1:n.thin.mcmc){

          R.Z <- st_cov_gaussian(delta = delta,
                                 Delta = Delta,
                                 lphis = log(phis[i.mcmc]),
                                 lphit = log(phit[i.mcmc]),
                                 lsig2 = 0)

          #############################################
          # Error Handling:: in case R.Z is singular #
          #############################################
          R.Z.in <- try(chol2inv(chol(R.Z)), silent = TRUE)
          if("try-error" %in% class(R.Z.in)) R.Z.in <- chol2inv(chol(R.Z + 1e-6 * diag(N)))


          grad.est <- matrix(NA, nrow = nrow(grid.points), ncol = ngrad)

          A.t <- (phit[i.mcmc]^2 * delta.t0^2 + 1)
          for(i in 1:nrow(grid.points)){
            # Data Object:: Matrix of: ||\Delta||, \Delta_1^2, \Delta_1\Delta_2, \Delta_2^2
            dist.delta.mat <- cbind(dist.s0[,i], delta.s0[[i]][,1]^2, delta.s0[[i]][,1] * delta.s0[[i]][,2], delta.s0[[i]][,2]^2)
            ###########################
            # gradients and curvature #
            ###########################
            # \nabla_s
            nabKs <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phis[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x^2/A.t)/A.t^2, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKs) <- NULL
            # \nabla_s^2
            nabK2s.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -2 * phis[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * (1 - 2 * phis[i.mcmc]^2 * x[2]/A.t)/A.t^2, simplify = FALSE)))
            nabK2s.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 4 * phis[i.mcmc]^4 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * x[3]/A.t^3, simplify = FALSE)))
            nabK2s.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -2 * phis[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * (1 - 2 * phis[i.mcmc]^2 * x[4]/A.t)/A.t^2, simplify = FALSE)))
            # \nabla_t
            nabKt <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x^2/A.t) * (1 - phis[i.mcmc]^2 * x^2/A.t) * delta.t0/A.t^2, simplify = FALSE)))
            # \nabla_t\nabla_s
            nabKst <- as.vector(unlist(sapply(dist.s0[,i], function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x^2/A.t) * (2 - phis[i.mcmc]^2 * x^2/A.t)/A.t^3 * delta.t0, simplify = F))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKst) <- NULL
            # \nabla_t\nabla_s^2
            nabK2st.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * (2 - phis[i.mcmc]^2 * x[1]^2/A.t - 6 * phis[i.mcmc]^2 * x[2]/A.t + 2 * phis[i.mcmc]^4 * x[2] * x[1]^2/A.t^2) * delta.t0/A.t^3, simplify = FALSE)))
            nabK2st.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -8 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * (3 - phis[i.mcmc]^2 * x[1]^2/A.t) * x[3] * delta.t0/A.t^4, simplify = FALSE)))
            nabK2st.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * (2 - phis[i.mcmc]^2 * x[1]^2/A.t - 6 * phis[i.mcmc]^2 * x[4]/A.t + 2 * phis[i.mcmc]^4 * x[4] * x[1]^2/A.t^2) * delta.t0/A.t^3, simplify = FALSE)))
            # \nabla_t^2
            nabK2t <- as.vector(unlist(sapply(dist.s0[,i], function(x) -2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x^2/A.t) * ((1 - phis[i.mcmc]^2 * x^2/A.t) -
                                                                                                                                 2 * (2 - 4 * phis[i.mcmc]^2 * x^2/A.t + phis[i.mcmc]^4 * x^4/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^2, simplify = FALSE)))
            # \nabla_t^2\nabla_s
            nabKs2t <- as.vector(unlist(sapply(dist.s0[,i], function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x^2/A.t) * ((2 - phis[i.mcmc]^2 * x^2/A.t) -
                                                                                                                                                2 * (6 - 6 * phis[i.mcmc]^2 * x^2/A.t + phis[i.mcmc]^4 * x^4/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3, simplify = FALSE))) * delta.s0[[i]][rep(seq_len(Ns), each = Nt),]; rownames(nabKs2t) <- NULL
            # \nabla_t^2\nabla_s^2
            nabK2s2t.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * ((2 - phis[i.mcmc]^2 * x[1]^2/A.t - 6 * phis[i.mcmc]^2 * x[2]/A.t + 2 * phis[i.mcmc]^4 * x[2] * x[1]^2/A.t^2) -
                                                                                                                                                              2 * (6 - 6 * phis[i.mcmc]^2 * x[1]^2/A.t - 24 * phis[i.mcmc]^2 * x[2]/A.t + phis[i.mcmc]^4 * x[1]^4/A.t^2 + 16 * phis[i.mcmc]^4 * x[2] * x[1]^2/A.t^2 - 2 * phis[i.mcmc]^6 * x[2] * x[1]^4/A.t^3) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3,
                                                  simplify = FALSE)))
            nabK2s2t.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -8 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * ((3 - phis[i.mcmc]^2 * x[1]^2/A.t) -
                                                                                                                                                               2 * (12 - 8 * phis[i.mcmc]^2 * x[1]^2/A.t + phis[i.mcmc]^4 * x[1]^4/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t) * x[3]/A.t^4,
                                                  simplify = FALSE)))
            nabK2s2t.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 4 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-phis[i.mcmc]^2 * x[1]^2/A.t) * ((2 - phis[i.mcmc]^2 * x[1]^2/A.t - 6 * phis[i.mcmc]^2 * x[4]/A.t + 2 * phis[i.mcmc]^4 * x[4] * x[1]^2/A.t^2) -
                                                                                                                                                              2 * (6 - 6 * phis[i.mcmc]^2 * x[1]^2/A.t - 24 * phis[i.mcmc]^2 * x[4]/A.t + phis[i.mcmc]^4 * x[1]^4/A.t^2 + 16 * phis[i.mcmc]^4 * x[4] * x[1]^2/A.t^2 - 2 * phis[i.mcmc]^6 * x[4] * x[1]^4/A.t^3) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3,
                                                  simplify = FALSE)))


            nabK0 <- as.matrix(cbind(nabKs,  nabK2s.11,  nabK2s.12,  nabK2s.22,
                                     nabKt, nabKst, nabK2st.11, nabK2st.12, nabK2st.22,
                                     nabK2t,  nabKs2t, nabK2s2t.11, nabK2s2t.12, nabK2s2t.22))
            nabK0.t <- as.matrix(cbind(-nabKs,  nabK2s.11,  nabK2s.12,  nabK2s.22,
                                       -nabKt, nabKst, -nabK2st.11, -nabK2st.12, -nabK2st.22,
                                       nabK2t,  -nabKs2t, nabK2s2t.11, nabK2s2t.12, nabK2s2t.22))
            colnames(nabK0) <- colnames(nabK0.t) <- c("nabKs.1", "nabKs.2",  "nabK2s.11",  "nabK2s.12",  "nabK2s.22",
                                                      "nabKt", "nabKst.1", "nabKst.2", "nabK2st.11", "nabK2st.12", "nabK2st.22",
                                                      "nabK2t",  "nabKs2t.1", "nabKs2t.2", "nabK2s2t.11", "nabK2s2t.12", "nabK2s2t.22")

            HK <- adiag(2 * phis[i.mcmc]^2 * I.d,
                        4 * phis[i.mcmc]^4 * I.tilde,
                        2 * phit[i.mcmc]^2,
                        8 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.d,
                        24 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * I.tilde,
                        24 * phit[i.mcmc]^4,
                        144 * phis[i.mcmc]^2 * phit[i.mcmc]^4 * I.d,
                        576 * phis[i.mcmc]^4 * phit[i.mcmc]^4 * I.tilde)

            # non-zero entries for HK
            HK[1:2, 13:14] <- HK[13:14, 1:2] <- -8 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.d
            HK[3:5, 12] <- HK[12, 3:5] <- 8 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.star
            HK[3:5, 15:17] <-  HK[15:17, 3:5] <- -24 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * I.tilde
            HK[6, 9:11] <- HK[9:11, 6] <- -8 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * I.star
            HK[12, 15:17] <- HK[15:17, 12] <- -144 * phis[i.mcmc]^2 * phit[i.mcmc]^4 * I.star

            # posterior mean
            m.tmp <- t(crossprod(nabK0.t, R.Z.in))
            muZ0 <- crossprod(m.tmp, Z[i.mcmc,])
            # posterior variance
            varZ0 <- HK - crossprod(m.tmp, nabK0)
            # sampling
            chol.varZ0 = try(chol(varZ0), silent = TRUE)
            if("try-error" %in% class(chol.varZ0)){
              # Alternative Approach
              # MASS::mvrnorm(1, mu = muZ0, Sigma = sig2[i.mcmc] * varZ0, tol = 1e-1)
              grad.est[i,] <- muZ0
            }else{
              grad.est[i,] <- sqrt(sig2[i.mcmc]) * crossprod(chol.varZ0, rnorm(ngrad)) + muZ0
            }
          }

          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }, mc.cores = ncores)
      return(results.grad)
    }, simplify = FALSE)
  }

  if(is.null(true)){
    plot.fn.cp <- lapply(t, function(x){
      results.grad = results.grad.st[[x]]
      grad.sx.mcmc.f <- grad.sy.mcmc.f <- curv.sxx.mcmc.f <- curv.sxy.mcmc.f <- curv.syy.mcmc.f <- grad.t.mcmc.f <- grad.sxt.mcmc.f <- grad.syt.mcmc.f <- curv.sxxt.mcmc.f <- curv.sxyt.mcmc.f <- curv.syyt.mcmc.f <- grad.tt.mcmc.f <- grad.sxtt.mcmc.f <- grad.sytt.mcmc.f <- curv.sxxtt.mcmc.f <- curv.sxytt.mcmc.f <- curv.syytt.mcmc.f <- c()
      for(i in 1:length(results.grad)){
        grad.sx.mcmc.f <- rbind(grad.sx.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,1])))
        grad.sy.mcmc.f <- rbind(grad.sy.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,2])))
        curv.sxx.mcmc.f <- rbind(curv.sxx.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,3])))
        curv.sxy.mcmc.f <- rbind(curv.sxy.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,4])))
        curv.syy.mcmc.f <- rbind(curv.syy.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,5])))
        grad.t.mcmc.f <- rbind(grad.t.mcmc.f,
                               do.call(rbind,lapply(results.grad[[i]],function(x) x[,6])))
        grad.sxt.mcmc.f <- rbind(grad.sxt.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,7])))
        grad.syt.mcmc.f <- rbind(grad.syt.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,8])))
        curv.sxxt.mcmc.f <- rbind(curv.sxxt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,9])))
        curv.sxyt.mcmc.f <- rbind(curv.sxyt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,10])))
        curv.syyt.mcmc.f <- rbind(curv.syyt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,11])))
        grad.tt.mcmc.f <- rbind(grad.tt.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,12])))
        grad.sxtt.mcmc.f <- rbind(grad.sxtt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,13])))
        grad.sytt.mcmc.f <- rbind(grad.sytt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,14])))
        curv.sxxtt.mcmc.f <- rbind(curv.sxxtt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,15])))
        curv.sxytt.mcmc.f <- rbind(curv.sxytt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,16])))
        curv.syytt.mcmc.f <- rbind(curv.syytt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,17])))
      }

      slist <- split(id.mcmc, ceiling(seq_along(id.mcmc)/(length(id.mcmc)/nbatch))) # change this to increase CP

      grad.sx.temp <- t(sapply(slist, function(x) apply(grad.sx.mcmc.f[x,], 2, median)))
      grad.sy.temp <- t(sapply(slist, function(x) apply(grad.sy.mcmc.f[x,], 2, median)))
      curv.sxx.temp <- t(sapply(slist, function(x) apply(curv.sxx.mcmc.f[x,], 2, median)))
      curv.sxy.temp <- t(sapply(slist, function(x) apply(curv.sxy.mcmc.f[x,], 2, median)))
      curv.syy.temp <- t(sapply(slist, function(x) apply(curv.syy.mcmc.f[x,], 2, median)))
      grad.t.temp <- t(sapply(slist, function(x) apply(grad.t.mcmc.f[x,], 2, median)))
      grad.sxt.temp <- t(sapply(slist, function(x) apply(grad.sxt.mcmc.f[x,], 2, median)))
      grad.syt.temp <- t(sapply(slist, function(x) apply(grad.syt.mcmc.f[x,], 2, median)))
      curv.sxxt.temp <- t(sapply(slist, function(x) apply(curv.sxxt.mcmc.f[x,], 2, median)))
      curv.sxyt.temp <- t(sapply(slist, function(x) apply(curv.sxyt.mcmc.f[x,], 2, median)))
      curv.syyt.temp <- t(sapply(slist, function(x) apply(curv.syyt.mcmc.f[x,], 2, median)))
      grad.tt.temp <- t(sapply(slist, function(x) apply(grad.tt.mcmc.f[x,], 2, median)))
      grad.sxtt.temp <- t(sapply(slist, function(x) apply(grad.sxtt.mcmc.f[x,], 2, median)))
      grad.sytt.temp <- t(sapply(slist, function(x) apply(grad.sytt.mcmc.f[x,], 2, median)))
      curv.sxxtt.temp <- t(sapply(slist, function(x) apply(curv.sxxtt.mcmc.f[x,], 2, median)))
      curv.sxytt.temp <- t(sapply(slist, function(x) apply(curv.sxytt.mcmc.f[x,], 2, median)))
      curv.syytt.temp <- t(sapply(slist, function(x) apply(curv.syytt.mcmc.f[x,], 2, median)))

      grad.sx.mcmc <- coda::as.mcmc(grad.sx.temp)
      grad.sy.mcmc <- coda::as.mcmc(grad.sy.temp)
      curv.sxx.mcmc <- coda::as.mcmc(curv.sxx.temp)
      curv.sxy.mcmc <- coda::as.mcmc(curv.sxy.temp)
      curv.syy.mcmc <- coda::as.mcmc(curv.syy.temp)
      grad.t.mcmc <- coda::as.mcmc(grad.t.temp)
      grad.sxt.mcmc <- coda::as.mcmc(grad.sxt.temp)
      grad.syt.mcmc <- coda::as.mcmc(grad.syt.temp)
      curv.sxxt.mcmc <- coda::as.mcmc(curv.sxxt.temp)
      curv.sxyt.mcmc <- coda::as.mcmc(curv.sxyt.temp)
      curv.syyt.mcmc <- coda::as.mcmc(curv.syyt.temp)
      grad.tt.mcmc <- coda::as.mcmc(grad.tt.temp)
      grad.sxtt.mcmc <- coda::as.mcmc(grad.sxtt.temp)
      grad.sytt.mcmc <- coda::as.mcmc(grad.sytt.temp)
      curv.sxxtt.mcmc <- coda::as.mcmc(curv.sxxtt.temp)
      curv.sxytt.mcmc <- coda::as.mcmc(curv.sxytt.temp)
      curv.syytt.mcmc <- coda::as.mcmc(curv.syytt.temp)

      grad.mcmc = list(grad.sx.mcmc = grad.sx.mcmc,
                       grad.sy.mcmc = grad.sy.mcmc,
                       curv.sxx.mcmc = curv.sxx.mcmc,
                       curv.sxy.mcmc = curv.sxy.mcmc,
                       curv.syy.mcmc = curv.syy.mcmc,
                       grad.t.mcmc = grad.t.mcmc,
                       grad.sxt.mcmc = grad.sxt.mcmc,
                       grad.syt.mcmc = grad.syt.mcmc,
                       curv.sxxt.mcmc = curv.sxxt.mcmc,
                       curv.sxyt.mcmc = curv.sxyt.mcmc,
                       curv.syyt.mcmc = curv.syyt.mcmc,
                       grad.tt.mcmc = grad.tt.mcmc,
                       grad.sxtt.mcmc = grad.sxtt.mcmc,
                       grad.sytt.mcmc = grad.sytt.mcmc,
                       curv.sxxtt.mcmc = curv.sxxtt.mcmc,
                       curv.sxytt.mcmc = curv.sxytt.mcmc,
                       curv.syytt.mcmc = curv.syytt.mcmc)

      grad.sx.est <- apply(grad.sx.mcmc, 2, median)
      grad.sy.est <- apply(grad.sy.mcmc, 2, median)
      curv.sxx.est <- apply(curv.sxx.mcmc, 2, median)
      curv.sxy.est <- apply(curv.sxy.mcmc, 2, median)
      curv.syy.est <- apply(curv.syy.mcmc, 2, median)
      grad.t.est <- apply(grad.t.mcmc, 2, median)
      grad.sxt.est <- apply(grad.sxt.mcmc, 2, median)
      grad.syt.est <- apply(grad.syt.mcmc, 2, median)
      curv.sxxt.est <- apply(curv.sxxt.mcmc, 2, median)
      curv.sxyt.est <- apply(curv.sxyt.mcmc, 2, median)
      curv.syyt.est <- apply(curv.syyt.mcmc, 2, median)
      grad.tt.est <- apply(grad.tt.mcmc, 2, median)
      grad.sxtt.est <- apply(grad.sxtt.mcmc, 2, median)
      grad.sytt.est <- apply(grad.sytt.mcmc, 2, median)
      curv.sxxtt.est <- apply(curv.sxxtt.mcmc, 2, median)
      curv.sxytt.est <- apply(curv.sxytt.mcmc, 2, median)
      curv.syytt.est <- apply(curv.syytt.mcmc, 2, median)

      col.grad <- c()
      grad.sx.hpd <- data.frame(round(cbind(grad.sx.est, coda::HPDinterval(grad.sx.mcmc)), 2))
      colnames(grad.sx.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.sx.hpd$signif <- apply(grad.sx.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.sx.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.sy.hpd <- data.frame(round(cbind(grad.sy.est, coda::HPDinterval(grad.sy.mcmc)), 2))
      colnames(grad.sy.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.sy.hpd$signif <- apply(grad.sy.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.sy.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxx.hpd <- data.frame(round(cbind(curv.sxx.est, coda::HPDinterval(curv.sxx.mcmc)), 2))
      colnames(curv.sxx.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxx.hpd$signif <- apply(curv.sxx.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxx.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxy.hpd <- data.frame(round(cbind(curv.sxy.est, coda::HPDinterval(curv.sxy.mcmc)), 2))
      colnames(curv.sxy.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxy.hpd$signif <- apply(curv.sxy.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxy.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.syy.hpd <- data.frame(round(cbind(curv.syy.est, coda::HPDinterval(curv.syy.mcmc)), 2))
      colnames(curv.syy.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.syy.hpd$signif <- apply(curv.syy.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.syy.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.t.hpd <- data.frame(round(cbind(grad.t.est, coda::HPDinterval(grad.t.mcmc)), 2))
      colnames(grad.t.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.t.hpd$signif <- apply(grad.t.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.t.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.sxt.hpd <- data.frame(round(cbind(grad.sxt.est, coda::HPDinterval(grad.sxt.mcmc)), 2))
      colnames(grad.sxt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.sxt.hpd$signif <- apply(grad.sxt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.sxt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.syt.hpd <- data.frame(round(cbind(grad.syt.est, coda::HPDinterval(grad.syt.mcmc)), 2))
      colnames(grad.syt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.syt.hpd$signif <- apply(grad.syt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.syt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxxt.hpd <- data.frame(round(cbind(curv.sxxt.est, coda::HPDinterval(curv.sxxt.mcmc)), 2))
      colnames(curv.sxxt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxxt.hpd$signif <- apply(curv.sxxt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxxt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxyt.hpd <- data.frame(round(cbind(curv.sxyt.est, coda::HPDinterval(curv.sxyt.mcmc)), 2))
      colnames(curv.sxyt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxyt.hpd$signif <- apply(curv.sxyt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxyt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.syyt.hpd <- data.frame(round(cbind(curv.syyt.est, coda::HPDinterval(curv.syyt.mcmc)), 2))
      colnames(curv.syyt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.syyt.hpd$signif <- apply(curv.syyt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.syyt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.tt.hpd <- data.frame(round(cbind(grad.tt.est, coda::HPDinterval(grad.tt.mcmc)), 2))
      colnames(grad.tt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.tt.hpd$signif <- apply(grad.tt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.tt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.sxtt.hpd <- data.frame(round(cbind(grad.sxtt.est, coda::HPDinterval(grad.sxtt.mcmc)), 2))
      colnames(grad.sxtt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.sxtt.hpd$signif <- apply(grad.sxtt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.sxtt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      grad.sytt.hpd <- data.frame(round(cbind(grad.sytt.est, coda::HPDinterval(grad.sytt.mcmc)), 2))
      colnames(grad.sytt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      grad.sytt.hpd$signif <- apply(grad.sytt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(grad.sytt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxxtt.hpd <- data.frame(round(cbind(curv.sxxtt.est, coda::HPDinterval(curv.sxxtt.mcmc)), 2))
      colnames(curv.sxxtt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxxtt.hpd$signif <- apply(curv.sxxtt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxxtt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.sxytt.hpd <- data.frame(round(cbind(curv.sxytt.est, coda::HPDinterval(curv.sxytt.mcmc)), 2))
      colnames(curv.sxytt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.sxytt.hpd$signif <- apply(curv.sxytt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.sxytt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      curv.syytt.hpd <- data.frame(round(cbind(curv.syytt.est, coda::HPDinterval(curv.syytt.mcmc)), 2))
      colnames(curv.syytt.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      curv.syytt.hpd$signif <- apply(curv.syytt.hpd,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.grad = cbind(col.grad, sapply(curv.syytt.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      }))

      t.seq <- seq(x, N, by = length(t))
      y <- model$y[t.seq]
      ci.y.z <- t(apply(model$z[id.mcmc,t.seq], 2, quantile, probs = c(0.5, 0.025, 0.975)))
      y.hpd <- data.frame(round(cbind(y, ci.y.z), 2))
      colnames(y.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
      y.hpd$signif <- apply(y.hpd, 1, function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      col.y = sapply(y.hpd$signif, function(sig){
        if(sig == 0) return("#FFFFFF")
        if(sig == 1) return("#00FF00")
        if(sig == -1) return("#00FFFF")
      })

      if(plots){
        plots.t = spt_gradients_est_plot(coords = coords,
                                         data_frame = y.hpd[, "est"],
                                         grid.points = grid.points,
                                         gradients = data.frame(nabla.sx = grad.sx.hpd[, "est"],
                                                                nabla.sy = grad.sy.hpd[, "est"],
                                                                nabla.sxx = curv.sxx.hpd[, "est"],
                                                                nabla.sxy = curv.sxy.hpd[, "est"],
                                                                nabla.syy = curv.syy.hpd[, "est"],
                                                                nabla.t = grad.t.hpd[, "est"],
                                                                nabla.sxt = grad.sxt.hpd[, "est"],
                                                                nabla.syt = grad.syt.hpd[, "est"],
                                                                nabla.sxxt = curv.sxxt.hpd[, "est"],
                                                                nabla.sxyt = curv.sxyt.hpd[, "est"],
                                                                nabla.syyt = curv.syyt.hpd[, "est"],
                                                                nabla.tt = grad.tt.hpd[, "est"],
                                                                nabla.sxtt = grad.sxtt.hpd[, "est"],
                                                                nabla.sytt = grad.sytt.hpd[, "est"],
                                                                nabla.sxxtt = curv.sxxtt.hpd[, "est"],
                                                                nabla.sxytt = curv.sxytt.hpd[, "est"],
                                                                nabla.syytt = curv.syytt.hpd[, "est"]),
                                         col.y = col.y, col.grad = col.grad,
                                         only.grad.no.curv = only.grad.no.curv)

        plots = plot_grid(plotlist = plots.t, labels = LETTERS[1:18],
                          label_size = 6, nrow = 6, ncol = 3, hjust = 0.1)

        return(list(plots = plots,
                    grad.mcmc = grad.mcmc,
                    y = y.hpd,
                    sx = grad.sx.hpd,
                    sy = grad.sy.hpd,
                    sxx = curv.sxx.hpd,
                    sxy = curv.sxy.hpd,
                    syy = curv.syy.hpd,
                    t = grad.t.hpd,
                    sxt = grad.sxt.hpd,
                    syt = grad.syt.hpd,
                    sxxt = curv.sxxt.hpd,
                    sxyt = curv.sxyt.hpd,
                    syyt = curv.syyt.hpd,
                    tt = grad.tt.hpd,
                    sxtt = grad.sxtt.hpd,
                    sytt = grad.sytt.hpd,
                    sxxtt = curv.sxxtt.hpd,
                    sxytt = curv.sxytt.hpd,
                    syytt = curv.syytt.hpd))
      }else{
        return(list(grad.mcmc = grad.mcmc,
                    y = y.hpd,
                    sx = grad.sx.hpd,
                    sy = grad.sy.hpd,
                    sxx = curv.sxx.hpd,
                    sxy = curv.sxy.hpd,
                    syy = curv.syy.hpd,
                    t = grad.t.hpd,
                    sxt = grad.sxt.hpd,
                    syt = grad.syt.hpd,
                    sxxt = curv.sxxt.hpd,
                    sxyt = curv.sxyt.hpd,
                    syyt = curv.syyt.hpd,
                    tt = grad.tt.hpd,
                    sxtt = grad.sxtt.hpd,
                    sytt = grad.sytt.hpd,
                    sxxtt = curv.sxxtt.hpd,
                    sxytt = curv.sxytt.hpd,
                    syytt = curv.syytt.hpd))
      }
    })
    return(plot.fn.cp = plot.fn.cp)
  }else{

    true.grad.sx = true$true.grad.sx
    true.grad.sy = true$true.grad.sy
    true.curv.sxx = true$true.curv.sxx
    true.curv.sxy = true$true.curv.sxy
    true.curv.syy = true$true.curv.syy
    true.grad.t = true$true.grad.t
    true.grad.sxt = true$true.grad.sxt
    true.grad.syt = true$true.grad.syt
    true.curv.sxxt = true$true.curv.sxxt
    true.curv.sxyt = true$true.curv.sxyt
    true.curv.syyt = true$true.curv.syyt
    true.grad.tt = true$true.grad.tt
    true.grad.sxtt = true$true.grad.sxtt
    true.grad.sytt = true$true.grad.sytt
    true.curv.sxxtt = true$true.curv.sxxtt
    true.curv.sxytt = true$true.curv.sxytt
    true.curv.syytt = true$true.curv.syytt

    plot.fn.cp <- lapply(t, function(x){
      results.grad = results.grad.st[[x]]
      grad.sx.mcmc.f <- grad.sy.mcmc.f <- curv.sxx.mcmc.f <- curv.sxy.mcmc.f <- curv.syy.mcmc.f <- grad.t.mcmc.f <- grad.sxt.mcmc.f <- grad.syt.mcmc.f <- curv.sxxt.mcmc.f <- curv.sxyt.mcmc.f <- curv.syyt.mcmc.f <- grad.tt.mcmc.f <- grad.sxtt.mcmc.f <- grad.sytt.mcmc.f <- curv.sxxtt.mcmc.f <- curv.sxytt.mcmc.f <- curv.syytt.mcmc.f <- c()
      for(i in 1:length(results.grad)){
        grad.sx.mcmc.f <- rbind(grad.sx.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,1])))
        grad.sy.mcmc.f <- rbind(grad.sy.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,2])))
        curv.sxx.mcmc.f <- rbind(curv.sxx.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,3])))
        curv.sxy.mcmc.f <- rbind(curv.sxy.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,4])))
        curv.syy.mcmc.f <- rbind(curv.syy.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,5])))
        grad.t.mcmc.f <- rbind(grad.t.mcmc.f,
                               do.call(rbind,lapply(results.grad[[i]],function(x) x[,6])))
        grad.sxt.mcmc.f <- rbind(grad.sxt.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,7])))
        grad.syt.mcmc.f <- rbind(grad.syt.mcmc.f,
                                 do.call(rbind,lapply(results.grad[[i]],function(x) x[,8])))
        curv.sxxt.mcmc.f <- rbind(curv.sxxt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,9])))
        curv.sxyt.mcmc.f <- rbind(curv.sxyt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,10])))
        curv.syyt.mcmc.f <- rbind(curv.syyt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,11])))
        grad.tt.mcmc.f <- rbind(grad.tt.mcmc.f,
                                do.call(rbind,lapply(results.grad[[i]],function(x) x[,12])))
        grad.sxtt.mcmc.f <- rbind(grad.sxtt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,13])))
        grad.sytt.mcmc.f <- rbind(grad.sytt.mcmc.f,
                                  do.call(rbind,lapply(results.grad[[i]],function(x) x[,14])))
        curv.sxxtt.mcmc.f <- rbind(curv.sxxtt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,15])))
        curv.sxytt.mcmc.f <- rbind(curv.sxytt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,16])))
        curv.syytt.mcmc.f <- rbind(curv.syytt.mcmc.f,
                                   do.call(rbind,lapply(results.grad[[i]],function(x) x[,17])))
      }

      slist <- split(id.mcmc, ceiling(seq_along(id.mcmc)/(length(id.mcmc)/nbatch))) # change this to increase CP

      grad.sx.temp <- t(sapply(slist, function(x) apply(grad.sx.mcmc.f[x,], 2, median)))
      grad.sy.temp <- t(sapply(slist, function(x) apply(grad.sy.mcmc.f[x,], 2, median)))
      curv.sxx.temp <- t(sapply(slist, function(x) apply(curv.sxx.mcmc.f[x,], 2, median)))
      curv.sxy.temp <- t(sapply(slist, function(x) apply(curv.sxy.mcmc.f[x,], 2, median)))
      curv.syy.temp <- t(sapply(slist, function(x) apply(curv.syy.mcmc.f[x,], 2, median)))
      grad.t.temp <- t(sapply(slist, function(x) apply(grad.t.mcmc.f[x,], 2, median)))
      grad.sxt.temp <- t(sapply(slist, function(x) apply(grad.sxt.mcmc.f[x,], 2, median)))
      grad.syt.temp <- t(sapply(slist, function(x) apply(grad.syt.mcmc.f[x,], 2, median)))
      curv.sxxt.temp <- t(sapply(slist, function(x) apply(curv.sxxt.mcmc.f[x,], 2, median)))
      curv.sxyt.temp <- t(sapply(slist, function(x) apply(curv.sxyt.mcmc.f[x,], 2, median)))
      curv.syyt.temp <- t(sapply(slist, function(x) apply(curv.syyt.mcmc.f[x,], 2, median)))
      grad.tt.temp <- t(sapply(slist, function(x) apply(grad.tt.mcmc.f[x,], 2, median)))
      grad.sxtt.temp <- t(sapply(slist, function(x) apply(grad.sxtt.mcmc.f[x,], 2, median)))
      grad.sytt.temp <- t(sapply(slist, function(x) apply(grad.sytt.mcmc.f[x,], 2, median)))
      curv.sxxtt.temp <- t(sapply(slist, function(x) apply(curv.sxxtt.mcmc.f[x,], 2, median)))
      curv.sxytt.temp <- t(sapply(slist, function(x) apply(curv.sxytt.mcmc.f[x,], 2, median)))
      curv.syytt.temp <- t(sapply(slist, function(x) apply(curv.syytt.mcmc.f[x,], 2, median)))

      grad.sx.mcmc <- coda::as.mcmc(grad.sx.temp)
      grad.sy.mcmc <- coda::as.mcmc(grad.sy.temp)
      curv.sxx.mcmc <- coda::as.mcmc(curv.sxx.temp)
      curv.sxy.mcmc <- coda::as.mcmc(curv.sxy.temp)
      curv.syy.mcmc <- coda::as.mcmc(curv.syy.temp)
      grad.t.mcmc <- coda::as.mcmc(grad.t.temp)
      grad.sxt.mcmc <- coda::as.mcmc(grad.sxt.temp)
      grad.syt.mcmc <- coda::as.mcmc(grad.syt.temp)
      curv.sxxt.mcmc <- coda::as.mcmc(curv.sxxt.temp)
      curv.sxyt.mcmc <- coda::as.mcmc(curv.sxyt.temp)
      curv.syyt.mcmc <- coda::as.mcmc(curv.syyt.temp)
      grad.tt.mcmc <- coda::as.mcmc(grad.tt.temp)
      grad.sxtt.mcmc <- coda::as.mcmc(grad.sxtt.temp)
      grad.sytt.mcmc <- coda::as.mcmc(grad.sytt.temp)
      curv.sxxtt.mcmc <- coda::as.mcmc(curv.sxxtt.temp)
      curv.sxytt.mcmc <- coda::as.mcmc(curv.sxytt.temp)
      curv.syytt.mcmc <- coda::as.mcmc(curv.syytt.temp)

      grad.mcmc = list(grad.sx.mcmc = grad.sx.mcmc,
                       grad.sy.mcmc = grad.sy.mcmc,
                       curv.sxx.mcmc = curv.sxx.mcmc,
                       curv.sxy.mcmc = curv.sxy.mcmc,
                       curv.syy.mcmc = curv.syy.mcmc,
                       grad.t.mcmc = grad.t.mcmc,
                       grad.sxt.mcmc = grad.sxt.mcmc,
                       grad.syt.mcmc = grad.syt.mcmc,
                       curv.sxxt.mcmc = curv.sxxt.mcmc,
                       curv.sxyt.mcmc = curv.sxyt.mcmc,
                       curv.syyt.mcmc = curv.syyt.mcmc,
                       grad.tt.mcmc = grad.tt.mcmc,
                       grad.sxtt.mcmc = grad.sxtt.mcmc,
                       grad.sytt.mcmc = grad.sytt.mcmc,
                       curv.sxxtt.mcmc = curv.sxxtt.mcmc,
                       curv.sxytt.mcmc = curv.sxytt.mcmc,
                       curv.syytt.mcmc = curv.syytt.mcmc)

      grad.sx.est <- apply(grad.sx.mcmc, 2, median)
      grad.sy.est <- apply(grad.sy.mcmc, 2, median)
      curv.sxx.est <- apply(curv.sxx.mcmc, 2, median)
      curv.sxy.est <- apply(curv.sxy.mcmc, 2, median)
      curv.syy.est <- apply(curv.syy.mcmc, 2, median)
      grad.t.est <- apply(grad.t.mcmc, 2, median)
      grad.sxt.est <- apply(grad.sxt.mcmc, 2, median)
      grad.syt.est <- apply(grad.syt.mcmc, 2, median)
      curv.sxxt.est <- apply(curv.sxxt.mcmc, 2, median)
      curv.sxyt.est <- apply(curv.sxyt.mcmc, 2, median)
      curv.syyt.est <- apply(curv.syyt.mcmc, 2, median)
      grad.tt.est <- apply(grad.tt.mcmc, 2, median)
      grad.sxtt.est <- apply(grad.sxtt.mcmc, 2, median)
      grad.sytt.est <- apply(grad.sytt.mcmc, 2, median)
      curv.sxxtt.est <- apply(curv.sxxtt.mcmc, 2, median)
      curv.sxytt.est <- apply(curv.sxytt.mcmc, 2, median)
      curv.syytt.est <- apply(curv.syytt.mcmc, 2, median)

      grad.sx.hpd <- data.frame(round(cbind(true.grad.sx[[x]], grad.sx.est, coda::HPDinterval(grad.sx.mcmc)), 2))
      colnames(grad.sx.hpd) <- c("true","est", "lower.ci.95", "upper.ci.95")
      grad.sx.hpd$signif <- apply(grad.sx.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.sx.hpd$cp <- apply(grad.sx.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      grad.sy.hpd <- data.frame(round(cbind(true.grad.sy[[x]], grad.sy.est, coda::HPDinterval(grad.sy.mcmc)), 2))
      colnames(grad.sy.hpd) <- c("true","est", "lower.ci.95", "upper.ci.95")
      grad.sy.hpd$signif <- apply(grad.sy.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.sy.hpd$cp <- apply(grad.sy.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      curv.sxx.hpd <- data.frame(round(cbind(true.curv.sxx[[x]], curv.sxx.est, coda::HPDinterval(curv.sxx.mcmc)), 2))
      colnames(curv.sxx.hpd) <- c("true","est", "lower.ci.95", "upper.ci.95")
      curv.sxx.hpd$signif <- apply(curv.sxx.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxx.hpd$cp <- apply(curv.sxx.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      curv.sxy.hpd <- data.frame(round(cbind(true.curv.sxy[[x]], curv.sxy.est, coda::HPDinterval(curv.sxy.mcmc)), 2))
      colnames(curv.sxy.hpd) <- c("true","est", "lower.ci.95", "upper.ci.95")
      curv.sxy.hpd$signif <- apply(curv.sxy.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxy.hpd$cp <- apply(curv.sxy.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      curv.syy.hpd <- data.frame(round(cbind(true.curv.syy[[x]], curv.syy.est, coda::HPDinterval(curv.syy.mcmc)), 2))
      colnames(curv.syy.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.syy.hpd$signif <- apply(curv.syy.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.syy.hpd$cp <- apply(curv.syy.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.t.hpd <- data.frame(round(cbind(true.grad.t[[x]], grad.t.est, coda::HPDinterval(grad.t.mcmc)), 2))
      colnames(grad.t.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.t.hpd$signif <- apply(grad.t.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.t.hpd$cp <- apply(grad.t.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.sxt.hpd <- data.frame(round(cbind(true.grad.sxt[[x]], grad.sxt.est, coda::HPDinterval(grad.sxt.mcmc)), 2))
      colnames(grad.sxt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.sxt.hpd$signif <- apply(grad.sxt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.sxt.hpd$cp <- apply(grad.sxt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.syt.hpd <- data.frame(round(cbind(true.grad.syt[[x]],grad.syt.est, coda::HPDinterval(grad.syt.mcmc)), 2))
      colnames(grad.syt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.syt.hpd$signif <- apply(grad.syt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.syt.hpd$cp <- apply(grad.syt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      curv.sxxt.hpd <- data.frame(round(cbind(true.curv.sxxt[[x]], curv.sxxt.est, coda::HPDinterval(curv.sxxt.mcmc)), 2))
      colnames(curv.sxxt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.sxxt.hpd$signif <- apply(curv.sxxt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxxt.hpd$cp <- apply(curv.sxxt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      curv.sxyt.hpd <- data.frame(round(cbind(true.curv.sxyt[[x]], curv.sxyt.est, coda::HPDinterval(curv.sxyt.mcmc)), 2))
      colnames(curv.sxyt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.sxyt.hpd$signif <- apply(curv.sxyt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxyt.hpd$cp <- apply(curv.sxyt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      curv.syyt.hpd <- data.frame(round(cbind(true.curv.syyt[[x]], curv.syyt.est, coda::HPDinterval(curv.syyt.mcmc)), 2))
      colnames(curv.syyt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.syyt.hpd$signif <- apply(curv.syyt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.syyt.hpd$cp <- apply(curv.syyt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.tt.hpd <- data.frame(round(cbind(true.grad.tt[[x]], grad.tt.est, coda::HPDinterval(grad.tt.mcmc)), 2))
      colnames(grad.tt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.tt.hpd$signif <- apply(grad.tt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.tt.hpd$cp <- apply(grad.tt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.sxtt.hpd <- data.frame(round(cbind(true.grad.sxtt[[x]], grad.sxtt.est, coda::HPDinterval(grad.sxtt.mcmc)), 2))
      colnames(grad.sxtt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.sxtt.hpd$signif <- apply(grad.sxtt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.sxtt.hpd$cp <- apply(grad.sxtt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      grad.sytt.hpd <- data.frame(round(cbind(true.grad.sytt[[x]], grad.sytt.est, coda::HPDinterval(grad.sytt.mcmc)), 2))
      colnames(grad.sytt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      grad.sytt.hpd$signif <- apply(grad.sytt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      grad.sytt.hpd$cp <- apply(grad.sytt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      curv.sxxtt.hpd <- data.frame(round(cbind(true.curv.sxxtt[[x]], curv.sxxtt.est, coda::HPDinterval(curv.sxxtt.mcmc)), 2))
      colnames(curv.sxxtt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.sxxtt.hpd$signif <- apply(curv.sxxtt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxxtt.hpd$cp <- apply(curv.sxxtt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      curv.sxytt.hpd <- data.frame(round(cbind(true.curv.sxytt[[x]], curv.sxytt.est, coda::HPDinterval(curv.sxytt.mcmc)), 2))
      colnames(curv.sxytt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.sxytt.hpd$signif <- apply(curv.sxytt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.sxytt.hpd$cp <- apply(curv.sxytt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      curv.syytt.hpd <- data.frame(round(cbind(true.curv.syytt[[x]], curv.syytt.est, coda::HPDinterval(curv.syytt.mcmc)), 2))
      colnames(curv.syytt.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      curv.syytt.hpd$signif <- apply(curv.syytt.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      curv.syytt.hpd$cp <- apply(curv.syytt.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })


      t.seq <- seq(x, Ns * Nt, by = Nt)
      y <- model$y[t.seq]
      ci.y.z <- t(apply(model$z[id.mcmc,t.seq], 2, quantile, probs = c(0.5, 0.025, 0.975)))
      y.hpd <- data.frame(round(cbind(y, ci.y.z), 2))
      colnames(y.hpd) <- c("true", "est", "lower.ci.95", "upper.ci.95")
      y.hpd$signif <- apply(y.hpd, 1, function(x){
        if(x[3]>0 & x[4]>0) return (1)
        if(x[3]<0 & x[4]<0) return (-1)
        else return(0)
      })
      y.hpd$cp <- apply(y.hpd,1,function(x){
        ifelse(x[1]>=x[3] & x[1]<=x[4], 0, 1)
      })

      if(plots){
        col.grad <- c()
        col.grad = cbind(col.grad, sapply(grad.sx.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.sy.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxx.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxy.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.syy.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.t.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.sxt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.syt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxxt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxyt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.syyt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.tt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.sxtt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(grad.sytt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxxtt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.sxytt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.grad = cbind(col.grad, sapply(curv.syytt.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        }))
        col.y = sapply(y.hpd$signif, function(sig){
          if(sig == 0) return("#FFFFFF")
          if(sig == 1) return("#00FF00")
          if(sig == -1) return("#00FFFF")
        })
        plots.t = spt_gradients_est_plot(coords = coords,
                                         data_frame = y.hpd[, "est"],
                                         grid.points = grid.points,
                                         gradients = data.frame(nabla.sx = grad.sx.hpd[, "est"],
                                                                nabla.sy = grad.sy.hpd[, "est"],
                                                                nabla.sxx = curv.sxx.hpd[, "est"],
                                                                nabla.sxy = curv.sxy.hpd[, "est"],
                                                                nabla.syy = curv.syy.hpd[, "est"],
                                                                nabla.t = grad.t.hpd[, "est"],
                                                                nabla.sxt = grad.sxt.hpd[, "est"],
                                                                nabla.syt = grad.syt.hpd[, "est"],
                                                                nabla.sxxt = curv.sxxt.hpd[, "est"],
                                                                nabla.sxyt = curv.sxyt.hpd[, "est"],
                                                                nabla.syyt = curv.syyt.hpd[, "est"],
                                                                nabla.tt = grad.tt.hpd[, "est"],
                                                                nabla.sxtt = grad.sxtt.hpd[, "est"],
                                                                nabla.sytt = grad.sytt.hpd[, "est"],
                                                                nabla.sxxtt = curv.sxxtt.hpd[, "est"],
                                                                nabla.sxytt = curv.sxytt.hpd[, "est"],
                                                                nabla.syytt = curv.syytt.hpd[, "est"]),
                                         col.y = col.y, col.grad = col.grad,
                                         only.grad.no.curv = FALSE)

        plots = plot_grid(plotlist = plots.t, labels = LETTERS[1:18],
                          label_size = 6, nrow = 6, ncol = 3, hjust = 0.1)

        return(list(CP = round(c(sum(grad.sx.hpd[,"cp"]),
                                 sum(grad.sy.hpd[,"cp"]),
                                 sum(curv.sxx.hpd[,"cp"]),
                                 sum(curv.sxy.hpd[,"cp"]),
                                 sum(curv.syy.hpd[,"cp"]),
                                 sum(grad.t.hpd[,"cp"]),
                                 sum(grad.sxt.hpd[,"cp"]),
                                 sum(grad.syt.hpd[,"cp"]),
                                 sum(curv.sxxt.hpd[,"cp"]),
                                 sum(curv.sxyt.hpd[,"cp"]),
                                 sum(curv.syyt.hpd[,"cp"]),
                                 sum(grad.tt.hpd[,"cp"]),
                                 sum(grad.sxtt.hpd[,"cp"]),
                                 sum(grad.sytt.hpd[,"cp"]),
                                 sum(curv.sxxtt.hpd[,"cp"]),
                                 sum(curv.sxytt.hpd[,"cp"]),
                                 sum(curv.syytt.hpd[,"cp"]))/nrow(grid.points),2),
                    plots = plots,
                    grad.mcmc = grad.mcmc,
                    y = y.hpd,
                    sx = grad.sx.hpd,
                    sy = grad.sy.hpd,
                    sxx = curv.sxx.hpd,
                    sxy = curv.sxy.hpd,
                    syy = curv.syy.hpd,
                    t = grad.t.hpd,
                    sxt = grad.sxt.hpd,
                    syt = grad.syt.hpd,
                    sxxt = curv.sxxt.hpd,
                    sxyt = curv.sxyt.hpd,
                    syyt = curv.syyt.hpd,
                    tt = grad.tt.hpd,
                    sxtt = grad.sxtt.hpd,
                    sytt = grad.sytt.hpd,
                    sxxtt = curv.sxxtt.hpd,
                    sxytt = curv.sxytt.hpd,
                    syytt = curv.syytt.hpd))
      }else{
        return(list(CP = round(c(sum(grad.sx.hpd[,"cp"]),
                                 sum(grad.sy.hpd[,"cp"]),
                                 sum(curv.sxx.hpd[,"cp"]),
                                 sum(curv.sxy.hpd[,"cp"]),
                                 sum(curv.syy.hpd[,"cp"]),
                                 sum(grad.t.hpd[,"cp"]),
                                 sum(grad.sxt.hpd[,"cp"]),
                                 sum(grad.syt.hpd[,"cp"]),
                                 sum(curv.sxxt.hpd[,"cp"]),
                                 sum(curv.sxyt.hpd[,"cp"]),
                                 sum(curv.syyt.hpd[,"cp"]),
                                 sum(grad.tt.hpd[,"cp"]),
                                 sum(grad.sxtt.hpd[,"cp"]),
                                 sum(grad.sytt.hpd[,"cp"]),
                                 sum(curv.sxxtt.hpd[,"cp"]),
                                 sum(curv.sxytt.hpd[,"cp"]),
                                 sum(curv.syytt.hpd[,"cp"]))/nrow(grid.points),2),
                    grad.mcmc = grad.mcmc,
                    y = y.hpd,
                    sx = grad.sx.hpd,
                    sy = grad.sy.hpd,
                    sxx = curv.sxx.hpd,
                    sxy = curv.sxy.hpd,
                    syy = curv.syy.hpd,
                    t = grad.t.hpd,
                    sxt = grad.sxt.hpd,
                    syt = grad.syt.hpd,
                    sxxt = curv.sxxt.hpd,
                    sxyt = curv.sxyt.hpd,
                    syyt = curv.syyt.hpd,
                    tt = grad.tt.hpd,
                    sxtt = grad.sxtt.hpd,
                    sytt = grad.sytt.hpd,
                    sxxtt = curv.sxxtt.hpd,
                    sxytt = curv.sxytt.hpd,
                    syytt = curv.syytt.hpd))
      }
    })

    ###############
    # Diagnostics #
    ###############
    mspe = c()
    if(plots) d.ggplot = list()
    for(i in 1:Nt){
      spt.grad.obj = plot.fn.cp[[i]]
      if(plots){
        obvf_grad = data.frame(rbind(cbind(spt.grad.obj$y, type = 0),
                                     cbind(spt.grad.obj$sx, type = 1),
                                     cbind(spt.grad.obj$sy, type = 2),
                                     cbind(spt.grad.obj$sxx, type = 3),
                                     cbind(spt.grad.obj$sxy, type = 4),
                                     cbind(spt.grad.obj$syy, type = 5),
                                     cbind(spt.grad.obj$t, type = 6),
                                     cbind(spt.grad.obj$sxt, type = 7),
                                     cbind(spt.grad.obj$syt, type = 8),
                                     cbind(spt.grad.obj$sxxt, type = 9),
                                     cbind(spt.grad.obj$sxyt, type = 10),
                                     cbind(spt.grad.obj$syyt, type = 11),
                                     cbind(spt.grad.obj$tt, type = 12),
                                     cbind(spt.grad.obj$sxtt, type = 13),
                                     cbind(spt.grad.obj$sytt, type = 14),
                                     cbind(spt.grad.obj$sxxtt, type = 15),
                                     cbind(spt.grad.obj$sxytt, type = 16),
                                     cbind(spt.grad.obj$syytt, type = 17)))
        obvf_grad$type = as.factor(obvf_grad$type)
        levels(obvf_grad$type) = c(TeX("$y$"), TeX("$\\partial_{s_x}$"), TeX("$\\partial_{s_y}$"), TeX("$\\partial^2_{s_x}$"), TeX("$\\partial^2_{s_xs_y}$"), TeX("$\\partial^2_{s_y}$"),
                                   TeX("$\\partial_t$"), TeX("$\\partial_t\\partial_{s_x}$"), TeX("$\\partial_t\\nabla_{s_y}$"), TeX("$\\partial_t\\partial^2_{s_x}$"), TeX("$\\partial_t\\partial^2_{s_xs_y}$"), TeX("$\\partial_t\\partial^2_{s_y}$"),
                                   TeX("$\\partial^2_t$"), TeX("$\\partial^2_t\\partial_{s_x}$"), TeX("$\\partial^2_t\\partial_{s_y}$"), TeX("$\\partial^2_t\\partial^2_{s_x}$"), TeX("$\\partial^2_t\\partial^2_{s_xs_y}$"), TeX("$\\partial^2_t\\partial^2_{s_y}$"))
        d.ggplot[[i]] <- ggplot(obvf_grad, aes_string(x = "est", y = "true")) +
          geom_point(size = 0.8) +
          labs(x = "estimated", y = "true") +
          theme_bw() +
          geom_ribbon(aes_string(ymin = "lower.ci.95", ymax = "upper.ci.95"), fill = "blue", alpha = 0.3) +
          geom_line(aes_string(y = "est"), col = "red") +
          facet_wrap(~type, ncol =  3, nrow = 6, scales = "free", labeller = label_parsed)
      }

      ########
      # MSPE #
      ########
      mspe = rbind(mspe, c(sqrt(mean((spt.grad.obj$sx[,1] - spt.grad.obj$sx[,2])^2)),
                           sqrt(mean((spt.grad.obj$sy[,1] - spt.grad.obj$sy[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxx[,1] - spt.grad.obj$sxx[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxy[,1] - spt.grad.obj$sxy[,2])^2)),
                           sqrt(mean((spt.grad.obj$syy[,1] - spt.grad.obj$syy[,2])^2)),
                           sqrt(mean((spt.grad.obj$t[,1] - spt.grad.obj$t[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxt[,1] - spt.grad.obj$sxt[,2])^2)),
                           sqrt(mean((spt.grad.obj$syt[,1] - spt.grad.obj$syt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxxt[,1] - spt.grad.obj$sxxt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxyt[,1] - spt.grad.obj$sxyt[,2])^2)),
                           sqrt(mean((spt.grad.obj$syyt[,1] - spt.grad.obj$syyt[,2])^2)),
                           sqrt(mean((spt.grad.obj$tt[,1] - spt.grad.obj$tt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxtt[,1] - spt.grad.obj$sxtt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sytt[,1] - spt.grad.obj$sytt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxxtt[,1] - spt.grad.obj$sxxtt[,2])^2)),
                           sqrt(mean((spt.grad.obj$sxytt[,1] - spt.grad.obj$sxytt[,2])^2)),
                           sqrt(mean((spt.grad.obj$syytt[,1] - spt.grad.obj$syytt[,2])^2))))
    }
    colnames(mspe) = c("nabs.1", "nabs.2",  "nab2s.11",  "nab2s.12",  "nab2s.22",
                       "nabt", "nab2st.1", "nab2st.2", "nab3st.11", "nab3st.12", "nab3st.22",
                       "nab2t",  "nab3s2t.1", "nab3s2t.2", "nab4s2t.11", "nab4s2t.12", "nab4s2t.22")

    if(plots){
      return(list(plot.fn.cp = plot.fn.cp,
                  mspe = mspe,
                  d.ggplot = d.ggplot))
    }else{
      return(list(plot.fn.cp = plot.fn.cp,
                  mspe = mspe))
    }
  }
}
