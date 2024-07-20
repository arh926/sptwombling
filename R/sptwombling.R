#' Function for performing Bayesian spatiotemporal Wombling on surfaces
#'
#' Performs boundary analysis for surfaces. This uses a two-dimensional Riemann sum approximation.
#'
#' @param model a list containing posterior samples from a collapsed or Gibbs sampler that fits a spatiotemporal Bayesian hierarchical model
#' @param wombling.df a data-frame consisting of points, areas, normals for a triangulated surface
#' @keywords sptwombling
#' @import coda parallel magic
#' @export
sptwombling <- function(model = NULL,
                        wombling.df = NULL){
  grid.points = data.frame(wombling.df[, c(1:3)])
  unit.normals = wombling.df[, c("n.x", "n.y", "n.t")]/wombling.df[,"norm"]

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

  I.d = diag(d)
  I.tilde = matrix(c(3, 0, 1,
                     0, 1, 0,
                     1, 0, 3),
                   nrow = 3, ncol = 3, byrow = TRUE)
  I.star = c(1, 0, 1)
  ngrad = 3 * (1 + d + d * (d + 1)/2) - 1

  stm = proc.time()
  results.grad <- mclapply(parallel.index, function(x){
    samp.x <- samp.list[[x]]

    sig2 <- model$sig2[samp.x]
    phis <- model$phis[samp.x]
    phit <- model$phit[samp.x]
    Z <- model$z[samp.x,]

    n.mcmc <- length(sig2)

    mcmc.grad <- list()
    for(i.mcmc in 1:n.mcmc){

      R.Z <- st_cov_matern2(delta = delta,
                            Delta = Delta,
                            lphis = log(phis[i.mcmc]),
                            lphit = log(phit[i.mcmc]),
                            lsig2 = 0)

      R.Z.in <- try(chol2inv(chol(R.Z)), silent = TRUE)
      if("try-error" %in% class(R.Z.in)) R.Z.in <- chol2inv(chol(R.Z + 1e-4 * diag(N)))

      grad.est <- matrix(NA, nrow = nrow(grid.points), ncol = ngrad)

      for(i in 1:nrow(grid.points)){
        # Data Object:: Matrix of: ||\Delta||, \Delta_1^2, \Delta_1\Delta_2, \Delta_2^2
        dist.s0 <- apply(coords, 1, function(x) sqrt(sum((x - grid.points[i, c(1:2)])^2)))
        delta.s0 <- do.call(rbind, t(apply(coords, 1, function(x) x - grid.points[i, c(1:2)])))
        delta.t0 <- t - grid.points[i, 3]
        dist.delta.mat <- cbind(dist.s0, delta.s0[,1]^2, delta.s0[,1] * delta.s0[,2], delta.s0[,2]^2)

        A.t <- (phit[i.mcmc]^2 * delta.t0^2 + 1)
        ###########################
        # gradients and curvature #
        ###########################
        # \nabla_s
        nabKs <- as.vector(unlist(sapply(dist.s0, function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x/sqrt(A.t))/A.t^2, simplify = FALSE))) * delta.s0[rep(seq_len(Ns), each = Nt),]; rownames(nabKs) <- NULL
        # \nabla_s^2
        nabK2s.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[2]/A.t)/A.t^2, simplify = FALSE)))
        nabK2s.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 25/3 * phis[i.mcmc]^4 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * x[3]/A.t^3, simplify = FALSE)))
        nabK2s.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -5/3 * phis[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[4]/A.t)/A.t^2, simplify = FALSE)))
        # \nabla_t
        nabKt <- as.vector(unlist(sapply(dist.s0, function(x) -2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (1 + sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) + 5/6 * phis[i.mcmc]^2 * x^2/A.t - 5 * sqrt(5)/6 * phis[i.mcmc]^3 * x^3/A.t^(3/2)) * delta.t0/A.t^2, simplify = FALSE)))
        # \nabla_t\nabla_s
        nabKst <- as.vector(unlist(sapply(dist.s0, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 5 * phis[i.mcmc]^2 * x^2/A.t)/A.t^3 * delta.t0, simplify = FALSE))) * delta.s0[rep(seq_len(Ns), each = Nt),]; rownames(nabKst) <- NULL
        # \nabla_t\nabla_s^2
        nabK2st.11 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[2]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[2] * x[1]/A.t^(3/2)) * delta.t0/A.t^3, simplify = FALSE)))
        nabK2st.12 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) -25/3 * phis[i.mcmc]^4 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (6 - sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * x[3] * delta.t0/A.t^4, simplify = FALSE)))
        nabK2st.22 <- as.vector(unlist(apply(dist.delta.mat, 1, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t)) * (4 + 4 * sqrt(5) * phis[i.mcmc] * x[1]/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x[1]^2/A.t - 30 * phis[i.mcmc]^2 * x[4]/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x[4] * x[1]/A.t^(3/2)) * delta.t0/A.t^3, simplify = FALSE)))

        # \nabla_t^2
        nabK2t <- as.vector(unlist(sapply(dist.s0, function(x) -2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * ((1 + sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) + 5/6 * phis[i.mcmc]^2 * x^2/A.t - 5 * sqrt(5)/6 * phis[i.mcmc]^3 * x^3/A.t^(3/2)) -
                                                                                                                                     (4 + 4 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 20 * sqrt(5)/3 * phis[i.mcmc]^3 * x^3/A.t^(3/2) + 25/6 * phis[i.mcmc]^4 * x^4/A.t^2) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^2, simplify = FALSE)))
        # \nabla_t^2\nabla_s
        nabKs2t <- as.vector(unlist(sapply(dist.s0, function(x) 5/3 * phis[i.mcmc]^2 * phit[i.mcmc]^2 * exp(-sqrt(5) * phis[i.mcmc] * x/sqrt(A.t)) * ((4 + 4 * sqrt(5) * phis[i.mcmc] * x/sqrt(A.t) - 5 * phis[i.mcmc]^2 * x^2/A.t) -
                                                                                                                                                        (24 + 24 * sqrt(5) * phis[i.mcmc] * x/A.t^(1/2) - 60 * phis[i.mcmc]^2 * x^2/A.t + 5 * sqrt(5) * phis[i.mcmc]^3 * x^3/A.t^(3/2)) * phit[i.mcmc]^2 * delta.t0^2/A.t)/A.t^3, simplify = FALSE))) * delta.s0[rep(seq_len(Ns), each = Nt),]; rownames(nabKs2t) <- NULL
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
  stp = proc.time() - stm

  results.grad = unlist(results.grad, recursive = FALSE)

  results.dir.grads = lapply(results.grad, function(x){
    dir.grads = matrix(NA, nrow = nrow(x), ncol = 8)
    for(i in 1:nrow(x)){
      # N_ijk,st
      u.mat = matrix(0, nrow = 8, ncol = ncol(x))
      u.mat[1, 1:2] = c(unit.normals[i,1], unit.normals[i,2])
      u.mat[2, 3:5] = c(unit.normals[i,1]^2, 2 * unit.normals[i,1] * unit.normals[i,2], unit.normals[i,2]^2)
      u.mat[3, 6] = unit.normals[i, 3]
      u.mat[4, 7:8] = unit.normals[i, 3] * c(unit.normals[i,1], unit.normals[i,2])
      u.mat[5, 9:11] = unit.normals[i, 3] * c(unit.normals[i,1]^2, 2 * unit.normals[i,1] * unit.normals[i,2], unit.normals[i,2]^2)
      u.mat[6, 12] = unit.normals[i, 3]^2
      u.mat[7, 13:14] = unit.normals[i, 3]^2 * c(unit.normals[i,1], unit.normals[i,2])
      u.mat[8, 15:17] = unit.normals[i, 3]^2 * c(unit.normals[i,1]^2, 2 * unit.normals[i,1] * unit.normals[i,2], unit.normals[i,2]^2)
      dir.grads[i,] = u.mat %*% x[i,]
    }
    dir.grads
  })


  results.riemann.sum.approx = lapply(results.dir.grads, function(x){
    x.area = apply(x, 2, function(a) a * wombling.df$norm * wombling.df$area) # total wombling measure
    # need to calculate area
    riemann.sum = aggregate.data.frame(x.area, by = list(wombling.df$tri.id), "sum")
    order.id = match(unique(wombling.df$tri.id), riemann.sum$Group.1)
    riemann.sum = riemann.sum[order.id,]
    cnames = paste("WM", 1: (Nt-1), sep = "."); rnames  = riemann.sum[,1]
    riemann.sum = riemann.sum[,-1]
    colnames(riemann.sum) = cnames; rownames(riemann.sum) = rnames
    riemann.sum
  })

  nbatch = 100
  slist <- split(id.mcmc, ceiling(seq_along(id.mcmc)/(length(id.mcmc)/nbatch)))

  WM = list()

  for(i in 1:8){
    wm = do.call(rbind, lapply(results.riemann.sum.approx, function(x) x[,i]))
    wm.tmp = t(sapply(slist, function(x) apply(wm[x,], 2, median)))
    wm.mcmc = as.mcmc(wm.tmp)
    wm.est = apply(wm.mcmc, 2, median); names(wm.est) = unique(wombling.df$tri.id)
    wm.hpd = data.frame(cbind(wm.est, HPDinterval(wm.mcmc)))
    colnames(wm.hpd) <- c("est", "lower.ci.95", "upper.ci.95")
    wm.hpd$signif <- apply(wm.hpd,1,function(x){
      if(x[2]>0 & x[3]>0) return (1)
      if(x[2]<0 & x[3]<0) return (-1)
      else return(0)
    })
    WM[[i]] = wm.hpd
  }
  names(WM) = c("s.g","s.c",
                "t.g", "st.gg", "st.cg",
                "t.c", "st.gc", "st.cc")
  return(list(WM = WM, time = stp))
}
