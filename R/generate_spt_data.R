#' Generates simulated space-time data
#'
#' @param Ns number of spatial points
#' @param Nt number of temporal points
#' @param pattern.type four different choices: 1, 2, 3 and 4
#' @param tau error std. deviation
#' @param gradients logical computes only gradients if TRUE
#' @param derived.geom logical computes differential geometric quantities like divergence and Laplacian
#' @param seed sets seed
#' @param grid.points grid points for computing true values
#' @keywords generate_spt_data
#' @export
generate_spt_data <- function(Ns = NULL,
                              Nt = NULL,
                              pattern.type = c("1", "2", "3", "4"),
                              tau = 1,
                              gradients = TRUE,
                              derived.geom = TRUE,
                              seed = NULL,
                              grid.points = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(is.null(Ns)) Ns = 100
  if(is.null(Nt)) Nt = 9

  if(gradients & is.null(grid.points)){
    grid.points = expand.grid(x = seq(0, 1, by = 0.1),
                              y = seq(0, 1, by = 0.1))
    colnames(grid.points) = c("x", "y")
  }

  N = Ns * Nt
  t = seq(1, Nt, length.out = Nt)
  coords = matrix(runif(2 * Ns), ncol = 2)

  colnames(coords) = c("x", "y")


  tau = 1

  if(pattern.type == "1"){
    sim.pattern = array(NA, c(Ns, Nt))
    ## create synthetic y
    for(j in 1:Nt){
      sim.pattern[,j] = rnorm(Ns,
                              mean = 10 * (sin(3 * pi * coords[,1]) + cos(3 * pi * coords[,2]) * cos(t[j] * pi/7)),
                              sd = tau)
    }
    y = c()
    for(i in 1:Ns){
      y = c(y, sim.pattern[i, 1:Nt])
    }
  }else if(pattern.type == "2"){
    sim.pattern = array(NA, c(Ns, Nt))
    ## create synthetic y
    for(j in 1:Nt){
      sim.pattern[,j] = rnorm(Ns,
                              mean = 10 * (sin(3 * pi * coords[,1]) * cos(3 * pi * coords[,2]) * cos(t[j] * pi/7)),
                              sd = tau)
    }
    y = c()
    for(i in 1:Ns){
      y = c(y, sim.pattern[i, 1:Nt])
    }
  }else if(pattern.type == "3"){
    sim.pattern = array(NA, c(Ns, Nt))
    ## create synthetic y
    for(j in 1:Nt){
      sim.pattern[,j] = rnorm(Ns,
                              mean = 5 * (sin(3 * pi * coords[,1]) + cos(3 * pi * coords[,2]) * cos(t[j] * pi/7)),
                              sd = tau)
    }
    y = c()
    for(i in 1:Ns){
      y = c(y, sim.pattern[i, 1:Nt])
    }
  }else if(pattern.type == "4"){
    sim.pattern = array(NA, c(Ns, Nt))
    ## create synthetic y
    for(j in 1:Nt){
      sim.pattern[,j] = rnorm(Ns,
                              mean = 5 * (sin(3 * pi * coords[,1]) * cos(3 * pi * coords[,2]) * cos(t[j] * pi/7)),
                              sd = tau)
    }
    y = c()
    for(i in 1:Ns){
      y = c(y, sim.pattern[i, 1:Nt])
    }
  }
  if(gradients){

    grid.points = expand.grid(x = seq(0, 1, by = 0.1),
                              y = seq(0, 1, by = 0.1))

    if(pattern.type == "1"){
      grad.curv.st.list = list()
      if(derived.geom){
        derived.geom.st = list()
        for(i in t){
          nabla.sx = 30 * pi * cos(3 * pi * grid.points[, 1])
          nabla.sy = -30 * pi * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sxx = -90 * pi^2 * sin(3 * pi * grid.points[,1])
          nabla.sxy = rep(0, nrow(grid.points))
          nabla.syy = -90 * pi^2 * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -10 * pi/7 * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.sx =  rep(0, nrow(grid.points))
          nabla.t.sy = 30 * pi^2/7 * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sxx =  rep(0, nrow(grid.points))
          nabla.t.sxy =  rep(0, nrow(grid.points))
          nabla.t.syy = 90 * pi^3/7 * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.tt = -10 * pi^2/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.tt.sx =  rep(0, nrow(grid.points))
          nabla.tt.sy = 30 * pi^3/49 * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.tt.sxx =  rep(0, nrow(grid.points))
          nabla.tt.sxy =  rep(0, nrow(grid.points))
          nabla.tt.syy = 90 * pi^4/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)

          div = nabla.sx + nabla.sy
          lap = nabla.sxx + nabla.syy
          E.mat = sapply(1:length(nabla.sxx), function(x){
            eval = eigen(matrix(c(nabla.sxx[x], nabla.sxy[x],
                                  nabla.sxy[x], nabla.syy[x]),
                                nrow = 2, ncol = 2, byrow = TRUE))
            return(eval$values)
            }, simplify = TRUE)
          eigen.1 = E.mat[1,]
          eigen.2 = E.mat[2,]
          gauss.curv = eigen.1 * eigen.2
          div.t = nabla.t.sx + nabla.t.sy
          lap.t = nabla.t.sxx + nabla.t.syy
          E.mat.t = sapply(1:length(nabla.t.sxx), function(x){
            eigen(matrix(c(nabla.t.sxx[x], nabla.t.sxy[x],
                                   nabla.t.sxy[x], nabla.t.syy[x]),
                                 nrow = 2, ncol = 2, byrow = TRUE))$values}, simplify = TRUE)
          eigen.1.t = E.mat.t[1,]
          eigen.2.t = E.mat.t[2,]
          gauss.curv.t = eigen.1.t * eigen.2.t
          div.tt = nabla.tt.sx + nabla.tt.sy
          lap.tt = nabla.tt.sxx + nabla.tt.syy
          E.mat.tt = sapply(1:length(nabla.tt.sxx), function(x){
            eigen(matrix(c(nabla.tt.sxx[x], nabla.tt.sxy[x],
                                    nabla.tt.sxy[x], nabla.tt.syy[x]),
                                  nrow = 2, ncol = 2, byrow = TRUE))$values}, simplify = TRUE)
          eigen.1.tt = E.mat.tt[1,]
          eigen.2.tt = E.mat.tt[2,]
          gauss.curv.tt = eigen.1.tt * eigen.2.tt


          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy, nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy),
                                          6)
          derived.geom.st[[i]] <- round(cbind(div, lap, eigen.1, eigen.2, gauss.curv,
                                              div.t, lap.t, eigen.1.t, eigen.2.t, gauss.curv.t,
                                              div.tt, lap.tt, eigen.1.tt, eigen.2.tt, gauss.curv.tt),
                                        6)
        }
      }else{
        for(i in t){
          nabla.sx = 30 * pi * cos(3 * pi * grid.points[, 1])
          nabla.sy = -30 * pi * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sxx = -90 * pi^2 * sin(3 * pi * grid.points[,1])
          nabla.sxy = 0
          nabla.syy = -90 * pi^2 * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -10 * pi/7 * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.sx = 0
          nabla.t.sy = 30 * pi^2/7 * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sxx = 0
          nabla.t.sxy = 0
          nabla.t.syy = 90 * pi^3/7 * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.tt = -10 * pi^2/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.tt.sx = 0
          nabla.tt.sy = 30 * pi^3/49 * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.tt.sxx = 0
          nabla.tt.sxy = 0
          nabla.tt.syy = 90 * pi^4/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)

          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy, nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy),
                                          6)
        }
      }
    }

    if(pattern.type == "2"){

      grad.curv.st.list = list()

      if(derived.geom){
        derived.geom.st = list()
        for(i in t){
          nabla.sx = 30 * pi * cos(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sy = -30 * pi * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.sxx = -90 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2])* cos(pi * i/7)
          nabla.sxy = -90 * pi^2 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2])* cos(pi * i/7)
          nabla.syy = -90 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -10 * pi/7 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sx = -30 * pi^2/7 * cos(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sy = 30 * pi^2/7 * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sxx = 90 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.sxy = 90 * pi^3/7 * cos(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.syy = 90 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.tt = -10 * pi^2/49 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sx = -30 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sy = 30 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sxx = 90 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sxy = 90 * pi^4/49 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.syy = 90 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)

          div = nabla.sx + nabla.sy
          lap = nabla.sxx + nabla.syy
          E.mat = sapply(1:length(nabla.sxx), function(x){
            eval = eigen(matrix(c(nabla.sxx[x], nabla.sxy[x],
                                  nabla.sxy[x], nabla.syy[x]),
                                nrow = 2, ncol = 2, byrow = TRUE))
            return(eval$values)
          }, simplify = TRUE)
          eigen.1 = E.mat[1,]
          eigen.2 = E.mat[2,]
          gauss.curv = eigen.1 * eigen.2
          div.t = nabla.t.sx + nabla.t.sy
          lap.t = nabla.t.sxx + nabla.t.syy
          E.mat.t = sapply(1:length(nabla.t.sxx), function(x){
            eigen(matrix(c(nabla.t.sxx[x], nabla.t.sxy[x],
                           nabla.t.sxy[x], nabla.t.syy[x]),
                         nrow = 2, ncol = 2, byrow = TRUE))$values}, simplify = TRUE)
          eigen.1.t = E.mat.t[1,]
          eigen.2.t = E.mat.t[2,]
          gauss.curv.t = eigen.1.t * eigen.2.t
          div.tt = nabla.tt.sx + nabla.tt.sy
          lap.tt = nabla.tt.sxx + nabla.tt.syy
          E.mat.tt = sapply(1:length(nabla.tt.sxx), function(x){
            eigen(matrix(c(nabla.tt.sxx[x], nabla.tt.sxy[x],
                           nabla.tt.sxy[x], nabla.tt.syy[x]),
                         nrow = 2, ncol = 2, byrow = TRUE))$values}, simplify = TRUE)
          eigen.1.tt = E.mat.tt[1,]
          eigen.2.tt = E.mat.tt[2,]
          gauss.curv.tt = eigen.1.tt * eigen.2.tt


          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy, nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy),
                                          6)

          derived.geom.st[[i]] <- round(cbind(div, lap, eigen.1, eigen.2, gauss.curv,
                                              div.t, lap.t, eigen.1.t, eigen.2.t, gauss.curv.t,
                                              div.tt, lap.tt, eigen.1.tt, eigen.2.tt, gauss.curv.tt),
                                        6)
        }
      }else{
        for(i in t){
          nabla.sx = 30 * pi * cos(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sy = -30 * pi * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.sxx = -90 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2])* cos(pi * i/7)
          nabla.sxy = -90 * pi^2 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2])* cos(pi * i/7)
          nabla.syy = -90 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -10 * pi/7 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sx = -30 * pi^2/7 * cos(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sy = 30 * pi^2/7 * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sxx = 90 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.sxy = 90 * pi^3/7 * cos(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.syy = 90 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.tt = -10 * pi^2/49 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sx = -30 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sy = 30 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sxx = 90 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.sxy = 90 * pi^4/49 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.tt.syy = 90 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)

          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy, nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy), 6)
        }
      }
    }
    if(pattern.type == "3"){
      grad.curv.st.list = list()
      if(derived.geom){
        derived.geom.st = list()
        for(i in t){
          nabla.sx = 15 * pi * cos(3 * pi * grid.points[, 1])
          nabla.sx = nabla.sx - mean(nabla.sx)
          nabla.sy = -15 * pi * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sy = nabla.sy - mean(nabla.sy)
          # nabla.sxx = -45 * pi^2 * sin(3 * pi * grid.points[,1])
          # nabla.sxy = 0
          # nabla.syy = -45 * pi^2 * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -5 * pi/7 * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t = nabla.t - mean(nabla.t)
          nabla.t.sx = 0
          nabla.t.sy = 15 * pi^2/7 * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sy = nabla.t.sy - mean(nabla.t.sy)
          # nabla.t.sxx = 0
          # nabla.t.sxy = 0
          # nabla.t.syy = 45 * pi^3/7 * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          # nabla.tt = -5 * pi^2/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          # nabla.tt.sx = 0
          # nabla.tt.sy = 15 * pi^3/49 * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          # nabla.tt.sxx = 0
          # nabla.tt.sxy = 0
          # nabla.tt.syy = 45 * pi^4/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)

          div = nabla.sx + nabla.sy
          # lap = nabla.sxx + nabla.syy
          # E.mat = eigen(matrix(c(nabla.sxx, nabla.sxy,
          #                        nabla.sxy, nabla.syy),
          #                      nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1 = E.mat[1]
          # eigen.2 = E.mat[2]
          # gauss.curv = eigen.1 * eigen.2
          div.t = nabla.t.sx + nabla.t.sy
          # lap.t = nabla.t.sxx + nabla.t.syy
          # E.mat.t = eigen(matrix(c(nabla.t.sxx, nabla.t.sxy,
          #                          nabla.t.sxy, nabla.t.syy),
          #                        nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1.t = E.mat.t[1]
          # eigen.2.t = E.mat.t[2]
          # gauss.curv.t = eigen.1.t * eigen.2.t
          # div.tt = nabla.tt.sx + nabla.tt.sy
          # lap.tt = nabla.tt.sxx + nabla.tt.syy
          # E.mat.tt = eigen(matrix(c(nabla.tt.sxx, nabla.tt.sxy,
          #                           nabla.tt.sxy, nabla.tt.syy),
          #                         nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1.tt = E.mat.tt[1]
          # eigen.2.tt = E.mat.tt[2]
          # gauss.curv.tt = eigen.1.tt * eigen.2.tt


          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, # nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy # nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                # nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy
          ), 6)

          derived.geom.st[[i]] <- round(cbind(div, #lap, eigen.1, eigen.2, gauss.curv,
                                              div.t #, lap.t, eigen.1.t, eigen.2.t, gauss.curv.t,
                                              # div.tt, lap.tt, eigen.1.tt, eigen.2.tt, gauss.curv.tt
                                              ), 6)
        }
      }else{
        for(i in t){
          nabla.sx = 15 * pi * cos(3 * pi * grid.points[, 1])
          nabla.sy = -15 * pi * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          # nabla.sxx = -45 * pi^2 * sin(3 * pi * grid.points[,1])
          # nabla.sxy = 0
          # nabla.syy = -45 * pi^2 * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -5 * pi/7 * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          nabla.t.sx = 0
          nabla.t.sy = 15 * pi^2/7 * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          # nabla.t.sxx = 0
          # nabla.t.sxy = 0
          # nabla.t.syy = 45 * pi^3/7 * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          # nabla.tt = -5 * pi^2/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          # nabla.tt.sx = 0
          # nabla.tt.sy = 15 * pi^3/49 * sin(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          # nabla.tt.sxx = 0
          # nabla.tt.sxy = 0
          # nabla.tt.syy = 45 * pi^4/49 * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)

          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, # nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy # nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                # nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy
          ), 6)
        }
      }
    }
    if(pattern.type == "4"){
      grad.curv.st.list = list()
      if(derived.geom){
        derived.geom.st = list()
        for(i in t){
          nabla.sx = 15 * pi * cos(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sy = -15 * pi * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.sxx = -45 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2])* cos(pi * i/7)
          # nabla.sxy = -45 * pi^2 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2])* cos(pi * i/7)
          # nabla.syy = -45 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -5 * pi/7 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sx = -15 * pi^2/7 * cos(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sy = 15 * pi^2/7 * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          # nabla.t.sxx = 45 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.t.sxy = 45 * pi^3/7 * cos(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.t.syy = 45 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.tt = -5 * pi^2/49 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sx = -15 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sy = 15 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sxx = 45 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sxy = 45 * pi^4/49 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.syy = 45 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)

          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, # nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy # nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                # nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy
          ), 6)
        }
      }else{
        for(i in t){
          nabla.sx = 15 * pi * cos(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[, 2]) * cos(pi * i/7)
          nabla.sy = -15 * pi * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.sxx = -45 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2])* cos(pi * i/7)
          # nabla.sxy = -45 * pi^2 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2])* cos(pi * i/7)
          # nabla.syy = -45 * pi^2 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          nabla.t = -5 * pi/7 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sx = -15 * pi^2/7 * cos(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          nabla.t.sy = 15 * pi^2/7 * sin(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[, 2]) * sin(pi * i/7)
          # nabla.t.sxx = 45 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.t.sxy = 45 * pi^3/7 * cos(3 * pi * grid.points[, 1]) * sin(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.t.syy = 45 * pi^3/7 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * sin(pi * i/7)
          # nabla.tt = -5 * pi^2/49 * sin(3 * pi * grid.points[, 1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sx = -15 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sy = 15 * pi^3/49 * sin(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sxx = 45 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.sxy = 45 * pi^4/49 * cos(3 * pi * grid.points[,1]) * sin(3 * pi * grid.points[,2]) * cos(pi * i/7)
          # nabla.tt.syy = 45 * pi^4/49 * sin(3 * pi * grid.points[,1]) * cos(3 * pi * grid.points[,2]) * cos(pi * i/7)

          div = nabla.sx + nabla.sy
          # lap = nabla.sxx + nabla.syy
          # E.mat = eigen(matrix(c(nabla.sxx, nabla.sxy,
          #                        nabla.sxy, nabla.syy),
          #                      nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1 = E.mat[1]
          # eigen.2 = E.mat[2]
          # gauss.curv = eigen.1 * eigen.2
          div.t = nabla.t.sx + nabla.t.sy
          # lap.t = nabla.t.sxx + nabla.t.syy
          # E.mat.t = eigen(matrix(c(nabla.t.sxx, nabla.t.sxy,
          #                          nabla.t.sxy, nabla.t.syy),
          #                        nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1.t = E.mat.t[1]
          # eigen.2.t = E.mat.t[2]
          # gauss.curv.t = eigen.1.t * eigen.2.t
          # div.tt = nabla.tt.sx + nabla.tt.sy
          # lap.tt = nabla.tt.sxx + nabla.tt.syy
          # E.mat.tt = eigen(matrix(c(nabla.tt.sxx, nabla.tt.sxy,
          #                           nabla.tt.sxy, nabla.tt.syy),
          #                         nrow = 2, ncol = 2, byrow = TRUE))$values
          # eigen.1.tt = E.mat.tt[1]
          # eigen.2.tt = E.mat.tt[2]
          # gauss.curv.tt = eigen.1.tt * eigen.2.tt


          grad.curv.st.list[[i]] <- round(cbind(nabla.sx, nabla.sy, # nabla.sxx, nabla.sxy, nabla.syy,
                                                nabla.t, nabla.t.sx, nabla.t.sy # nabla.t.sxx, nabla.t.sxy, nabla.t.syy,
                                                # nabla.tt, nabla.tt.sx, nabla.tt.sy, nabla.tt.sxx, nabla.tt.sxy, nabla.tt.syy
          ), 6)

          derived.geom.st[[i]] <- round(cbind(div, #lap, eigen.1, eigen.2, gauss.curv,
                                              div.t #, lap.t, eigen.1.t, eigen.2.t, gauss.curv.t,
                                              # div.tt, lap.tt, eigen.1.tt, eigen.2.tt, gauss.curv.tt
          ), 6)
        }
      }
    }
  }

  if(gradients){
    if(derived.geom){
      return(list(coords = coords,
                  t = t,
                  sim.pattern = sim.pattern,
                  y = y,
                  grid.points = grid.points,
                  gradients = grad.curv.st.list,
                  derived.geom = derived.geom.st))
    }else{
      return(list(coords = coords,
                  t = t,
                  sim.pattern = sim.pattern,
                  y = y,
                  grid.points = grid.points,
                  gradients = grad.curv.st.list))
    }
  }else{return(list(coords = coords,
                t = t,
                sim.pattern = sim.pattern,
                y = y))
  }
}
