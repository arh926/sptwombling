#' Triangulates the wombling surface
#'
#' @param curves_part a list of partitioned curves
#' @param length.p length of partition for planar partitioned curves. Defaults to 10.
#' @param tr.points length of partition for triangular points. Defaults to 10
#' @param plot.it logical for whether to plot 3d points
#' @import rgl
#' @keywords surf_triangulate
#' @export
surf_triangulate <- function(curves_part = NULL,
                             tr.points = 10, length.p = 10,
                             plot.it = TRUE){
  nt = unique(curves_part[,"t"])[length(unique(curves_part[,"t"]))]
  om_values <- seq(0, 1, length.out = tr.points) # omega
  up_values <- seq(0, 1, length.out = tr.points) # upsilon

  # Create a grid using expand.grid
  grid <- expand.grid(om = om_values, up = up_values)

  # Filter grid based on the condition u <= 1 - v
  filtered_grid <- grid[grid$om <= 1 - grid$up, ]

  surf.triangle = list()
  normals.triangle = list()

  for(i in 1:(nt - 1)){
    points = rbind(curves_part[(length.p * (i - 1) + 1):(length.p * (i - 1) + length.p),],
                   curves_part[(length.p * i + 1):((i + 1)  * length.p),])

    normals.triangle.time = surf.triangle.time = c()
    for(j in 1:(length.p - 1)){
      points.sub = rbind(points[j:(j + 1),],
                         points[(length.p + j):(length.p + j + 1),])

      # s_0
      a.1 = as.numeric(points.sub[1,])
      # (u, 0)
      u.1 = as.numeric(points.sub[2,] - a.1)
      # (v, 1)
      v.1 = as.numeric(points.sub[3,] - a.1)

      # normal to plane 1
      n.1 = vcrossprod(u = v.1, v = u.1)

      a.2 = as.numeric(points.sub[4,])
      # (u - w, -1)
      u.2 = as.numeric(points.sub[2,] - a.2)
      # (v - w, 0)
      v.2 = as.numeric(points.sub[3,] - a.2)

      # normal to plane 2
      n.2 = vcrossprod(u = v.2, v = -u.2)

      surf.pts.1 = t(apply(filtered_grid, 1, function(x){
        a.1 + x[1] * u.1 + x[2] * v.1
      }))

      surf.pts.2 = t(apply(filtered_grid, 1, function(x){
        a.2 + x[1] * v.2 + x[2] * u.2
      }))

      surf.triangle.time = rbind(surf.triangle.time,
                                 surf.pts.1,
                                 surf.pts.2)
      normals.triangle.time = rbind(normals.triangle.time,
                                    n.1, n.2)

      # text3d(x = points[,1], y = points[,2], z = points[,3],
      #        texts = c("a.t1", "b.t1", "a.t2", "b.t2"), add = T)
      if(plot.it){
        points3d(x = surf.pts.1[,1],
                 y = surf.pts.1[,2],
                 z = surf.pts.1[,3], add = T, col = rep("yellow", nrow(surf.pts.1)))
        points3d(x = surf.pts.2[,1],
                 y = surf.pts.2[,2],
                 z = surf.pts.2[,3], add = T, col = rep("darkgreen", nrow(surf.pts.2)))
      }
    }
    normals.triangle[[i]] = normals.triangle.time
    surf.triangle[[i]] = surf.triangle.time
  }
  normals.triangle = lapply(normals.triangle, function(x){
    x = cbind(x, norm = sqrt(x[,1]^2 + x[,2]^2 + x[,3]^2))
  })

  area.mid = list()

  for(i.1 in 1:(nt-1)){
    grid.box = surf.triangle[[i.1]]
    seq.rect = seq(1, nrow(grid.box), by = tr.points * (tr.points + 1))

    area.mid.df = c()
    tri.df = c()
    for(i.2 in 1:length(seq.rect)){
      id.tr.low = seq.rect[i.2] + 0:54

      frame.low = cbind(c(id.tr.low[1:10]),
                        c(c(id.tr.low[11:19]), NA),
                        c(c(id.tr.low[20:27]), rep(NA, 2)),
                        c(c(id.tr.low[28:34]), rep(NA, 3)),
                        c(c(id.tr.low[35:40]), rep(NA, 4)),
                        c(c(id.tr.low[41:45]), rep(NA, 5)),
                        c(c(id.tr.low[46:49]), rep(NA, 6)),
                        c(c(id.tr.low[50:52]), rep(NA, 7)),
                        c(c(id.tr.low[53:54]), rep(NA, 8)),
                        c(c(id.tr.low[55]), rep(NA, 9)))


      frame.high = t(frame.low) + 55

      for(j.1 in 1:(tr.points - 1)){
        for(j.2 in 1:(tr.points - j.1)){
          d.1 = diff(grid.box[frame.low[j.2, j.1]:frame.low[(j.2 + 1), j.1],])
          d.2 = diff(grid.box[c(frame.low[j.2, j.1], frame.low[j.2, (j.1 + 1)]),])
          n.d1 = sqrt(sum(d.1^2)); n.d2 = sqrt(sum(d.2^2))
          # theta = acos((d.1 %*% t(d.2))/ n.d1/ n.d2)
          if(j.1 > 1){
            if(is.na(frame.low[j.2 + 2, j.1])) area = n.d1 * n.d2 # * sin(theta)/2 # area of parallelogram (not rectangle)
            else area = n.d1 * n.d2 # * sin(theta)
          }else area = n.d1 * n.d2 # * sin(theta)
          midpt = grid.box[frame.low[j.2, j.1],] + d.1/2 + d.2/2
          area.mid.df = rbind(area.mid.df,
                              c(midpt, area, normals.triangle[[i.1]][2 * (i.2 - 1) + 1,]))
          tri.df = c(tri.df, paste("t1", i.1, i.2, sep = "-"))
        }
        colnames(area.mid.df) = c("m.x", "m.y", "m.t", "area", "n.x", "n.y", "n.t", "norm")
      }

      for(j.1 in 1:(tr.points - 1)){
        for(j.2 in 1:(tr.points - j.1)){
          d.1 = diff(grid.box[c(frame.high[j.2, j.1], frame.high[(j.2 + 1), j.1]),])
          d.2 = diff(grid.box[c(frame.high[j.2, j.1]:frame.high[j.2, (j.1 + 1)]),])
          n.d1 = sqrt(sum(d.1^2)); n.d2 = sqrt(sum(d.2^2))
          # theta = acos((d.1 %*% t(d.2))/ n.d1/ n.d2)
          if(j.1 > 1){
            if(is.na(frame.low[j.2 + 2, j.1])) area = n.d1 * n.d2 # * sin(theta)/2 # area of parallelogram (not rectangle)
            else area = n.d1 * n.d2 # * sin(theta)
          }else area = n.d1 * n.d2 # * sin(theta)
          midpt = grid.box[frame.high[j.2, j.1],] + d.1/2 + d.2/2
          area.mid.df = rbind(area.mid.df,
                              c(midpt, area, normals.triangle[[i.1]][2 * (i.2 - 1) + 2,]))
          tri.df = c(tri.df, paste("t2", i.1, i.2, sep = "-"))
        }
        colnames(area.mid.df) = c("m.x", "m.y", "m.t", "area", "n.x", "n.y", "n.t", "norm")
      }
    }
    area.mid[[i.1]] = cbind(data.frame(area.mid.df), tri.id = tri.df)
  }

  return(list(surf.triangle = surf.triangle,
              normals.triangle = normals.triangle,
              wombling.df = area.mid))
}
