#' Plots significance at the triangular region level
#'
#' @param curves_part partitioned curves
#' @param tr.points number of points in triangle (should be fixed at 10 for now)
#' @param length.p length of partition (should be fixed at 10 for now)
#' @param wm.signif wombling measure estimated at the traingular level (should be a data.frame)
#' @keywords plot.surf.wm.sig
#' @import coda
#' @export
plot_surf_tri_sig <- function(curves_part = NULL,
                              tr.points = 10, length.p = 10,
                              wm.signif = NULL){
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
      a.1 = points.sub[1,]
      # (u, 0)
      u.1 = points.sub[2,] - a.1
      # (v, 1)
      v.1 = points.sub[3,] - a.1

      # normal to plane 1
      n.1 = vcrossprod(u = v.1, v = u.1)

      a.2 = points.sub[4,]
      # (u - w, -1)
      u.2 = points.sub[2,] - a.2
      # (v - w, 0)
      v.2 = points.sub[3,] - a.2

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

      points3d(x = surf.pts.1[,1],
               y = surf.pts.1[,2],
               z = surf.pts.1[,3], add = TRUE, col = wm.signif[paste0("t1-",i,"-",j), "signif"])
      points3d(x = surf.pts.2[,1],
               y = surf.pts.2[,2],
               z = surf.pts.2[,3], add = TRUE, col = wm.signif[paste0("t2-",i,"-",j), "signif"])
    }
  }
  return(0)
}
