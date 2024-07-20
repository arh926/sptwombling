#' Plots estimated differential geometric quantities
#'
#' Uses ggplot2 to generate estimated plots for spatotemporal divergence and Laplacian
#' @param coords coordinates of grid points
#' @param diff.qty differential geometric quantity
#' @param col.vec color vector for significance
#' @import MBA ggplot2 cowplot metR
#' @keywords plot_diff.geo
#' @export
plot_diff.geo <- function(coords = NULL,
                          diff.qty = NULL,
                          col.vec = NULL){
  key.height = 0.6
  key.width = 0.2
  text.size = 8
  contour.text.size = 2.5

  surf = mba.surf(cbind(coords, diff.qty),
                  no.X = 100, no.Y = 100, h = 5, m = 2, extend = TRUE)$xyz.est
  gg.grid = expand.grid(surf$x, surf$y)
  colnames(gg.grid) = c("s.x", "s.y")
  df.p = cbind(gg.grid, z = as.vector(t(surf$z)))

  if(is.null(col.vec)){
    plot.f = ggplot(df.p, aes(x = s.y, y = s.x)) +
      theme_cowplot(12) +
      geom_raster(aes(fill =  z)) +
      labs(x = "", y = "y", fill = "") +
      scale_fill_distiller(palette = "PuOr",
                           label = function(x) sprintf("%.2f", x)) +
      geom_contour2(aes(z = z, label = after_stat(level)), size = 0.1, label_size = contour.text.size) +
      geom_point(data = data.frame(coords),
                 aes(x = x, y = y),
                 size = 1.5) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line = element_line(linewidth = 0.15),
            axis.title = element_text(size = 5),
            axis.text = element_text(size = 5),
            legend.key.height = unit(key.height, 'cm'),
            legend.key.width = unit(key.width, 'cm'),
            legend.text = element_text(size = text.size),
            plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
  }else{
    plot.f = ggplot(df.p, aes(x = s.y, y = s.x)) +
      theme_cowplot(12) +
      geom_raster(aes(fill =  z)) +
      labs(x = "", y = "y", fill = "") +
      scale_fill_distiller(palette = "PuOr",
                           label = function(x) sprintf("%.2f", x)) +
      geom_contour2(aes(z = z, label = after_stat(level)), size = 0.1, label_size = contour.text.size) +
      geom_point(data = data.frame(grid.points),
                 aes(x = x, y = y),
                 fill = col.vec,
                 size = 1.5, stroke = 0.5, pch = 21) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line = element_line(linewidth = 0.15),
            axis.title = element_text(size = 5),
            axis.text = element_text(size = 5),
            legend.key.height = unit(key.height, 'cm'),
            legend.key.width = unit(key.width, 'cm'),
            legend.text = element_text(size = text.size),
            plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
  }

  return(plot.f)
}

# plot.list = list()
#
# for(t.id in 1:9){
#   plot.list.t = list()
#   plot.list.t[[1]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc[[t.id]][,1])
#   plot.list.t[[2]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc.t[[t.id]][,1])
#   plot.list.t[[3]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc.tt[[t.id]][,1])
#   plot.list.t[[4]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc[[t.id]][,2])
#   plot.list.t[[5]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc.t[[t.id]][,2])
#   plot.list.t[[6]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.eig.mcmc.tt[[t.id]][,2])
#   plot.list.t[[7]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.gauss.curv[[t.id]])
#   plot.list.t[[8]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.gauss.curv.t[[t.id]])
#   plot.list.t[[9]] = plot_diff.geo(coords = grid.points,
#                                    diff.qty = true.gauss.curv.tt[[t.id]])
#   plot.list.t[[10]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.div[[t.id]])
#   plot.list.t[[11]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.divt[[t.id]])
#   plot.list.t[[12]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.divtt[[t.id]])
#   plot.list.t[[13]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.lap[[t.id]])
#   plot.list.t[[14]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.lapt[[t.id]])
#   plot.list.t[[15]] = plot_diff.geo(coords = grid.points,
#                                     diff.qty = true.laptt[[t.id]])
#
#   plot.list[[t.id]] = plot_grid(plotlist = plot.list.t, labels = LETTERS[1:15],
#                                 label_size = 6, nrow = 5, ncol = 3, hjust = 0.1)
# }
#
#
# pdf("plots/diff-ge-true-1.pdf", width = 7, height = 9)
# plot.list
# dev.off()
#
# plot.list = list()
#
# for(t.id in 1:9){
#   plot.list.t = list()
#   plot.list.t[[1]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.1.hpd[[t.id]][,1],
#                                  col.vec = col.eig.1[[t.id]])
#   plot.list.t[[2]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.t.1.hpd[[t.id]][,1],
#                                  col.vec = col.eig.t.1[[t.id]])
#   plot.list.t[[3]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.tt.1.hpd[[t.id]][,1],
#                                  col.vec = col.eig.tt.1[[t.id]])
#   plot.list.t[[4]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.2.hpd[[t.id]][,1],
#                                  col.vec = col.eig.2[[t.id]])
#   plot.list.t[[5]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.t.2.hpd[[t.id]][,1],
#                                  col.vec = col.eig.t.2[[t.id]])
#   plot.list.t[[6]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = eig.mcmc.tt.2.hpd[[t.id]][,1],
#                                  col.vec = col.eig.tt.2[[t.id]])
#   plot.list.t[[7]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = gauss.curv.hpd[[t.id]][,1],
#                                  col.vec = col.gauss.curv[[t.id]])
#   plot.list.t[[8]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = gauss.curv.t.hpd[[t.id]][,1],
#                                  col.vec = col.gauss.curv.t[[t.id]])
#   plot.list.t[[9]] = plot_diff.geo(coords = grid.points,
#                                  diff.qty = gauss.curv.tt.hpd[[t.id]][,1],
#                                  col.vec = col.gauss.curv.tt[[t.id]])
#   plot.list.t[[10]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = div.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.div[[t.id]])
#   plot.list.t[[11]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = divt.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.div.t[[t.id]])
#   plot.list.t[[12]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = divtt.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.div.tt[[t.id]])
#   plot.list.t[[13]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = lap.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.lap[[t.id]])
#   plot.list.t[[14]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = lapt.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.lap.t[[t.id]])
#   plot.list.t[[15]] = plot_diff.geo(coords = grid.points,
#                                   diff.qty = laptt.mcmc.hpd[[t.id]][,1],
#                                   col.vec = col.lap.tt[[t.id]])
#
#   plot.list[[t.id]] = plot_grid(plotlist = plot.list.t, labels = LETTERS[1:15],
#                                 label_size = 6, nrow = 5, ncol = 3, hjust = 0.1)
# }
# pdf("plots/diff-ge-1.pdf", width = 7, height = 9)
# plot.list
# dev.off()
