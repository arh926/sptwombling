#' Plots estimated spatiotempoal differential processes
#'
#' Uses ggplot2 to generate estimated plots for spatotemporal divergence and Laplacian
#' @param coords spatial coordinates
#' @param Nt number of time points
#' @param data_frame response values: supply as a matrix of ncol = Nt
#' @import MBA ggplot2 cowplot metR
#' @keywords
#' @export
#' @examples
spt_plot <- function(coords = NULL,
                     Nt = NULL,
                     data_frame = NULL, # supplied as a matrix of ncol = Nt
                     shape = NULL){
  if(is.null(shape)){
    par(mfcol = c(1, 1))
    df.gg = list()
    for(i.t in 1:Nt){
      surf = mba.surf(cbind(coords, data_frame[, i.t]),
                      no.X = 100, no.Y = 100, h = 5, m = 2, extend = TRUE)$xyz.est
      gg.grid = expand.grid(surf$x, surf$y)
      colnames(gg.grid) = c("s.x", "s.y")
      df.gg[[i.t]] = cbind(gg.grid, z = as.vector(t(surf$z)), t = i.t)

    }
    df.gg.full = data.frame(do.call(rbind, df.gg))
    df.gg.full$t = as.factor(df.gg.full$t)
    levels(df.gg.full$t) = paste("t = ", 1:Nt, sep = "")

    plot.gg = ggplot(df.gg.full, aes(x = s.y, y = s.x)) +
      theme_cowplot(12) +
      geom_raster(aes(fill =  z)) +
      labs(x = "x", y = "y", fill = "") +
      # other color scaels
      # scico::scale_fill_scico(palette = "vik", label = function(x) sprintf("%.2f", x)) +
      # scale_fill_viridis_c(option = "magma", label = function(x) sprintf("%.2f", x)) +
      scale_fill_distiller(palette = "PiYG", label = function(x) sprintf("%.2f", x)) +
      geom_contour2(aes(z = z, label = after_stat(level))) +
      geom_point(data = data.frame(coords),
                 aes(x = x, y = y),
                 size = 0.8) +
      facet_wrap(~t, ncol = 3) +
      theme(legend.key.height = unit((Nt/3) * 1.3, 'cm'),
            legend.key.width = unit(0.6, 'cm'),
            legend.text = element_text(size = 20),
            strip.text.x = element_text(size = 20))
  }
  return(plot.gg)
}
