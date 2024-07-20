#' Plots estimated spatiotemporal differential processes
#'
#' Uses `ggplot2` to generate estimated plots for spatiotemporal divergence and Laplacian
#' @param coords spatial coordinates
#' @param data_frame response values: supply only one time point at a time
#' @param grid.points grid points over which gradients/curvature is estimated
#' @param gradients data frame of estimated gradients
#' @param only.grad.no.curv logical for only gradients
#' @param point.size size of points to be plotted
#' @param shape should be NULL for now
#' @import MBA ggplot2 cowplot metR
#' @keywords spt_gradients_plot
#' @export
spt_gradients_plot <- function(coords = NULL,
                               data_frame = NULL, # supply only one time point at a time
                               grid.points = NULL, # supply as "data.frame()"
                               gradients = NULL,
                               shape = NULL,
                               point.size = 1,
                               only.grad.no.curv = FALSE){
  if(only.grad.no.curv & ncol(gradients) > 5) stop("Please provide only spatial and temporal gradients in gradients argument")

  key.height = 0.6
  key.width = 0.2
  text.size = 8
  contour.text.size = 2.5

  # if(only.grad.no.curv){
  #   key.height = 1.2
  #   key.width = 0.3
  #   text.size = 12
  #   contour.text.size = 4
  # }else{
  # }

  if(is.null(shape)){
    plots.time = list()
    # Process
    surf = mba.surf(cbind(coords, data_frame),
                    no.X = 100, no.Y = 100, h = 5, m = 2, extend = TRUE)$xyz.est
    gg.grid = expand.grid(surf$x, surf$y)
    colnames(gg.grid) = c("s.x", "s.y")
    df.p = cbind(gg.grid, z = as.vector(t(surf$z)))

    plots.time[[1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
      theme_cowplot(12) +
      geom_raster(aes_string(fill =  "z")) +
      labs(x = "x", y = "y", fill = "") +
      scale_fill_distiller(palette = "PiYG",
                           label = function(x) sprintf("%.2f", x)) +
      geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
      geom_point(data = data.frame(coords),
                 aes_string(x = "x", y = "y"),
                 size = point.size) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line = element_line(linewidth = 0.1),
            axis.title = element_text(size = 5),
            axis.text = element_text(size = 5),
            legend.key.height = unit(key.height, 'cm'),
            legend.key.width = unit(key.width, 'cm'),
            legend.text = element_text(size = text.size),
            plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
  }

  if(only.grad.no.curv){
    for(i.g in 1:ncol(gradients)){
      surf = mba.surf(cbind(grid.points, gradients[,i.g]),
                      no.X = 100, no.Y = 100, h = 5, m = 2, extend = TRUE)$xyz.est
      gg.grid = expand.grid(surf$x, surf$y)
      colnames(gg.grid) = c("s.x", "s.y")
      df.p = cbind(gg.grid, z = as.vector(t(surf$z)))

      if(i.g <= 2){
        # Spatial Gradients
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "Spectral", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g == 3){
        # Temporal Gradients of process
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "RdYlBu", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g > 3){
        # Temporal Gradients of spatial gradients
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "PRGn", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
    }
  }else{
    for(i.g in 1:ncol(gradients)){
      surf = MBA::mba.surf(cbind(grid.points, gradients[,i.g]),
                           no.X = 100, no.Y = 100, h = 5, m = 2, extend = TRUE)$xyz.est
      gg.grid = expand.grid(surf$x, surf$y)
      colnames(gg.grid) = c("s.x", "s.y")
      df.p = cbind(gg.grid, z = as.vector(t(surf$z)))

      if(i.g <= 5){
        # Spatial Gradients & Curvature
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "Spectral", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g == 6){
        # Temporal Gradients of process
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "RdYlBu", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g > 6 & i.g <= 11){
        # Temporal Gradients of spatial gradients & curvature
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "PRGn", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g == 12){
        # Temporal Curvature of process
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "RdYlBu", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
      if(i.g > 12){
        # Temporal Curvature of spatial gradients & curvature
        plots.time[[i.g + 1]] = ggplot(df.p, aes_string(x = "s.y", y = "s.x")) +
          theme_cowplot(12) +
          geom_raster(aes_string(fill =  "z")) +
          labs(x = "x", y = "y", fill = "") +
          scale_fill_distiller(palette = "BrBG", label = function(x) sprintf("%.2f", x)) +
          geom_contour2(aes_string(z = "z", label = "after_stat(level)"), size = 0.1, label_size = contour.text.size) +
          geom_point(data = data.frame(grid.points),
                     aes_string(x = "x", y = "y"),
                     size = point.size) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.line = element_line(linewidth = 0.1),
                axis.title = element_text(size = 5),
                axis.text = element_text(size = 5),
                legend.key.height = unit(key.height, 'cm'),
                legend.key.width = unit(key.width, 'cm'),
                legend.text = element_text(size = text.size),
                plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))
      }
    }
  }
  plots.time
}
