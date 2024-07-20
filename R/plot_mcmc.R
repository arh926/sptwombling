# install.packages("ggmcmc")
# require(ggmcmc)
# require(coda)
# theta = coda::as.mcmc.list(coda::as.mcmc(cbind(beta.0 = results[[1]]$recover.beta[thin_id],
#                                                sigma.2 = results[[1]]$sig2[thin_id],
#                                                tau.2 = results[[1]]$tau2[thin_id],
#                                                phis = results[[1]]$phis[thin_id],
#                                                phit = results[[1]]$phit[thin_id])))
# theta.mcmc = ggs(theta)
# g1 = ggs_traceplot(theta.mcmc) + theme_bw()
# g2 = ggs_compare_partial(theta.mcmc) + theme_bw()
# g3 = ggs_autocorrelation(theta.mcmc) + theme_bw()
#
# cowplot::plot_grid(g1, g2, g3, labels = LETTERS[1:3], ncol = 3)
#
# ggs_pairs(theta.mcmc, lower = list(continuous = "density"))

############################################
# A plotting function for your MCMC output #
############################################
require(latex2exp)
plot_mcmc <- function(samples = NULL, # can be a vector, or a matrix--vector means 1 parameter, matrix means multiple parameters
                      true = NULL, # true values if available
                      col = "darkgreen",
                      acf.lag.max = 100, # no need to change
                      cnames = NULL, # names of parameters
                      diagnostic.measures = TRUE){

  # Detecting number of parameters in your input
  nc = ncol(samples)
  if(is.null(nc)) nc = 1
  N = nrow(samples)
  if(is.null(cnames)) cnames = colnames(samples) # give column names to matrix of samples
  if(nc == 1){
    # Converting samples into MCMC objects
    if(!coda::is.mcmc(samples)) samples = coda::as.mcmc(samples)
  }else{
    # Converting samples into MCMC objects
    if(!sum(apply(samples,2, function(x) coda::is.mcmc(x)))) samples = apply(samples,2, function(x) coda::as.mcmc(x))
  }



  if(nc == 1){
    # par(mfcol = c(1, 3))
    ##############
    # Trace Plot #
    ##############
    ts.plot(samples, main = TeX(paste0("Trace of ", cnames)),
            ylab = "", xlab = "Iterations", col = "blue")
    abline(h = quantile(samples, probs = 0.025), col = "orange", lty = "dashed", lwd = 1.5)
    abline(h = median(samples), lty = "dotted", col = col, lwd = 2)
    abline(h = quantile(samples, probs = 0.975), col = "orange", lty = "dashed", lwd = 1.5)
    if(!is.null(true)) abline(h = true, lty = "dotted", col = "darkred", lwd = 2)

    ################
    # Density Plot #
    ################
    plot(density(samples), main = TeX(paste0("Density of ", cnames)), ylab = "",
         xlab = paste("N =", N, "Bandwidth =", sprintf("%.3f", round(density(samples)$bw, 3)), sep = " "),
         col="blue")
    rug(jitter(samples))
    abline(v = quantile(samples, probs = 0.025), col = "orange", lty = "dashed", lwd = 1.5)
    abline(v = median(samples), lty = "dotted", col = col, lwd = 2)
    abline(v = quantile(samples, probs = 0.975), col = "orange", lty = "dashed", lwd = 1.5)
    if(!is.null(true)) abline(v = true, lty = "dotted", col = "darkred", lwd = 2)

    #######
    # ACF #
    #######
    acf(samples, lag.max = acf.lag.max, main = "") # not using main since margin is wrong
    title(TeX(paste0("ACF for ", cnames)))


  }else{

    if(sum(par()$mfrow==c(6, 3)) < 2){
      par(mfrow=c(6, 3))
      par(mar = rep(2, 4))
    }
    for(i in 1:nc){

      ##############
      # Trace Plot #
      ##############
      ts.plot(samples[,i], main = TeX(paste0("Trace of ", cnames[i])), ylab = "", xlab = "Iterations", col="blue")
      abline(h = quantile(samples[,i], probs = 0.025), col = "orange", lty = "dashed", lwd = 1.5)
      abline(h = median(samples[,i]), lty = "dotted", col = col, lwd = 2)
      abline(h = quantile(samples[,i], probs = 0.975), col = "orange", lty = "dashed", lwd = 1.5)
      if(!is.null(true)) abline(h = true[i], lty="dotted", col = "darkred", lwd = 2)

      ################
      # Density Plot #
      ################
      plot(density(samples[,i]), main = TeX(paste0("Density of ", cnames[i])), ylab = "",
           xlab = paste("N =", N, "Bandwidth =", sprintf("%.3f", round(density(samples[,i])$bw,3)), sep = " "), col="blue")
      rug(jitter(samples[,i]))
      abline(v = quantile(samples[,i], probs = 0.025), col = "orange", lty = "dashed", lwd = 1.5)
      abline(v = median(samples[,i]), lty = "dotted", col = col, lwd = 2)
      abline(v = quantile(samples[,i], probs = 0.975), col = "orange", lty = "dashed", lwd = 1.5)
      if(!is.null(true)) abline(v = true[i], lty="dotted", col = "darkred", lwd = 2)

      #######
      # ACF #
      #######
      acf(samples[,i], lag.max = acf.lag.max, main = "") # not using main since margin is wrong
      title(TeX(paste0("ACF for ", cnames[i])))

    }
  }

  if(diagnostic.measures){
    if(nc == 1){
      cat("Geweke-Brooks Diagnostic::",
          pnorm(abs(coda::geweke.diag(samples)$z), lower.tail = FALSE) * 2,
          "(values should be greater than 0.05 for convergence)")
    }else{
      gd.vals = apply(samples, 2, function(x){
        pnorm(abs(coda::geweke.diag(x)$z), lower.tail = FALSE) * 2
      })
      cat("Geweke-Brooks Diagnostic::", median(gd.vals),"(",
          quantile(gd.vals, probs = 0.025),
          ",",
          quantile(gd.vals, probs = 0.975), ")", "\n",
          "Min.::", min(gd.vals),
          "(value should be greater than 0.05 for convergence)")
    }
  }
  if(nc != 1) return(gd.vals)
}
