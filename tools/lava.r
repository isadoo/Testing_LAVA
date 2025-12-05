###########################################################################
#' @title lava: Log ratios of ancestral variances
#'
#' @description Given coancestry matrices for between and within populations and a trait data frame, 
#' this function estimates the log ratio of ancestral variances using a Bayesian mixed-effects model.
#' 
#' @usage log_av(Theta.P, The.M, trait_dataframe)
#'
#' @param Theta.P A square matrix representing the coancestry matrix between populations.
#'
#' @param The.M A square matrix representing the kinship-based relatedness matrix for individuals.
#'
#' @param trait_dataframe A data frame containing individual IDs and trait values. 
#' The first column should be individual IDs, and the second column should be trait values.
#'
#' @return A lava type object containing:
#' \item{posterior_samples}{A list with posterior samples of variance components and residuals.}
#' \item{BRMS_stats}{A list with the median, lower and upper bounds of the 95% credible interval for the variance difference.}
#' \item{log_ratio}{A list with the p-value of the log ratios, the mean of the log ratio of between- and within-population variance, and confidence intervals.}
#' \item{BRMS_model}{A list with the fitted Bayesian model and hypothesis test results.}
#'
#'
#' @details This function standardizes trait data, constructs a Bayesian mixed-effects model 
#' using `brms`, and estimates ancestral variances. The function assumes no crosses between populations 
#' and analyzes one trait at a time.
#'
#' @author Isabela do O \email{isabela.doo@@unil.ch}
#'
#' @references 
#' 
#' - Goudet & Weir (2023)
#' - do O et al (2025)
#'
lava <- function(Theta.P, 
                The.M, 
                trait_dataframe, 
                column_individual = "id", 
                column_trait = "trait", ...) {
  
  #check input types and dimensions ------------------------
  if (!is.matrix(Theta.P) || !is.matrix(The.M)) {
    stop("Theta.P and The.M must be matrices.")
  }
  
  if (!is.data.frame(trait_dataframe) || ncol(trait_dataframe) < 2) {
    stop("trait_dataframe must be a data frame with at least two columns (ID and trait values).")
  }
  #------------------------------------------------------------


  #extract columns by name if provided ------------------------------------
  id_col <- if(is.numeric(column_individual)) column_individual else which(names(trait_dataframe) == column_individual)
  trait_col <- if(is.numeric(column_trait)) column_trait else which(names(trait_dataframe) == column_trait)

  #Identify populations per individual
  population_blocks_df <- counting_blocks_matrix(The.M) #this function counts the number of blocks of non-zero rows in a matrix
  individuals_per_population_F1 <- population_blocks_df$rows
  number_of_blocks <- length(population_blocks_df$block) 
  pop_ids <- rep(1:number_of_blocks, individuals_per_population_F1[1:number_of_blocks]) 
  number_populations <- nrow(Theta.P)

  
  if (number_of_blocks != number_populations) {
    warning(paste0("Mismatch between detected populations based on The.M matrix (", number_of_blocks,") and Theta.P dimensions (", number_populations, ")."))
  }
  
  #Standardize trait data
  Y <- trait_dataframe[,trait_col]
  Y <- Y - mean(Y)
  var_Y <- var(Y)
  Y <- Y / sqrt(var_Y)
  
  #labels
  row.names(Theta.P) <- colnames(Theta.P) <- paste("pop", 1:number_populations, sep = "_")
  row.names(The.M) <- colnames(The.M) <- trait_dataframe[,id_col]
  
  #From VB = VA*2FST
  two.Theta.P <- 2 * Theta.P

  pop_labels <- paste0("pop_", pop_ids)

  dat <- data.frame(pop = pop_labels, ind = trait_dataframe[,id_col], Y = Y)
  
  #Bayesian model - using brms package
  brms_mf <- brm(Y ~ 1 + (1 | gr(pop, cov = two.Theta.P)) + (1 | gr(ind, cov = The.M)), 
                 data = dat, data2 = list(two.Theta.P = two.Theta.P, The.M = The.M), 
                 family = gaussian(), chains = 8, cores = 4, iter = 3000, warmup = 1000, thin = 2, ...)
  
  #variance components
  var_components <- lapply(VarCorr(brms_mf, summary = FALSE), function(x) x$sd^2)
  var_df <- as.data.frame(do.call(cbind, var_components))
  
  quant_med <- quantile(var_df$pop - var_df$ind, c(0.5, 0.025, 0.975))
  mean_diff <- mean(var_df$pop - var_df$ind)
  
  #Hypothesis testing
  hyp <- "sd_pop__Intercept^2 - sd_ind__Intercept^2 = 0"
  the_hyp <- hypothesis(brms_mf, hyp, class = NULL)
  
  #Posteriors: VA,B and VA,A - estimated ancestral variances
  post_samples <- as_draws_df(brms_mf, variable = c("sd_pop__Intercept", "sd_ind__Intercept"))
  post_samples$log_ratio <- log(post_samples$sd_pop__Intercept^2 / post_samples$sd_ind__Intercept^2)
  
  mean_log_ratio <- mean(post_samples$log_ratio)
  quant_log_ratio <- quantile(post_samples$log_ratio, probs = c(0.025, 0.975))
  
  
  p_value <- 2 * mean(sign(post_samples$log_ratio) != sign(median(post_samples$log_ratio)))
  
  #S3 object of class "lava"
  results <- list(
    posteriors_samples = post_samples,
    
    #basic statistics together
    BRMS_stats = list(
      median = quant_med["50%"],
      median_lower = quant_med["2.5%"],
      median_upper = quant_med["97.5%"],
      mean_diff = mean_diff
    ),
    
    #log ratio statistics together
    log_ratio = list(
      p_value = p_value,
      mean_log_ratio = mean_log_ratio,
      log_ratio_ci_lower = quant_log_ratio[1],
      log_ratio_ci_upper = quant_log_ratio[2]
    ),
    
    #model and hypothesis test results together
    BRMS_model = list(
      model = brms_mf,
      hypothesis = the_hyp$hypothesis[2:5]
    ),
    
    trait_name = names(trait_dataframe)[trait_col]
  )
  
  #Setting the class attribute to create an S3 object
  class(results) <- "lava"
  
  return(results)
}

#Defining print method for the lava S3 object
#' @export
#' @method print lava
print.lava <- function(x, ...) {
  #header
  cat("\n===============================\n")
  cat("Log Ancestral Variance Analysis (LAVA)\n")
  cat("===============================\n\n")
  
  #key findings
  cat("Log ratio of estimated ancestral variances (between equation/ within equation):\n")
  cat(sprintf("  Mean: %.4f (95%% CI: %.4f to %.4f)\n", 
              x$log_ratio$mean_log_ratio, 
              x$log_ratio$log_ratio_ci_lower, 
              x$log_ratio$log_ratio_ci_upper))
  
  cat("\nDifference in variance components (between equation/ within equation):\n")
  cat(sprintf("  Mean: %.4f\n", x$BRMS_stats$mean_diff))
  cat(sprintf("  Median: %.4f (95%% CI: %.4f to %.4f)\n", 
              x$BRMS_stats$median, 
              x$BRMS_stats$median_lower, 
              x$BRMS_stats$median_upper))
  
  cat("\nStatistical tests:\n")
  cat(sprintf("  Hypothesis test p-value: %.4f\n", x$log_ratio$p_value))
  
  cat("\n===============================\n")
  
  #~~invisible return~~
  invisible(x)
}


#plotting method for lava S3 object
#' @export
#' @method plot lava
plot.lava <- function(x, which = "both", 
                     main_density = "Posterior Distribution of Log Ratio of Ancestral Variances",
                     main_scatter = "Posterior Samples of Variance Components",
                     ...) {
  
  #check if we print just one of the plots - 'both' is the default
  which <- match.arg(which, choices = c("both", "density", "scatter"))
  
  #setting up the figure layout
  if (which == "both") {
    old_par <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
    on.exit(par(old_par))  # Restore original parameters when function exits
  } else {
    old_par <- par(mar = c(4, 4, 3, 2))
    on.exit(par(old_par))
  }
  
  #plot the density function ------------------------------------
  plot_density <- function() {
    log_ratio <- x$post_samples$log_ratio
    
    
    dens <- density(log_ratio)
    
    
    plot(dens, main = main_density,
         xlab = "Log Ratio (Between-population / Within-population)",
         ylab = "Density", col = "blue", lwd = 2)
    
    #pretty shading
    polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), 
            col = rgb(0, 0, 1, 0.3), border = NA)
    
    #reference line at 0 (equal variances)
    abline(v = 0, lty = 2, col = "red", lwd = 2)
    
    #mean 
    abline(v = x$BRMS_log$mean_log_ratio, col = "darkblue", lwd = 2)
    
    #CI lines
    abline(v = x$BRMS_log$log_ratio_ci_lower, col = "darkblue", lty = 3, lwd = 1.5)
    abline(v = x$BRMS_log$log_ratio_ci_upper, col = "darkblue", lty = 3, lwd = 1.5)
    
    
    legend("topright", 
           legend = c("Density", "Equal variances", "Mean", "95% CI"),
           lty = c(1, 2, 1, 3), 
           col = c("blue", "red", "darkblue", "darkblue"),
           lwd = c(2, 2, 2, 1.5),
           bty = "n")
  }
  
  #scatter function ------------------------------------
  plot_scatter <- function() {
    between_var <- x$post_samples$between_pop_variance
    within_var <- x$post_samples$within_pop_variance
    
    
    xlim <- ylim <- range(c(between_var, within_var))
    
    
    plot(within_var, between_var, 
         main = main_scatter,
         xlab = "Within-population Variance", 
         ylab = "Between-population Variance",
         pch = 16, col = rgb(0, 0, 1, 0.3),
         xlim = xlim, ylim = ylim)
    
    #equality line (1:1)
    abline(0, 1, lty = 2, col = "red", lwd = 2)
    
    #mean point
    points(mean(within_var), mean(between_var), 
           pch = 16, col = "darkred", cex = 1.5)
    
    #legend
    legend("topleft", 
           legend = c("Posterior samples", "Equal variances", "Mean"),
           pch = c(16, NA, 16), 
           lty = c(NA, 2, NA),
           col = c(rgb(0, 0, 1, 0.3), "red", "darkred"),
           lwd = c(NA, 2, NA),
           pt.cex = c(1, NA, 1.5),
           bty = "n")
  }
  
  #plotting only what was asked
  if (which == "density" || which == "both") {
    plot_density()
  }
  
  if (which == "scatter" || which == "both") {
    plot_scatter()
  }
  
  invisible(x)
}