# Suppress R CMD check NOTEs for ggplot2 aes variables
utils::globalVariables(c("est", "lo", "hi", "label"))
utils::globalVariables(c("value", "covariate", "trt_group"))

#' Diagnostic: Covariate Distribution by Treatment Group
#'
#' Plots covariate distributions stratified by treatment group.
#' Numeric covariates with more than \code{cat_threshold} unique values are
#' shown as overlapping density curves; all other covariates (binary,
#' ordinal, character, factor) are shown as side-by-side bar charts.
#'
#' @param data A data.frame.
#' @param treatment Character. Name of binary treatment column.
#' @param covariates Character vector of covariate column names.
#' @param labels Character vector of length 2. Labels for treatment = 0
#'   and treatment = 1. Default \code{c("Control", "Treated")}.
#' @param cat_threshold Integer. Numeric covariates with \emph{fewer than or
#'   equal to} this many unique values are treated as categorical.
#'   Default \code{5}.
#'
#' @return A named list with elements \code{$continuous} and/or
#'   \code{$categorical}, each a ggplot2 object (or \code{NULL} if there
#'   are no covariates of that type). Both plots are printed as a side effect.
#' @export
plot_covariate_diagnostics <- function(data,
                                       treatment = "A",
                                       covariates,
                                       labels = c("Control", "Treated"),
                                       cat_threshold = 5L) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Install with install.packages('tidyr').")
  }

  data <- as.data.frame(data)

  # Label treatment groups
  trt_vals   <- sort(unique(data[[treatment]]))
  trt_labels <- if (length(trt_vals) == 2 && all(trt_vals %in% c(0, 1, 0L, 1L))) {
    labels
  } else {
    paste0("Group (A=", trt_vals, ")")
  }
  data[[treatment]] <- factor(data[[treatment]],
                              levels = trt_vals, labels = trt_labels)

  # Classify covariates
  is_cont <- sapply(data[, covariates, drop = FALSE], function(x) {
    is.numeric(x) && length(unique(x)) > cat_threshold
  })
  cont_covars <- covariates[is_cont]
  cat_covars  <- covariates[!is_cont]

  COLORS <- c("#E74C3C", "#3498DB")
  out <- list(continuous = NULL, categorical = NULL)

  # --- Continuous: density plots ---
  if (length(cont_covars) > 0) {
    long_cont <- tidyr::pivot_longer(
      data[, c(treatment, cont_covars), drop = FALSE],
      cols      = -tidyr::all_of(treatment),
      names_to  = "covariate",
      values_to = "value"
    )
    p_cont <- ggplot2::ggplot(long_cont,
                              ggplot2::aes(x = value,
                                           color = .data[[treatment]],
                                           fill  = .data[[treatment]])) +
      ggplot2::geom_density(alpha = 0.15, linewidth = 1) +
      ggplot2::facet_wrap(~ covariate, scales = "free") +
      ggplot2::scale_color_manual(values = COLORS) +
      ggplot2::scale_fill_manual(values  = COLORS) +
      ggplot2::labs(title = "Continuous Covariates",
                    color = NULL, fill = NULL, x = NULL) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(legend.position = "top",
                     strip.text = ggplot2::element_text(face = "bold"))
    print(p_cont)
    out$continuous <- p_cont
  }

  # --- Categorical: dodged bar charts ---
  if (length(cat_covars) > 0) {
    long_cat <- tidyr::pivot_longer(
      data[, c(treatment, cat_covars), drop = FALSE],
      cols      = -tidyr::all_of(treatment),
      names_to  = "covariate",
      values_to = "value"
    )
    long_cat$value <- as.character(long_cat$value)
    p_cat <- ggplot2::ggplot(long_cat,
                             ggplot2::aes(x    = value,
                                          fill = .data[[treatment]])) +
      ggplot2::geom_bar(position = "dodge", alpha = 0.85) +
      ggplot2::facet_wrap(~ covariate, scales = "free_x") +
      ggplot2::scale_fill_manual(values = COLORS) +
      ggplot2::labs(title = "Categorical Covariates",
                    fill = NULL, x = NULL, y = "Count") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(legend.position = "top",
                     strip.text  = ggplot2::element_text(face = "bold"),
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    print(p_cat)
    out$categorical <- p_cat
  }

  invisible(out)
}

#' List Available Causal Contrasts
#'
#' @return A data.frame describing each supported contrast.
#' @export
available_contrasts <- function() {
  data.frame(
    contrast = c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE", "E[Y1]", "E[Y0]"),
    name     = c("Risk Difference", "Excess Relative Risk", "Risk Ratio",
                 "Odds Ratio", "Number Needed to Treat", "Survival Ratio",
                 "Vaccine Efficacy",
                 "Mean Potential Outcome under Treatment",
                 "Mean Potential Outcome under Control"),
    formula  = c("P(Y=1|A=1) - P(Y=1|A=0)",
                 "[P(Y=1|A=1) - P(Y=1|A=0)] / P(Y=1|A=0)",
                 "P(Y=1|A=1) / P(Y=1|A=0)",
                 "[P(Y=1|A=1)/P(Y=0|A=1)] / [P(Y=1|A=0)/P(Y=0|A=0)]",
                 "1 / |RD|",
                 "[1-P(Y=1|A=1)] / [1-P(Y=1|A=0)]",
                 "1 - RR",
                 "E[Y(1)]",
                 "E[Y(0)]"),
    stringsAsFactors = FALSE
  )
}

#' List Available SuperLearner Libraries
#'
#' @return A data.frame of algorithm names, categories, and
#'   required packages.
#' @export
available_sl_libraries <- function() {
  data.frame(
    algorithm        = c(
      # Parametric / GLM family
      "SL.glm", "SL.glm.interaction", "SL.bayesglm",
      "SL.step", "SL.step.forward", "SL.step.interaction", "SL.stepAIC",
      "SL.lm", "SL.speedglm", "SL.speedlm",
      "SL.ridge", "SL.glmnet",
      # GAM / Splines
      "SL.gam", "SL.earth", "SL.polymars", "SL.loess",
      # Tree-based
      "SL.rpart", "SL.rpartPrune", "SL.cforest",
      "SL.randomForest", "SL.ranger", "SL.ipredbagg",
      "SL.gbm", "SL.xgboost", "SL.bartMachine",
      # Kernel / KNN / SVM
      "SL.svm", "SL.ksvm", "SL.knn", "SL.kernelKnn",
      # Neural network
      "SL.nnet",
      # Discriminant analysis
      "SL.lda", "SL.qda",
      # Other
      "SL.biglasso", "SL.logreg", "SL.leekasso",
      "SL.caret", "SL.caret.rpart",
      "SL.mean", "SL.nnls"
    ),
    package_required = c(
      # Parametric / GLM family
      "stats", "stats", "arm",
      "stats", "stats", "stats", "MASS",
      "stats", "speedglm", "speedglm",
      "MASS", "glmnet",
      # GAM / Splines
      "mgcv", "earth", "polspline", "stats",
      # Tree-based
      "rpart", "rpart", "party",
      "randomForest", "ranger", "ipred",
      "gbm", "xgboost", "bartMachine",
      # Kernel / KNN / SVM
      "e1071", "kernlab", "class", "KernelKnn",
      # Neural network
      "nnet",
      # Discriminant analysis
      "MASS", "MASS",
      # Other
      "biglasso", "LogicReg", "leekasso",
      "caret", "caret",
      "SuperLearner", "nnls"
    )
  )
}

#' AIPW Bootstrap Confidence Intervals
#'
#' Nonparametric bootstrap wrapper around \code{\link{aipw_estimate}}.
#'
#' @inheritParams aipw_estimate
#' @param n_boot Number of bootstrap replicates. Default 200.
#'
#' @return A data.frame with columns: measure, estimate, lower, upper.
#' @export
aipw_bootstrap <- function(data,
                           outcome = "Y",
                           treatment = "A",
                           covariates,
                           sl_library_PS = NULL,
                           sl_library_OR = NULL,
                           n_folds = 5,
                           g_bounds = c(0.01, 0.99),
                           alpha = 0.05,
                           contrast = "all",
                           n_boot = 200) {

  # Original estimate
  orig <- aipw_estimate(
    data = data, outcome = outcome, treatment = treatment,
    covariates = covariates, sl_library_PS = sl_library_PS, sl_library_OR  = sl_library_OR,
    n_folds = n_folds, g_bounds = g_bounds, alpha = alpha,
    contrast = "all"
  )

  n <- nrow(data)
  measures <- c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE")
  boot_mat <- matrix(NA, nrow = n_boot, ncol = 7,
                     dimnames = list(NULL, measures))

  for (b in seq_len(n_boot)) {
    boot_data <- data[sample(n, n, replace = TRUE), ]
    boot_res <- tryCatch(
      aipw_estimate(
        data = boot_data, outcome = outcome, treatment = treatment,
        covariates = covariates, sl_library_PS = sl_library_PS, sl_library_OR  = sl_library_OR,
        n_folds = n_folds, g_bounds = g_bounds, alpha = alpha,
        contrast = "all"
      ),
      error = function(e) NULL
    )
    if (!is.null(boot_res)) {
      boot_mat[b, ] <- boot_res$estimate[match(measures, boot_res$measure)]
    }
  }

  probs <- c(alpha / 2, 1 - alpha / 2)

  # RD: quantile on original scale
  RD_ci <- stats::quantile(boot_mat[, "RD"], probs, na.rm = TRUE)

  # Log-scale CIs for RR, OR, NNT, SR
  log_measures <- c("RR", "OR", "NNT", "SR")
  log_boot <- log(boot_mat[, log_measures])
  log_ci <- apply(log_boot, 2, stats::quantile, probs = probs, na.rm = TRUE)

  # ERR = RR - 1
  ERR_ci <- exp(log_ci[, "RR"]) - 1

  # VE = 1 - RR (flip bounds)
  VE_ci <- 1 - exp(rev(log_ci[, "RR"]))

  results <- data.frame(
    measure  = measures,
    estimate = orig$estimate[match(measures, orig$measure)],
    lower    = c(RD_ci[1], ERR_ci[1], exp(log_ci[1, "RR"]),
                 exp(log_ci[1, "OR"]), exp(log_ci[1, "NNT"]),
                 exp(log_ci[1, "SR"]), VE_ci[1]),
    upper    = c(RD_ci[2], ERR_ci[2], exp(log_ci[2, "RR"]),
                 exp(log_ci[2, "OR"]), exp(log_ci[2, "NNT"]),
                 exp(log_ci[2, "SR"]), VE_ci[2]),
    stringsAsFactors = FALSE
  )
  rownames(results) <- NULL

  # Filter contrasts
  if (!identical(contrast, "all")) {
    results <- results[results$measure %in% contrast, ]
  }

  class(results) <- c("aipw_result", "data.frame")
  results
}

#' Plot AIPW Causal Estimate with Confidence Interval
#'
#' Generates a point estimate with vertical error bar for one or more
#' causal contrasts. Each measure is plotted as a separate figure.
#'
#' @param x An \code{aipw_result} object returned by \code{\link{aipw_estimate}}
#'   or \code{\link{aipw_bootstrap}}.
#' @param measures Character vector. Which measures to plot.
#'   Default \code{"all"} plots each measure as a separate figure.
#' @param null_line Logical. If \code{TRUE} (default), draw a dashed
#'   reference line at the null value.
#' @param save_dir Character or \code{NULL}. If provided, saves each plot
#'   to this directory as \code{<measure>_estimate.png}.
#' @param width Numeric. Width in inches. Default 4.
#' @param height Numeric. Height in inches. Default 5.
#'
#' @return A named list of ggplot2 objects (invisibly).
#' @export
plot_estimates <- function(x,
                           measures = "all",
                           null_line = TRUE,
                           save_dir = NULL,
                           width = 4,
                           height = 5) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  }

  df <- as.data.frame(x)

  # Filter measures
  if (!identical(measures, "all")) {
    valid <- df$measure
    bad <- setdiff(measures, valid)
    if (length(bad) > 0) {
      stop("Measure(s) not found in results: ", paste(bad, collapse = ", "))
    }
    df <- df[df$measure %in% measures, ]
  }

  # Null values
  additive <- c("RD", "VE", "ERR")
  get_null <- function(m) ifelse(m %in% additive, 0, 1)

  # Full names for titles
  full_names <- c(
    RD  = "Risk Difference",
    ERR = "Excess Relative Risk",
    RR  = "Risk Ratio",
    OR  = "Odds Ratio",
    NNT = "Number Needed to Treat",
    SR  = "Survival Ratio",
    VE  = "Vaccine Efficacy"
  )

  # Create save directory if needed
  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  plots <- list()

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    m <- row$measure
    null_val <- get_null(m)
    label <- sprintf("%.4f\n[%.4f, %.4f]", row$estimate, row$lower, row$upper)
    title <- ifelse(m %in% names(full_names), full_names[m], m)

    plot_df <- data.frame(
      x = 1,
      est = row$estimate,
      lo = row$lower,
      hi = row$upper
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = est)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lo, ymax = hi),
        width = 0.15, linewidth = 0.8
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        hjust = -0.3, size = 3.5
      ) +
      ggplot2::labs(
        title = title,
        y = m,
        x = NULL
      ) +
      ggplot2::xlim(0.5, 2) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )

    if (null_line) {
      p <- p + ggplot2::geom_hline(
        yintercept = null_val,
        linetype = "dashed", color = "red", alpha = 0.6
      )
    }

    print(p)
    plots[[m]] <- p

    # Save
    if (!is.null(save_dir)) {
      filepath <- file.path(save_dir, paste0(m, "_estimate.png"))
      ggplot2::ggsave(filepath, plot = p, width = width, height = height, dpi = 300)
      message("Saved: ", filepath)
    }
  }

  invisible(plots)
}
