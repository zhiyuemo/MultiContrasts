#' @import SuperLearner
NULL


#' @keywords internal
estimate_ps_sl <- function(train_data, test_data, covariates, sl_library) {
  X_train <- train_data[, covariates, drop = FALSE]
  X_test <- test_data[, covariates, drop = FALSE]
  sl_fit <- SuperLearner::SuperLearner(
    Y = train_data$A, X = X_train,
    family = stats::binomial(),
    SL.library = sl_library,
    verbose = FALSE
  )
  as.numeric(SuperLearner::predict.SuperLearner(sl_fit, newdata = X_test)$pred)
}

#' @keywords internal
estimate_outcome_sl <- function(train_data, test_data, covariates, sl_library,
                                outcome_type = "binary") {
  train_1 <- train_data[train_data$A == 1, ]
  train_0 <- train_data[train_data$A == 0, ]

  X_1 <- train_1[, covariates, drop = FALSE]
  X_0 <- train_0[, covariates, drop = FALSE]
  X_test <- test_data[, covariates, drop = FALSE]

  fam <- if (outcome_type == "binary") stats::binomial() else stats::gaussian()

  sl_fit_1 <- SuperLearner::SuperLearner(
    Y = train_1$Y, X = X_1,
    family = fam,
    SL.library = sl_library, method = "method.NNLS", verbose = FALSE
  )
  sl_fit_0 <- SuperLearner::SuperLearner(
    Y = train_0$Y, X = X_0,
    family = fam,
    SL.library = sl_library, method = "method.NNLS", verbose = FALSE
  )

  mu_1 <- as.numeric(SuperLearner::predict.SuperLearner(sl_fit_1, newdata = X_test)$pred)
  mu_0 <- as.numeric(SuperLearner::predict.SuperLearner(sl_fit_0, newdata = X_test)$pred)
  list(mu_1 = mu_1, mu_0 = mu_0)
}

#' @keywords internal
compute_contrasts <- function(aipw_1_n, aipw_0_n, alpha = 0.05) {
  n <- length(aipw_1_n)
  psi_1 <- mean(aipw_1_n)
  psi_0 <- mean(aipw_0_n)

  RD  <- psi_1 - psi_0
  RR  <- psi_1 / psi_0
  ERR <- (psi_1 - psi_0) / psi_0
  OR  <- (psi_1 / (1 - psi_1)) / (psi_0 / (1 - psi_0))
  NNT <- abs(1 / RD)
  SR  <- (1 - psi_1) / (1 - psi_0)
  VE  <- 1 - RR

  IF_psi1 <- aipw_1_n - psi_1
  IF_psi0 <- aipw_0_n - psi_0

  IF_RD     <- IF_psi1 - IF_psi0
  IF_logRR  <- (1 / psi_1) * IF_psi1 - (1 / psi_0) * IF_psi0
  IF_logOR  <- (1 / psi_1 + 1 / (1 - psi_1)) * IF_psi1 +
    (-1 / (1 - psi_0) - 1 / psi_0) * IF_psi0
  IF_logNNT <- (1 / (psi_0 - psi_1)) * IF_psi1 +
    (-1 / (psi_0 - psi_1)) * IF_psi0
  IF_logSR  <- (-1 / (1 - psi_1)) * IF_psi1 +
    (1 / (1 - psi_0)) * IF_psi0

  se_RD  <- sqrt(stats::var(IF_RD) / n)
  se_psi1 <- sqrt(stats::var(IF_psi1) / n)
  se_psi0 <- sqrt(stats::var(IF_psi0) / n)
  se_log <- sapply(
    list(IF_logRR, IF_logOR, IF_logNNT, IF_logSR),
    function(x) sqrt(stats::var(x) / n)
  )

  z <- stats::qnorm(1 - alpha / 2)

  # RD: Wald CI on identity scale
  RD_lower <- RD - z * se_RD
  RD_upper <- RD + z * se_RD

  # Psi1
  psi1_lower <- psi_1 - z * se_psi1
  psi1_upper <- psi_1 + z * se_psi1

  # Psi0
  psi0_lower <- psi_0 - z * se_psi0
  psi0_upper <- psi_0 + z * se_psi0

  # Log-scale CIs for RR, OR, NNT, SR
  log_est <- log(c(RR = RR, OR = OR, NNT = NNT, SR = SR))
  log_lower <- log_est - z * se_log
  log_upper <- log_est + z * se_log

  # NNT = 1/RD
  if(RD_lower*RD_upper < 0){
    NNT_lower <- exp(log_lower["NNT"])
    NNT_upper <- exp(log_upper["NNT"])
  }else if(RD_lower < 0){
    NNT_lower <- abs(1/RD_lower)
    NNT_upper <- abs(1/RD_upper)
  }else{
    NNT_lower <- abs(1/RD_upper)
    NNT_upper <- abs(1/RD_lower)
  }

  # ERR = RR - 1
  ERR_lower <- exp(log_lower[1]) - 1
  ERR_upper <- exp(log_upper[1]) - 1

  # VE = 1 - RR (bounds flip)
  VE_lower <- 1 - exp(log_upper[1])
  VE_upper <- 1 - exp(log_lower[1])

  data.frame(
    measure   = c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE","E[Y1]","E[Y0]"),
    estimate  = c(RD,   ERR,   RR,   OR,   NNT,   SR,   VE , psi_1, psi_0),
    lower     = c(RD_lower, ERR_lower, exp(log_lower["RR"]),
                  exp(log_lower["OR"]), NNT_lower,
                  exp(log_lower["SR"]), VE_lower, psi1_lower, psi0_lower),
    upper     = c(RD_upper, ERR_upper, exp(log_upper["RR"]),
                  exp(log_upper["OR"]), NNT_upper,
                  exp(log_upper["SR"]), VE_upper,psi1_upper, psi0_upper),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
compute_contrasts_continuous <- function(aipw_1_n, aipw_0_n, alpha = 0.05) {
  n <- length(aipw_1_n)
  psi_1 <- mean(aipw_1_n)
  psi_0 <- mean(aipw_0_n)

  mean_diff  <- psi_1 - psi_0
  mean_ratio <- psi_1 / psi_0

  IF_psi1 <- aipw_1_n - psi_1
  IF_psi0 <- aipw_0_n - psi_0
  IF_diff <- IF_psi1 - IF_psi0
  IF_log_ratio <- (1 / psi_1) * IF_psi1 - (1 / psi_0) * IF_psi0

  se_diff  <- sqrt(stats::var(IF_diff) / n)
  se_psi1  <- sqrt(stats::var(IF_psi1) / n)
  se_psi0  <- sqrt(stats::var(IF_psi0) / n)
  se_log_ratio <- sqrt(stats::var(IF_log_ratio) / n)

  z <- stats::qnorm(1 - alpha / 2)

  log_ratio <- log(mean_ratio)
  ratio_lower <- exp(log_ratio - z * se_log_ratio)
  ratio_upper <- exp(log_ratio + z * se_log_ratio)

  data.frame(
    measure  = c("mean_diff", "mean_ratio", "E[Y1]", "E[Y0]"),
    estimate = c(mean_diff, mean_ratio, psi_1, psi_0),
    lower    = c(mean_diff - z * se_diff, ratio_lower,
                 psi_1 - z * se_psi1, psi_0 - z * se_psi0),
    upper    = c(mean_diff + z * se_diff, ratio_upper,
                 psi_1 + z * se_psi1, psi_0 + z * se_psi0),
    stringsAsFactors = FALSE
  )
}

#' AIPW Estimator with Cross-Fitting
#'
#' @param data A data.frame with columns named by \code{outcome},
#'   \code{treatment}, and \code{covariates}.
#' @param outcome Character. Name of binary outcome column (0/1).
#' @param treatment Character. Name of binary treatment column (0/1).
#' @param covariates Character vector of covariate column names.
#' @param sl_library_PS Character vector of SuperLearner algorithms
#'   for the propensity score model.
#' @param sl_library_OR Character vector of SuperLearner algorithms
#'   for the outcome regression model.
#' @param n_folds Number of cross-fitting folds. Default 5.
#' @param g_bounds Numeric vector of length 2. Truncation bounds
#'   for propensity scores. Default \code{c(0.01, 0.99)}.
#' @param alpha Significance level for CIs. Default 0.05.
#' @param contrast Which contrasts to report. Default \code{"all"}.
#'   Can also be a character vector subset of
#'   \code{c("RD","ERR","RR","OR","NNT","SR","VE")}.
#' @param population A logical vector of length \code{nrow(data)}, or
#'   \code{NULL} (default). When supplied, nuisance models
#'   (\eqn{\hat{\mu}_1}, \eqn{\hat{\mu}_0}, \eqn{\hat{g}}) are fitted via
#'   cross-fitting on the \strong{full} dataset, yielding individual-level
#'   AIPW pseudo-outcomes for every observation. Causal contrasts are then
#'   computed by averaging those pseudo-outcomes only within the subgroup
#'   \code{population == TRUE}, i.e.\cr
#'   \deqn{\hat{\psi}^S = \frac{1}{n_S}\sum_{i: S_i = 1} \phi(O_i)}
#'   This gives the conditional ATE (CATE) / average treatment effect in
#'   the target subpopulation, using the efficiency of the full-data nuisance
#'   estimates. Influence-function variance is computed within the subgroup.
#'
#' @return A data.frame with columns: measure, estimate, lower, upper.
#' @export
aipw_estimate <- function(data,
                          outcome = "Y",
                          treatment = "A",
                          covariates,
                          sl_library_PS = NULL,
                          sl_library_OR = NULL,
                          n_folds = 5,
                          g_bounds = c(0.01, 0.99),
                          alpha = 0.05,
                          contrast = "all",
                          population = NULL,
                          outcome_type = c("binary", "continuous")) {

  outcome_type <- match.arg(outcome_type)

  data <- as.data.frame(data)

  # Validate population argument
  if (!is.null(population)) {
    if (!is.logical(population) || length(population) != nrow(data)) {
      stop("`population` must be a logical vector with length equal to nrow(data).")
    }
    if (sum(population) < 5) {
      stop("Fewer than 5 observations in the target `population`. Expand the subgroup.")
    }
  }

  # --- Rename columns internally to A and Y ---
  data$A <- data[[treatment]]
  data$Y <- as.numeric(data[[outcome]])

  n <- nrow(data)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  aipw_1_n <- numeric(n)
  aipw_0_n <- numeric(n)

  for (k in 1:n_folds) {
    train_data <- data[fold_ids != k, ]
    test_data  <- data[fold_ids == k, ]
    test_idx   <- which(fold_ids == k)

    # Propensity score
    ps_hat <- estimate_ps_sl(train_data, test_data, covariates, sl_library_PS)
    ps_hat <- pmax(pmin(ps_hat, g_bounds[2]), g_bounds[1])

    # Outcome model
    mu_hats <- estimate_outcome_sl(train_data, test_data, covariates, sl_library_OR,
                                    outcome_type = outcome_type)

    # AIPW pseudo-outcomes for every observation in this fold
    aipw_1_n[test_idx] <- mu_hats$mu_1 +
      (test_data$A / ps_hat) * (test_data$Y - mu_hats$mu_1)
    aipw_0_n[test_idx] <- mu_hats$mu_0 +
      ((1 - test_data$A) / (1 - ps_hat)) * (test_data$Y - mu_hats$mu_0)
  }

  # If a target population is specified, restrict pseudo-outcomes to that
  # subgroup before computing contrasts.  Nuisance models were fitted on the
  # full data, so their predictions already benefit from the full sample.
  if (!is.null(population)) {
    aipw_1_n <- aipw_1_n[population]
    aipw_0_n <- aipw_0_n[population]
  }

  # Compute contrasts with influence-function-based CIs
  if (outcome_type == "binary") {
    results <- compute_contrasts(aipw_1_n, aipw_0_n, alpha = alpha)
    valid <- c("RD", "ERR", "RR", "OR", "NNT", "SR", "VE", "E[Y1]", "E[Y0]")
  } else {
    results <- compute_contrasts_continuous(aipw_1_n, aipw_0_n, alpha = alpha)
    valid <- c("mean_diff", "mean_ratio", "E[Y1]", "E[Y0]")
  }

  # Filter contrasts if requested
  if (!identical(contrast, "all")) {
    bad <- setdiff(contrast, valid)
    if (length(bad) > 0) {
      stop("Unknown contrast(s): ", paste(bad, collapse = ", "),
           ". Choose from: ", paste(valid, collapse = ", "))
    }
    results <- results[results$measure %in% contrast, ]
  }

  rownames(results) <- NULL
  class(results) <- c("aipw_result", "data.frame")
  results
}

#' Print method for aipw_result
#' @param x An aipw_result object.
#' @param digits Number of digits. Default 4.
#' @param ... Ignored.
#' @export
print.aipw_result <- function(x, digits = 4, ...) {
  cat("\n=== AIPW Causal Estimates ===\n\n")
  fmt <- x
  fmt$estimate <- round(fmt$estimate, digits)
  fmt$lower    <- round(fmt$lower, digits)
  fmt$upper    <- round(fmt$upper, digits)
  print.data.frame(fmt, row.names = FALSE)
  cat("\n")
  invisible(x)
}
