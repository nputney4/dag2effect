#' Estimate intervention means using the parametric g-formula
#'
#' @param data Data frame containing the exposure and outcome.
#' @param exposure Name of the exposure column.
#' @param outcome Name of the outcome column.
#' @param adjustment_set Covariates used in the outcome model.
#' @param exposure_type Exposure type.
#' @param outcome_type Outcome type.
#' @param exposure_quantiles Quantiles used for non-binary exposures.
#' @param exposure_levels Optional named vector with fixed `low` and `high`
#'   intervention values.
#'
#' @returns A list containing the estimated low- and high-intervention means.
#' @export
estimate_gformula <- function(
    data,
    exposure,
    outcome,
    adjustment_set,
    exposure_type,
    outcome_type,
    exposure_quantiles = c(0.20, 0.80),
    exposure_levels = NULL
) {

  levels <- if (is.null(exposure_levels)) {
    get_exposure_levels(
      data = data,
      exposure = exposure,
      exposure_type = exposure_type,
      exposure_quantiles = exposure_quantiles
    )
  } else {
    if (!all(c("low", "high") %in% names(exposure_levels))) {
      stop("exposure_levels must contain named low and high values.")
    }
    exposure_levels[c("low", "high")]
  }

  # Outcome model:
  # outcome is a specific full-time node, and adjustment_set
  # contains specific full-time nodes.
  adjustment_set <- setdiff(
    adjustment_set,
    c(exposure, outcome)
  )

  rhs <- unique(
    c(exposure, adjustment_set)
  )

  form <- stats::reformulate(
    termlabels = rhs,
    response = outcome
  )

  fit <- if (outcome_type == "continuous") {

    stats::lm(
      formula = form,
      data = data
    )

  } else {

    stats::glm(
      formula = form,
      data = data,
      family = stats::binomial()
    )
  }

  # Two copies of the SAME longitudinal dataset.
  data_low <- data
  data_high <- data

  # Intervene only on the selected full-time node.
  data_low[[exposure]] <- levels["low"]
  data_high[[exposure]] <- levels["high"]

  if (outcome_type == "continuous") {

    mu_low <- mean(
      stats::predict(fit, newdata = data_low),
      na.rm = TRUE
    )

    mu_high <- mean(
      stats::predict(fit, newdata = data_high),
      na.rm = TRUE
    )

  } else {

    mu_low <- mean(
      stats::predict(
        fit,
        newdata = data_low,
        type = "response"
      ),
      na.rm = TRUE
    )

    mu_high <- mean(
      stats::predict(
        fit,
        newdata = data_high,
        type = "response"
      ),
      na.rm = TRUE
    )
  }

  list(
    mu_low = mu_low,
    mu_high = mu_high
  )
}

#' Bootstrap intervention means from the parametric g-formula
#'
#' @param data Data frame containing the exposure and outcome.
#' @param exposure Name of the exposure column.
#' @param outcome Name of the outcome column.
#' @param exposure_type Exposure type.
#' @param outcome_type Outcome type.
#' @param adjustment_sets Named list of adjustment sets.
#' @param exposure_quantiles Quantiles used for non-binary exposures.
#' @param exposure_levels Optional named vector with fixed `low` and `high`
#'   intervention values shared by all bootstrap samples.
#' @param B Number of bootstrap samples.
#'
#' @returns A named list containing point estimates, bootstrap intervention
#'   means, and confidence intervals for each intervention mean.
#' @export
bootstrap_gformula <- function(
    data,
    exposure,
    outcome,
    exposure_type,
    outcome_type,
    adjustment_sets,
    exposure_quantiles = c(0.20, 0.80),
    exposure_levels = NULL,
    B = 100
) {
  # ---------------------------------------------------------
  # Point estimates
  # ---------------------------------------------------------

  point_estimates <- lapply(
    adjustment_sets,
    function(adjustment_set) {

      tryCatch(
        estimate_gformula(
          data = data,
          exposure = exposure,
          outcome = outcome,
          adjustment_set = adjustment_set,
          exposure_type = exposure_type,
          outcome_type = outcome_type,
          exposure_quantiles = exposure_quantiles,
          exposure_levels = exposure_levels
        ),
        error = function(e) {
          list(mu_low = NA_real_, mu_high = NA_real_)
        }
      )
    }
  )

  # ---------------------------------------------------------
  # Bootstrap estimates
  #
  # One bootstrap sample is generated and then used for ALL
  # adjustment strategies.
  # ---------------------------------------------------------

  boot_mu_low <- lapply(adjustment_sets, function(x) rep(NA_real_, B))
  boot_mu_high <- lapply(adjustment_sets, function(x) rep(NA_real_, B))

  for (b in seq_len(B)) {

    # Same participant resample for every adjustment strategy
    idx <- sample.int(
      n = nrow(data),
      size = nrow(data),
      replace = TRUE
    )

    boot_data <- data[idx, , drop = FALSE]

    for (strategy in names(adjustment_sets)) {

      boot_estimate <- tryCatch({

        estimate_gformula(
          data = boot_data,
          exposure = exposure,
          outcome = outcome,
          adjustment_set = adjustment_sets[[strategy]],
          exposure_type = exposure_type,
          outcome_type = outcome_type,
          exposure_quantiles = exposure_quantiles,
          exposure_levels = exposure_levels
        )

      }, error = function(e) {
        list(mu_low = NA_real_, mu_high = NA_real_)
      })

      boot_mu_low[[strategy]][b] <- boot_estimate$mu_low
      boot_mu_high[[strategy]][b] <- boot_estimate$mu_high
    }
  }

  # Confidence intervals and summary for each intervention mean
  results <- lapply(names(adjustment_sets), function(strategy) {

    valid_estimates <- is.finite(boot_mu_low[[strategy]]) &
      is.finite(boot_mu_high[[strategy]])
    low_estimates <- boot_mu_low[[strategy]][valid_estimates]
    high_estimates <- boot_mu_high[[strategy]][valid_estimates]

    low_ci <- if (length(low_estimates) > 0) {
      stats::quantile(
        low_estimates,
        probs = c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      )
    } else {
      c(NA_real_, NA_real_)
    }

    high_ci <- if (length(high_estimates) > 0) {
      stats::quantile(
        high_estimates,
        probs = c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      )
    } else {
      c(NA_real_, NA_real_)
    }

    list(
      adjustment = strategy,
      mu_low = point_estimates[[strategy]]$mu_low,
      mu_high = point_estimates[[strategy]]$mu_high,
      bootstrap_mu_low = low_estimates,
      bootstrap_mu_high = high_estimates,
      mu_low_ci_lower = low_ci[1],
      mu_low_ci_upper = low_ci[2],
      mu_high_ci_lower = high_ci[1],
      mu_high_ci_upper = high_ci[2],
      n_boot_success = min(length(low_estimates), length(high_estimates))
    )
  })

  names(results) <- names(adjustment_sets)

  results
}
