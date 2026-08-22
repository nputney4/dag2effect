estimate_gformula <- function(
    data,
    exposure,
    outcome,
    adjustment_set,
    exposure_type,
    outcome_type,
    exposure_quantiles = c(0.20, 0.80)
) {

  levels <- get_exposure_levels(
    data = data,
    exposure = exposure,
    exposure_type = exposure_type,
    exposure_quantiles = exposure_quantiles
  )

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
    estimate = mu_high - mu_low,
    mu_low = mu_low,
    mu_high = mu_high,
    a_low = levels["low"],
    a_high = levels["high"]
  )
}

bootstrap_gformula <- function(
    data,
    exposure,
    outcome,
    exposure_type,
    outcome_type,
    adjustment_sets,
    exposure_quantiles = c(0.20, 0.80),
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
          exposure_quantiles = exposure_quantiles
        )$estimate,
        error = function(e) NA_real_
      )
    }
  )

  # ---------------------------------------------------------
  # Bootstrap estimates
  #
  # One bootstrap sample is generated and then used for ALL
  # adjustment strategies.
  # ---------------------------------------------------------

  boot_estimates <- lapply(
    adjustment_sets,
    function(x) rep(NA_real_, B)
  )

  for (b in seq_len(B)) {

    # Same participant resample for every adjustment strategy
    idx <- sample.int(
      n = nrow(data),
      size = nrow(data),
      replace = TRUE
    )

    boot_data <- data[idx, , drop = FALSE]

    for (strategy in names(adjustment_sets)) {

      boot_estimates[[strategy]][b] <- tryCatch({

        estimate_gformula(
          data = boot_data,
          exposure = exposure,
          outcome = outcome,
          adjustment_set = adjustment_sets[[strategy]],
          exposure_type = exposure_type,
          outcome_type = outcome_type,
          exposure_quantiles = exposure_quantiles
        )$estimate

      }, error = function(e) {
        NA_real_
      })
    }
  }

  # Confidence intervals and summary for each strategy
  results <- lapply(names(adjustment_sets), function(strategy) {

    estimates <- boot_estimates[[strategy]]
    estimates <- estimates[is.finite(estimates)]

    ci <- if (length(estimates) > 0) {
      stats::quantile(
        estimates,
        probs = c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      )
    } else {
      c(NA_real_, NA_real_)
    }

    list(
      adjustment = strategy,
      estimate = point_estimates[[strategy]],
      ci_lower = ci[1],
      ci_upper = ci[2],
      ci_width = ci[2] - ci[1],
      n_boot_success = length(estimates)
    )
  })

  names(results) <- names(adjustment_sets)

  results
}
