get_exposure_levels <- function(
    data,
    exposure,
    exposure_type,
    exposure_quantiles = c(0.20, 0.80)
) {

  x <- data[[exposure]]

  if (exposure_type == "binary") {

    levels <- sort(unique(x[!is.na(x)]))

    if (length(levels) != 2) {
      stop("Binary exposure does not have exactly two observed levels.")
    }

    return(
      c(
        low = levels[1],
        high = levels[2]
      )
    )
  }

  q <- quantile(
    x,
    probs = exposure_quantiles,
    na.rm = TRUE,
    names = FALSE
  )

  if (!all(is.finite(q)) || q[1] >= q[2]) {
    stop("Continuous exposure has unusable intervention quantiles.")
  }

  c(
    low = q[1],
    high = q[2]
  )
}
