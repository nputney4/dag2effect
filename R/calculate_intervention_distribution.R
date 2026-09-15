#' Simulate data under an intervention
#'
#' @param N Number of observations to simulate.
#' @param T Number of time points.
#' @param H Baseline-variable adjacency or coefficient matrix.
#' @param B Within-time exposure-variable matrix.
#' @param Gamma Baseline-to-exposure coefficient matrix.
#' @param A Lagged exposure-variable coefficient array.
#' @param intervention Named intervention values.
#' @param type_X Named or ordered exposure-variable types.
#' @param type_G Named or ordered baseline-variable types.
#' @param L Number of lags.
#' @param seed Optional random seed.
#' @param add_time_suffix Whether to append time suffixes to output columns.
#' @param ... Additional arguments passed to the data-generating functions.
#'
#' @returns A simulated wide data frame.
#' @export
calculate_intervention_distribution <- function(N,
                                                T,
                                                H,
                                                B,
                                                Gamma,
                                                A = NULL,
                                                intervention = NULL,
                                                type_X,
                                                type_G,
                                                L = NULL,
                                                seed = NULL,
                                                add_time_suffix = NULL,
                                                ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  strip_suffix <- function(x) sub("_(l|t)[0-9]+$", "", x)

  lookup_type <- function(node_name, type_vector, idx_fallback) {
    if (!is.null(names(type_vector))) {
      if (node_name %in% names(type_vector)) {
        return(as.character(type_vector[[node_name]]))
      }
      base_name <- strip_suffix(node_name)
      if (base_name %in% names(type_vector)) {
        return(as.character(type_vector[[base_name]]))
      }
    }

    if (length(type_vector) >= idx_fallback) {
      return(as.character(type_vector[[idx_fallback]]))
    }

    stop(paste0("Could not resolve type for node: ", node_name))
  }

  p <- if (is.null(B)) 0L else nrow(B)
  q <- if (is.null(H)) 0L else nrow(H)

  if (p > 0 && is.null(rownames(B))) stop("B must have row names.")
  if (q > 0 && is.null(rownames(H))) stop("H must have row names.")

  if (p > 0 && q > 0) {
    if (is.null(Gamma) || nrow(Gamma) != p || ncol(Gamma) != q) {
      stop("Gamma must be nrow(B) x nrow(H).")
    }
  }

  if (is.null(A)) {
    L <- 0L
  } else {
    if (length(dim(A)) != 3) stop("A must be a 3D array.")
    if (dim(A)[1] != p || dim(A)[2] != p) {
      stop("A must be nrow(B) x nrow(B) x L.")
    }
    if (is.null(L)) L <- dim(A)[3]
    if (L > dim(A)[3]) stop("Provided L is larger than dim(A)[3].")
  }

  x_are_time_explicit <- p > 0 && all(grepl("_t[[:digit:]]+$", rownames(B)))
  if (x_are_time_explicit && T != 1) {
    stop("When B contains X*_t* nodes, use T = 1 and A = NULL.")
  }

  if (is.null(add_time_suffix)) {
    add_time_suffix <- !x_are_time_explicit
  }

  intervention_lookup <- NULL
  if (!is.null(intervention)) {
    if (is.null(names(intervention)) || any(names(intervention) == "")) {
      stop("intervention must be a named list/vector.")
    }
    intervention_lookup <- intervention
  }

  dist_args <- list(...)
  dist_params <- dist_args$dist_params
  dist_args$dist_params <- NULL
  G <- matrix(0, nrow = N, ncol = q)
  X <- array(0, dim = c(N, T, p))

  if (q > 0) {
    for (k in seq_len(q)) {
      g_name <- rownames(H)[k]
      g_type <- lookup_type(g_name, type_G, k)

      if (!is.null(intervention_lookup) && g_name %in% names(intervention_lookup)) {
        G[, k] <- intervention_lookup[[g_name]]
        next
      }

      parent_idx <- which(H[k, ] != 0)
      if (length(parent_idx) > 0 && any(parent_idx >= k)) {
        stop("H is not topologically ordered.")
      }

      if (length(parent_idx) == 0) {
        G[, k] <- do.call(
          randMixDAG::eta_to_obs,
          c(list(root = TRUE, type = g_type, N = N,
                 node_name = g_name, dist_params = dist_params), dist_args)
        )
      } else {
        parent_values <- matrix(G[, parent_idx, drop = FALSE], nrow = N, ncol = length(parent_idx))
        eta_k <- parent_values %*% matrix(H[k, parent_idx], ncol = 1)
        G[, k] <- do.call(
          randMixDAG::eta_to_obs,
          c(list(root = FALSE, type = g_type, eta = eta_k, N = N,
                 node_name = g_name, dist_params = dist_params), dist_args)
        )
      }
    }
  }

  if (p > 0) {
    for (t in seq_len(T)) {
      for (j in seq_len(p)) {
        x_base <- rownames(B)[j]
        x_node <- if (x_are_time_explicit) x_base else paste0(x_base, "_t", t)

        if (!is.null(intervention_lookup)) {
          if (x_node %in% names(intervention_lookup)) {
            X[, t, j] <- intervention_lookup[[x_node]]
            next
          }
          if (!x_are_time_explicit && x_base %in% names(intervention_lookup)) {
            X[, t, j] <- intervention_lookup[[x_base]]
            next
          }
        }

        eta_tj <- matrix(0, nrow = N, ncol = 1)

        parent_idx <- which(B[j, ] != 0)
        if (length(parent_idx) > 0) {
          if (any(parent_idx >= j)) {
            stop("B is not topologically ordered.")
          }
          parent_values <- matrix(X[, t, parent_idx, drop = FALSE], nrow = N, ncol = length(parent_idx))
          eta_tj <- eta_tj + parent_values %*% matrix(B[j, parent_idx], ncol = 1)
        }

        if (q > 0) {
          eta_tj <- eta_tj + G %*% matrix(Gamma[j, ], ncol = 1)
        }

        if (!is.null(A) && t > 1 && L > 0) {
          for (l in seq_len(min(t - 1, L))) {
            lag_values <- matrix(X[, t - l, , drop = FALSE], nrow = N, ncol = p)
            eta_tj <- eta_tj + lag_values %*% matrix(A[j, , l], ncol = 1)
          }
        }

        x_type <- lookup_type(x_base, type_X, j)
        X[, t, j] <- do.call(
          randMixDAG::eta_to_obs,
          c(list(root = FALSE, type = x_type, eta = eta_tj, N = N,
                 node_name = x_node, dist_params = dist_params), dist_args)
        )
      }
    }
  }

  randMixDAG::data_to_wide_df(
    X,
    G,
    H,
    B,
    Gamma,
    A,
    add_time_suffix = add_time_suffix
  )
}
