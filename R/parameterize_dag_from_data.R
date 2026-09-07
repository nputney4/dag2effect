get_dist_from_var_type <- function(var_type) {
  if (var_type == "continuous") return("gaussian")
  if (var_type == "binary") return("binomial")
  if (var_type == "poisson") return("poisson")
  stop(paste0("Unknown variable type: ", var_type))
}

subset_adj_block <- function(adj_coeffs, row_nodes, col_nodes) {
  out <- matrix(
    0,
    nrow = length(row_nodes),
    ncol = length(col_nodes),
    dimnames = list(row_nodes, col_nodes)
  )

  row_present <- intersect(row_nodes, rownames(adj_coeffs))
  col_present <- intersect(col_nodes, colnames(adj_coeffs))

  if (length(row_present) > 0 && length(col_present) > 0) {
    out[row_present, col_present] <- adj_coeffs[row_present, col_present, drop = FALSE]
  }

  out
}

fit_conditional_model <- function(df, response, rhs, var_type) {
  fma <- stats::as.formula(paste(response, "~", rhs))

  if (var_type == "ordinal") {
    model_df <- df
    response_values <- model_df[[response]]
    if (is.factor(response_values)) {
      observed_levels <- levels(droplevels(response_values))
    } else {
      observed_levels <- sort(unique(response_values[!is.na(response_values)]))
    }
    if (length(observed_levels) < 2) {
      stop(paste0("Ordinal node ", response, " must have at least 2 observed levels."))
    }
    model_df[[response]] <- factor(response_values, levels = observed_levels, ordered = TRUE)
    return(MASS::polr(fma, data = model_df, method = "probit", Hess = TRUE))
  }

  stats::glm(fma, family = get_dist_from_var_type(var_type), data = df)
}

extract_model_coeffs <- function(mdl, var_type) {
  coeffs <- stats::coef(mdl)
  if (var_type != "ordinal") {
    coeffs <- coeffs[names(coeffs) != "(Intercept)"]
  }
  coeffs <- coeffs[!is.na(coeffs)]
  coeffs
}

extract_model_dist_params <- function(mdl, var_type) {
  out <- list(
    node_intercept = 0,
    cont_sd = NULL,
    ord_breaks = NULL
  )

  if (var_type == "ordinal") {
    if (is.null(mdl$zeta) || length(mdl$zeta) == 0) {
      stop("Failed to extract ordinal thresholds (zeta) from fitted model.")
    }
    out$ord_breaks <- c(-Inf, as.numeric(mdl$zeta), Inf)
    return(out)
  }

  coeffs <- stats::coef(mdl)
  intercept <- unname(coeffs["(Intercept)"])
  if (!is.na(intercept)) {
    out$node_intercept <- intercept
  }

  if (var_type == "continuous") {
    dispersion <- summary(mdl)$dispersion
    if (!is.na(dispersion) && is.finite(dispersion) && dispersion > 0) {
      out$cont_sd <- sqrt(dispersion)
    }
  }

  out
}

fit_adj_coeffs_from_dag <- function(df, dag, var_types) {
  vars_in_order <- caugi::topological_sort(dag)
  adj <- t(caugi::as_adjacency(dag)) # rows = child, cols = parent

  adj_coeffs <- matrix(
    0,
    nrow = nrow(adj),
    ncol = ncol(adj),
    dimnames = dimnames(adj)
  )

  mdl_list <- vector("list", length(vars_in_order))
  names(mdl_list) <- vars_in_order

  dist_params <- list(
    node_intercepts = setNames(numeric(0), character(0)),
    cont_sd = setNames(numeric(0), character(0)),
    ord_breaks = list()
  )

  for (vbl in vars_in_order) {
    vbl_base <- stringr::str_remove(vbl, "_(l|t)[[:digit:]]+$")
    var_type <- var_types[[vbl_base]]
    if (is.null(var_type)) {
      stop(paste0("No variable type found for node: ", vbl))
    }

    pa_var <- rownames(adj)[which(adj[vbl, ] > 0)]
    rhs <- if (length(pa_var) == 0) "1" else paste(pa_var, collapse = " + ")

    mdl <- fit_conditional_model(df = df, response = vbl, rhs = rhs, var_type = var_type)
    mdl_list[[vbl]] <- mdl

    node_dist <- extract_model_dist_params(mdl, var_type)
    dist_params$node_intercepts[vbl] <- node_dist$node_intercept
    if (!is.null(node_dist$cont_sd)) {
      dist_params$cont_sd[vbl] <- node_dist$cont_sd
    }
    if (!is.null(node_dist$ord_breaks)) {
      dist_params$ord_breaks[[vbl]] <- node_dist$ord_breaks
    }

    coeffs <- extract_model_coeffs(mdl, var_type)

    if (length(coeffs) > 0) {
      matched_coeffs <- intersect(names(coeffs), colnames(adj_coeffs))
      if (length(matched_coeffs) > 0) {
        adj_coeffs[vbl, matched_coeffs] <- coeffs[matched_coeffs]
      }
    }
  }

  list(
    adj = adj,
    adj_coeffs = adj_coeffs,
    models = mdl_list,
    dist_params = dist_params,
    vars_in_order = vars_in_order
  )
}

extract_matrices_lagged <- function(adj_coeffs, vars_in_order, g_nodes, x_nodes) {
  g_order <- intersect(vars_in_order, g_nodes)
  if (length(g_order) == 0) g_order <- g_nodes

  x_order <- intersect(vars_in_order, x_nodes)
  if (length(x_order) == 0) x_order <- x_nodes

  H <- subset_adj_block(adj_coeffs, g_order, g_order)
  B <- subset_adj_block(adj_coeffs, x_order, x_order)
  Gamma <- subset_adj_block(adj_coeffs, x_order, g_order)

  lag_candidates <- unique(c(
    grep("_l[[:digit:]]+$", rownames(adj_coeffs), value = TRUE),
    grep("_l[[:digit:]]+$", colnames(adj_coeffs), value = TRUE)
  ))
  lag_candidates <- lag_candidates[
    stringr::str_remove(lag_candidates, "_l[[:digit:]]+$") %in% x_order
  ]

  if (length(lag_candidates) == 0) {
    A <- NULL
    L <- 0L
  } else {
    lag_ids <- as.integer(stringr::str_match(lag_candidates, "_l([[:digit:]]+)$")[, 2])
    L <- max(lag_ids)

    A <- array(
      0,
      dim = c(length(x_order), length(x_order), L),
      dimnames = list(x_order, x_order, paste0("lag_", seq_len(L)))
    )

    row_present <- intersect(x_order, rownames(adj_coeffs))

    for (l in seq_len(L)) {
      lag_parent_nodes <- paste0(x_order, "_l", l)
      col_present <- intersect(lag_parent_nodes, colnames(adj_coeffs))
      if (length(col_present) == 0) next

      parent_base <- stringr::str_remove(col_present, "_l[[:digit:]]+$")
      A[row_present, parent_base, l] <- adj_coeffs[row_present, col_present, drop = FALSE]
    }
  }

  list(H = H, Gamma = Gamma, B = B, A = A, L = L, g_order = g_order, x_order = x_order)
}

extract_matrices_full_time <- function(adj_coeffs, vars_in_order, g_nodes, x_nodes) {
  g_order <- intersect(vars_in_order, g_nodes)
  if (length(g_order) == 0) g_order <- g_nodes

  x_time_order <- vars_in_order[grepl("_t[[:digit:]]+$", vars_in_order)]
  x_time_order <- x_time_order[
    stringr::str_remove(x_time_order, "_t[[:digit:]]+$") %in% x_nodes
  ]

  H <- subset_adj_block(adj_coeffs, g_order, g_order)
  B <- subset_adj_block(adj_coeffs, x_time_order, x_time_order)
  Gamma <- subset_adj_block(adj_coeffs, x_time_order, g_order)

  list(H = H, Gamma = Gamma, B = B, A = NULL, L = 0L, g_order = g_order, x_order = x_time_order)
}

data2dag <- function(df,
                                       dag,
                                       var_G,
                                       var_X,
                                       graph_type = c("lagged", "full_time")) {
  graph_type <- match.arg(graph_type)
  var_types <- c(var_G, var_X)

  fit <- fit_adj_coeffs_from_dag(df = df, dag = dag, var_types = var_types)

  matrices <- switch(
    graph_type,
    lagged = extract_matrices_lagged(
      adj_coeffs = fit$adj_coeffs,
      vars_in_order = fit$vars_in_order,
      g_nodes = names(var_G),
      x_nodes = names(var_X)
    ),
    full_time = extract_matrices_full_time(
      adj_coeffs = fit$adj_coeffs,
      vars_in_order = fit$vars_in_order,
      g_nodes = names(var_G),
      x_nodes = names(var_X)
    )
  )
  matrices$G <- matrices$H
  matrices$dist_params <- fit$dist_params

  list(
    adj = fit$adj,
    adj_coeffs = fit$adj_coeffs,
    models = fit$models,
    dist_params = fit$dist_params,
    matrices = matrices
  )
}

# Example: lagged and full-time ATEs from parameterized matrices.
if (FALSE) {
  path_to_repo <- "C:\\Users\\putne001\\OneDrive - Wageningen University & Research\\Desktop\\WageningenPhD_NPutney\\phd-core\\thesis\\"
  path_to_code <- paste0(path_to_repo, "chapter1\\code\\hrs-causal-discovery\\R\\")
  path_to_data <- paste0(path_to_repo, "chapter1\\data\\simulated-datasets\\mixed-data-random-dgp\\")
  path_to_results <- paste0(path_to_repo, "chapter1\\results\\sim-study\\estimated-graphs")

  # Required by calculate_intervention_distribution()
  source(file.path(path_to_code, "simulate_mixed_data_from_adjacency_matrices.R", fsep = "\\"))
  source(file.path(path_to_code, "calculate_intervention_distribution.R", fsep = "\\"))

  dataset_object <- readRDS(file.path(path_to_data, "datasets_N.rds", fsep = "\\"))
  df_lagged <- dataset_object$`500`$replicates[[1]]$df$data_l
  df_full <- dataset_object$`500`$replicates[[1]]$df$data_w

  var_G <- dataset_object$`500`$replicates[[1]]$var_types$type_G
  names(var_G) <- paste0("G", seq_along(var_G))
  var_X <- dataset_object$`500`$replicates[[1]]$var_types$type_X
  names(var_X) <- paste0("X", seq_along(var_X))

  dag_lagged <- caugi::dag_from_pdag(
    caugi::read_caugi(file.path(
      path_to_results,
      "datasetsN\\datasetsN_5000N_10d_5T_2L_lagged_rep2__boss_cg.caugi",
      fsep = "\\"
    ))
  )
  dag_full <- caugi::dag_from_pdag(
    caugi::read_caugi(file.path(
      path_to_results,
      "datasetsN\\datasetsN_5000N_10d_5T_2L_full_rep2__boss_cg.caugi",
      fsep = "\\"
    ))
  )

  dag_gt <- caugi::read_caugi(
    "C:/Users/putne001/OneDrive - Wageningen University & Research/Desktop/WageningenPhD_NPutney/phd-core/thesis/chapter1/results/sim-study/ground-truth-graphs/datasetsN/datasetsN_5000N_10d_5T_2L_full_rep2.caugi"
  )

  dag_lagged_unrolled <- unroll_lagged_graph(dag_lagged, T = 5)
  shd(dag_gt, dag_lagged_unrolled)
  shd(dag_gt, dag_full)
  precision(dag_gt, dag_lagged_unrolled)
  recall(dag_gt, dag_lagged_unrolled)
  precision(dag_gt, dag_full)
  recall(dag_gt, dag_full)
  aid(dag_gt, dag_lagged_unrolled, "ancestor")
  aid(dag_gt, dag_full, "ancestor")

  lagged_fit <- parameterize_dag_from_data(
    df = df_lagged,
    dag = dag_lagged,
    var_G = var_G,
    var_X = var_X,
    graph_type = "lagged"
  )
  full_fit <- parameterize_dag_from_data(
    df = df_full,
    dag = dag_full,
    var_G = var_G,
    var_X = var_X,
    graph_type = "full_time"
  )

  # Lagged ATE: do(X1_t1 = high) - do(X1_t1 = low) on an outcome at t2.
  lag_x <- rownames(lagged_fit$matrices$B)
  lag_do <- paste0(lag_x[1], "_t1")
  lag_y <- paste0(lag_x[min(2, length(lag_x))], "_t", min(2, dataset_object$`500`$params$T))

  lag_low <- calculate_intervention_distribution(
    N = 10000,
    T = dataset_object$`500`$params$T,
    H = lagged_fit$matrices$H,
    B = lagged_fit$matrices$B,
    Gamma = lagged_fit$matrices$Gamma,
    A = lagged_fit$matrices$A,
    L = lagged_fit$matrices$L,
    intervention = setNames(list(-2), lag_do),
    type_X = var_X,
    type_G = var_G,
    dist_params = lagged_fit$dist_params,
    seed = 123
  )
  lag_high <- calculate_intervention_distribution(
    N = 10000,
    T = dataset_object$`500`$params$T,
    H = lagged_fit$matrices$H,
    B = lagged_fit$matrices$B,
    Gamma = lagged_fit$matrices$Gamma,
    A = lagged_fit$matrices$A,
    L = lagged_fit$matrices$L,
    intervention = setNames(list(2), lag_do),
    type_X = var_X,
    type_G = var_G,
    dist_params = lagged_fit$dist_params,
    seed = 123
  )
  ate_lagged <- mean(lag_high[[lag_y]]) - mean(lag_low[[lag_y]])

  # Full-time ATE: B already contains all X*_t* nodes, so T = 1 and A = NULL.
  full_x <- rownames(full_fit$matrices$B)
  full_do <- full_x[1]
  full_y <- full_x[min(2, length(full_x))]

  full_low <- calculate_intervention_distribution(
    N = 10000,
    T = 1,
    H = full_fit$matrices$H,
    B = full_fit$matrices$B,
    Gamma = full_fit$matrices$Gamma,
    A = NULL,
    L = 0L,
    intervention = setNames(list(-2), full_do),
    type_X = var_X,
    type_G = var_G,
    dist_params = full_fit$dist_params,
    seed = 123
  )
  full_high <- calculate_intervention_distribution(
    N = 10000,
    T = 1,
    H = full_fit$matrices$H,
    B = full_fit$matrices$B,
    Gamma = full_fit$matrices$Gamma,
    A = NULL,
    L = 0L,
    intervention = setNames(list(2), full_do),
    type_X = var_X,
    type_G = var_G,
    dist_params = full_fit$dist_params,
    seed = 123
  )
  ate_full <- mean(full_high[[full_y]]) - mean(full_low[[full_y]])

  print(data.frame(
    graph = c("lagged", "full_time"),
    intervention = c(lag_do, full_do),
    outcome = c(lag_y, full_y),
    ate_high_minus_low = c(ate_lagged, ate_full)
  ))
}
