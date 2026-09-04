#' Factor Analysis Clustering
#'
#' Assigns respondents to clusters based on factor scores from
#' \code{psych::fa}. Multiple solutions (different numbers of factors)
#' can be requested in a single call.
#'
#' If the first column is an ID variable (e.g. \code{resp_id}), it is
#' excluded from the analysis and appended to \code{Cluster_Assignments}.
#'
#' Column names in \code{Cluster_Assignments} follow the pattern
#' \code{facts\{CorType\}\{Fm\}_k}, e.g. \code{factsCorPa_3}, so that
#' downstream processes that parse \code{NAME_N} continue to work.
#'
#' @param data A data frame or matrix. If the first column is an ID
#'   (e.g. respondent ID), it is preserved and returned with the
#'   cluster assignments. All remaining columns are used as inputs.
#' @param num_solutions Integer vector of the numbers of clusters
#'   (factors) to extract. Default is \code{3:8}.
#' @param cor_type Character. Correlation type passed to \code{psych::fa}.
#'   One of \code{"cor"}, \code{"poly"}, or \code{"tet"}. Default \code{"cor"}.
#' @param rotate Character. Rotation method. Default \code{"varimax"}.
#' @param fm Character. Factoring method. Default \code{"pa"}
#'   (principal axis). See \code{?psych::fa} for options.
#' @param scoring Character. How to assign a cluster from the factor scores.
#'   One of \code{"max"} (maximum score), \code{"min"} (minimum score),
#'   or \code{"maxabs"} (maximum absolute score). Default \code{"max"}.
#' @param standardize Logical. If \code{TRUE}, input variables (excluding
#'   any ID column) are z-scored before analysis. Default \code{FALSE}.
#' @param max_iter Integer. Maximum iterations for \code{psych::fa}.
#'   Default \code{50}. Increase if you see "maximum iteration exceeded".
#' @param check_nonnegative Logical. If \code{TRUE}, stop when any value
#'   is negative. Default \code{TRUE}.
#' @param id_col Logical or character. If \code{TRUE} (default), treat the
#'   first column as an ID. If a character string, treat the column with
#'   that name as the ID. If \code{FALSE}, no ID column is assumed.
#'
#' @return A list with:
#' \describe{
#'   \item{Cluster_Assignments}{Data frame of cluster memberships
#'     (ID column first, if present, then one column per solution).
#'     Column names look like \code{factsCorPa_3}.}
#'   \item{Cluster_Sizes}{Long data frame of cluster sizes.}
#'   \item{Rotated_Component_Matrix}{Named list of loadings matrices,
#'     one per requested solution.}
#'   \item{Correlation_Matrix}{Correlation matrix used in the analysis.}
#'   \item{Input_Value_Frequencies}{Frequency table of the (converted)
#'     input values.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cluster_fa(
#'   data          = for_clust,
#'   num_solutions = 3:8,
#'   cor_type      = "cor",
#'   rotate        = "varimax",
#'   fm            = "pa",
#'   scoring       = "max",
#'   standardize   = TRUE,
#'   max_iter      = 100
#' )
#' names(result$Cluster_Assignments)
#' # e.g. "resp_id" "factsCorPa_3" "factsCorPa_4" ...
#' }
#'
#' @importFrom psych fa polychoric tetrachoric
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_fa <- function(data,
                       num_solutions = 3:8,
                       cor_type = c("cor", "poly", "tet"),
                       rotate = "varimax",
                       fm = "pa",
                       scoring = c("max", "min", "maxabs"),
                       standardize = FALSE,
                       max_iter = 50L,
                       check_nonnegative = TRUE,
                       id_col = TRUE) {
  
  cor_type <- match.arg(cor_type)
  scoring  <- match.arg(scoring)
  max_iter <- as.integer(max_iter)
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two analysis variables (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  if (check_nonnegative && any(data < 0, na.rm = TRUE)) {
    stop("All input values must be non-negative.")
  }
  
  # ---- optional standardization ----
  if (isTRUE(standardize)) {
    data <- as.data.frame(scale(data))
  }
  
  # ---- correlation-type validation ----
  n_unique <- vapply(data, function(x) length(unique(stats::na.omit(x))),
                     integer(1))
  
  if (cor_type == "tet" && any(n_unique > 2L)) {
    stop("Tetrachoric correlations ('tet') require dichotomous (binary) data.\n",
         "For continuous data use 'cor'. For Likert / ordinal use 'poly'.")
  }
  if (cor_type == "poly" && any(n_unique > 15L)) {
    warning("Polychoric correlations work best with ordinal data ",
            "that have relatively few categories.", call. = FALSE)
  }
  
  # ---- solution prefix: factsCorPa, factsPolyMinres, etc. ----
  # Capitalise first letter of cor_type and fm for camelCase
  cap1 <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  sol_prefix <- paste0("facts", cap1(cor_type), cap1(fm))
  # e.g. "factsCorPa", "factsPolyPa", "factsTetMl"
  sol_names  <- paste0(sol_prefix, "_", num_solutions)
  
  # ---- scoring helper ----
  get_cluster <- function(scores, method) {
    switch(method,
           max    = max.col(scores),
           min    = max.col(-scores),
           maxabs = max.col(abs(scores)))
  }
  
  # ---- main clustering + loadings for every solution ----
  fa_results <- purrr::map(num_solutions, function(k) {
    psych::fa(
      r        = data,
      nfactors = k,
      n.iter   = 1,
      rotate   = rotate,
      fm       = fm,
      SMC      = FALSE,
      cor      = cor_type,
      max.iter = max_iter
    )
  })
  names(fa_results) <- sol_names
  
  clusters <- purrr::map(fa_results, function(fa_res) {
    get_cluster(fa_res$scores, scoring)
  }) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  # Prepend ID if present
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- rotated loadings: one matrix per solution ----
  loadings_list <- purrr::map(fa_results, function(fa_res) {
    k <- ncol(fa_res$loadings)
    mat <- as.data.frame(unclass(fa_res$loadings))
    colnames(mat) <- paste0("RC", seq_len(k))
    mat$Variable <- rownames(mat)
    mat[, c("Variable", paste0("RC", seq_len(k)))]
  })
  names(loadings_list) <- sol_names
  
  # ---- cluster sizes ----
  assign_cols <- sol_names
  cluster_sizes <- lapply(assign_cols, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- correlation matrix ----
  corr_mat <- switch(
    cor_type,
    cor  = stats::cor(data, use = "pairwise.complete.obs"),
    poly = psych::polychoric(data)$rho,
    tet  = psych::tetrachoric(data)$rho
  )
  corr_mat <- as.data.frame(corr_mat)
  corr_mat$Variable <- rownames(corr_mat)
  corr_mat <- corr_mat[, c("Variable", setdiff(names(corr_mat), "Variable"))]
  
  # ---- value frequencies ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  # ---- return ----
  list(
    Cluster_Assignments      = clusters,
    Cluster_Sizes            = cluster_sizes,
    Rotated_Component_Matrix = loadings_list,
    Correlation_Matrix       = corr_mat,
    Input_Value_Frequencies  = value_freqs
  )
}

#' NMF Clustering
#'
#' Assigns respondents to clusters via non-negative matrix factorization
#' using the \pkg{NMF} package. Multiple ranks can be requested in one call.
#'
#' If the first column is an ID variable (e.g. \code{resp_id}), it is
#' excluded from the analysis and appended to \code{Cluster_Assignments}.
#'
#' Column names follow the pattern \code{nmf\{Method\}_k}, e.g.
#' \code{nmfLee_3}, \code{nmfBrunet_4}.
#'
#' @param data A data frame or matrix. First column may be an ID.
#'   Remaining columns must be non-negative.
#' @param num_solutions Integer vector of ranks. Default \code{3:8}.
#' @param method Character. NMF algorithm (e.g. \code{"lee"}, \code{"brunet"},
#'   \code{"offset"}, \code{"nsNMF"}, \code{"snmf/r"}, \code{"snmf/l"},
#'   \code{"ls-nmf"}). Default \code{"lee"}.
#' @param nrun Integer. Number of random runs per rank. Default \code{30}.
#' @param seed Integer. Random seed. Default \code{123}.
#' @param scoring Character. \code{"max"} or \code{"min"} basis value.
#'   Default \code{"max"}.
#' @param check_nonnegative Logical. Default \code{TRUE}.
#' @param id_col Logical or character. If \code{TRUE}, first column is ID.
#'   If a string, that column name is the ID. If \code{FALSE}, no ID.
#'
#' @return A list with Cluster_Assignments, Cluster_Sizes, Basis_W,
#'   Coefficient_H, Fitted, NMF_Method, Nrun, and Input_Value_Frequencies.
#'
#' @importFrom NMF nmf basis coef fitted
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_nmf <- function(data,
                        num_solutions = 3:8,
                        method = "lee",
                        nrun = 30L,
                        seed = 123L,
                        scoring = c("max", "min"),
                        check_nonnegative = TRUE,
                        id_col = TRUE) {
  
  scoring <- match.arg(scoring)
  nrun    <- as.integer(nrun)
  seed    <- as.integer(seed)
  
  # NMF must be attached; NMF:: alone often triggers path.package() errors
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("Package 'NMF' is required. Install it with install.packages(\"NMF\").")
  }
  # Force full load so algorithm registry and path.package work
  suppressPackageStartupMessages(library(NMF))
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 1L) {
    stop("Please supply at least one analysis variable (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  if (check_nonnegative && any(data < 0, na.rm = TRUE)) {
    stop("NMF requires all input values to be non-negative.")
  }
  
  if (anyNA(data)) {
    warning("Missing values detected; they will be replaced with 0.", call. = FALSE)
    data[is.na(data)] <- 0
  }
  
  inputs_mat <- as.matrix(data)
  
  # ---- solution names: nmfLee_3, nmfBrunet_4, ... ----
  method_clean <- method
  method_clean <- gsub("/", "", method_clean)
  method_clean <- gsub("-", "", method_clean)
  method_clean <- gsub("nmf", "", method_clean, ignore.case = TRUE)
  cap1 <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  sol_prefix <- paste0("nmf", cap1(method_clean))
  sol_names  <- paste0(sol_prefix, "_", num_solutions)
  
  # ---- scoring helper ----
  get_cluster <- function(basis_mat, method) {
    switch(method,
           max = max.col(basis_mat),
           min = max.col(-basis_mat))
  }
  
  # ---- fit NMF for each rank ----
  set.seed(seed)
  
  nmf_results <- purrr::map(num_solutions, function(k) {
    if (identical(method, "ls-nmf")) {
      NMF::nmf(
        inputs_mat,
        rank   = k,
        method = method,
        nrun   = nrun,
        weight = rep(1, nrow(inputs_mat)),
        seed   = seed
      )
    } else {
      NMF::nmf(
        inputs_mat,
        rank   = k,
        method = method,
        nrun   = nrun,
        seed   = seed
      )
    }
  })
  names(nmf_results) <- sol_names
  
  # ---- cluster assignments ----
  clusters <- purrr::map(nmf_results, function(fit) {
    get_cluster(NMF::basis(fit), scoring)
  }) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- cluster sizes ----
  cluster_sizes <- lapply(sol_names, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- Basis_W, Coefficient_H, Fitted (one matrix per rank) ----
  basis_list  <- list()
  coef_list   <- list()
  fitted_list <- list()
  
  for (i in seq_along(nmf_results)) {
    k       <- num_solutions[i]
    nmf_obj <- nmf_results[[i]]
    key     <- sol_names[i]
    
    # Basis (W) — respondents x rank
    bmat <- as.data.frame(NMF::basis(nmf_obj))
    colnames(bmat) <- paste0("Basis_", seq_len(k))
    bmat$Case <- seq_len(nrow(bmat))
    bmat <- bmat[, c("Case", paste0("Basis_", seq_len(k)))]
    if (!is.null(id_vector)) {
      bmat <- dplyr::bind_cols(!!id_name := id_vector, bmat[, -1, drop = FALSE])
    }
    basis_list[[key]] <- bmat
    
    # Coefficient (H) — rank x variables (NMF stores H as rank x features)
    cmat <- as.data.frame(NMF::coef(nmf_obj))
    # coef is typically rank x n_features; transpose for variables as rows if needed
    # Keep as returned by NMF; name columns by component
    if (nrow(cmat) == k) {
      colnames(cmat) <- colnames(inputs_mat)
      cmat$Component <- paste0("Comp_", seq_len(k))
      cmat <- cmat[, c("Component", colnames(inputs_mat))]
    } else {
      colnames(cmat) <- paste0("Coef_", seq_len(ncol(cmat)))
      cmat$Variable <- rownames(cmat)
      cmat <- cmat[, c("Variable", setdiff(names(cmat), "Variable"))]
    }
    coef_list[[key]] <- cmat
    
    # Fitted
    fmat <- as.data.frame(NMF::fitted(nmf_obj))
    colnames(fmat) <- colnames(inputs_mat)
    fmat$Case <- seq_len(nrow(fmat))
    fmat <- fmat[, c("Case", colnames(inputs_mat))]
    if (!is.null(id_vector)) {
      fmat <- dplyr::bind_cols(!!id_name := id_vector, fmat[, -1, drop = FALSE])
    }
    fitted_list[[key]] <- fmat
  }
  
  # ---- value frequencies ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  # ---- return ----
  list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Basis_W                 = basis_list,
    Coefficient_H           = coef_list,
    Fitted                  = fitted_list,
    NMF_Method              = method,
    Nrun                    = nrun,
    Input_Value_Frequencies = value_freqs,
    Note                    = paste0(
      "Basis_W, Coefficient_H and Fitted are lists named by solution ",
      "(e.g. Basis_W$", sol_names[1], ")."
    )
  )
}

#' Archetypal Clustering
#'
#' Assigns respondents to clusters via archetypal analysis using the
#' \pkg{archetypes} package. Multiple numbers of archetypes can be
#' requested in a single call.
#'
#' If the first column is an ID variable (e.g. \code{resp_id}), it is
#' excluded from the analysis and appended to \code{Cluster_Assignments}.
#'
#' Column names in \code{Cluster_Assignments} follow the pattern
#' \code{arch_k}, e.g. \code{arch_3}, \code{arch_4}, so that downstream
#' processes that parse \code{NAME_N} continue to work.
#'
#' @param data A data frame or matrix. If the first column is an ID
#'   (e.g. respondent ID), it is preserved and returned with the
#'   cluster assignments. All remaining columns are used as inputs.
#' @param num_solutions Integer vector of the numbers of archetypes
#'   to extract. Default is \code{3:8}.
#' @param scoring Character. How to assign a cluster from the alpha
#'   (soft membership) matrix. One of \code{"max"} (maximum alpha) or
#'   \code{"min"} (minimum alpha). Default \code{"max"}.
#' @param seed Integer. Random seed for reproducibility. Default \code{123}.
#' @param id_col Logical or character. If \code{TRUE} (default), treat the
#'   first column as an ID. If a character string, treat the column with
#'   that name as the ID. If \code{FALSE}, no ID column is assumed.
#'
#' @return A list with:
#' \describe{
#'   \item{Cluster_Assignments}{Data frame of cluster memberships
#'     (ID column first, if present, then one column per solution).
#'     Column names look like \code{arch_3}.}
#'   \item{Cluster_Sizes}{Long data frame of cluster sizes.}
#'   \item{Alphas}{Named list of soft-membership matrices, one per solution.}
#'   \item{Archetypes}{Named list of archetype profile matrices, one per solution.}
#'   \item{RSS}{Named list of residual sum of squares, one per solution.}
#'   \item{Input_Value_Frequencies}{Frequency table of the (converted)
#'     input values.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cluster_arch(
#'   data          = for_clust,   # resp_id in first column
#'   num_solutions = 3:5,
#'   scoring       = "max",
#'   seed          = 123
#' )
#' names(result$Cluster_Assignments)
#' # e.g. "resp_id" "arch_3" "arch_4" "arch_5"
#' result$Archetypes$arch_3
#' result$Alphas$arch_4
#' }
#'
#' @importFrom archetypes archetypes
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_arch <- function(data,
                         num_solutions = 3:8,
                         scoring = c("max", "min"),
                         seed = 123L,
                         id_col = TRUE) {
  
  scoring <- match.arg(scoring)
  seed    <- as.integer(seed)
  
  if (!requireNamespace("archetypes", quietly = TRUE)) {
    stop("Package 'archetypes' is required. Install it with install.packages(\"archetypes\").")
  }
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two analysis variables (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  # ---- value frequencies (before matrix conversion) ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  inputs_mat <- as.matrix(data)
  
  # ---- solution names: arch_3, arch_4, ... ----
  # No method subtypes in base archetypes::archetypes
  sol_names <- paste0("arch_", num_solutions)
  
  # ---- scoring helper ----
  get_cluster <- function(alphas, method) {
    switch(method,
           max = max.col(alphas),
           min = max.col(-alphas))
  }
  
  # ---- fit archetypal analysis for each k ----
  set.seed(seed)
  
  arch_results <- purrr::map(num_solutions, function(k) {
    archetypes::archetypes(inputs_mat, k)
  })
  names(arch_results) <- sol_names
  
  # ---- cluster assignments ----
  clusters <- purrr::map(arch_results, function(arch_obj) {
    get_cluster(arch_obj$alphas, scoring)
  }) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- cluster sizes ----
  cluster_sizes <- lapply(sol_names, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- Alphas, Archetypes, RSS (one entry per solution) ----
  alphas_list     <- list()
  archetypes_list <- list()
  rss_list        <- list()
  
  for (i in seq_along(arch_results)) {
    k        <- num_solutions[i]
    arch_obj <- arch_results[[i]]
    key      <- sol_names[i]
    
    # Alphas (soft memberships): respondents x k
    amat <- as.data.frame(arch_obj$alphas)
    colnames(amat) <- paste0("Alpha_", seq_len(k))
    if (!is.null(id_vector)) {
      amat <- dplyr::bind_cols(!!id_name := id_vector, amat)
    } else {
      amat$Case <- seq_len(nrow(amat))
      amat <- amat[, c("Case", paste0("Alpha_", seq_len(k)))]
    }
    alphas_list[[key]] <- amat
    
    # Archetypes (prototype profiles): k x variables
    archmat <- as.data.frame(arch_obj$archetypes)
    colnames(archmat) <- colnames(inputs_mat)
    archmat$Archetype <- paste0("AT_", seq_len(k))
    archmat <- archmat[, c("Archetype", colnames(inputs_mat))]
    archetypes_list[[key]] <- archmat
    
    # Residual sum of squares
    rss_list[[key]] <- arch_obj$rss
  }
  
  # ---- return ----
  list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Alphas                  = alphas_list,
    Archetypes              = archetypes_list,
    RSS                     = rss_list,
    Input_Value_Frequencies = value_freqs,
    Note                    = paste0(
      "Alphas and Archetypes are lists named by solution ",
      "(e.g. Alphas$", sol_names[1], ", Archetypes$", sol_names[1], ")."
    )
  )
}

#' Model-Based Clustering (mclust)
#'
#' Assigns respondents to clusters via Gaussian mixture models using
#' \pkg{mclust}. Multiple numbers of components (\code{G}) can be
#' requested in a single call.
#'
#' If the first column is an ID variable (e.g. \code{resp_id}), it is
#' excluded from the analysis and appended to \code{Cluster_Assignments}.
#'
#' Column names follow the pattern \code{mod\{Model\}_k}, e.g.
#' \code{modVvv_3}, \code{modAll_4}, so that downstream processes that
#' parse \code{NAME_N} continue to work.
#'
#' @section Covariance models (\code{model_names}):
#' Each mixture component is a multivariate normal. The three-letter
#' \code{modelNames} code in \pkg{mclust} controls the covariance
#' structure:
#' \describe{
#'   \item{1st letter – Volume}{\code{E} = equal across clusters;
#'     \code{V} = variable.}
#'   \item{2nd letter – Shape}{\code{E} = equal; \code{V} = variable;
#'     \code{I} = spherical (identity).}
#'   \item{3rd letter – Orientation}{\code{E} = equal; \code{V} = variable;
#'     \code{I} = axis-aligned.}
#' }
#' Common options include \code{"EII"} (equal spherical), \code{"VII"}
#' (variable spherical), \code{"EEE"} (equal ellipsoids), \code{"VVV"}
#' (unconstrained), and \code{"ALL"} (try a suite of models and pick the
#' best BIC for each \code{G}). Different models can produce different
#' assignments when clusters differ in size, shape, or orientation.
#'
#' @param data A data frame or matrix. If the first column is an ID
#'   (e.g. respondent ID), it is preserved and returned with the
#'   cluster assignments. All remaining columns are used as inputs.
#' @param num_solutions Integer vector of the numbers of clusters
#'   (\code{G}) to fit. Default is \code{3:8}.
#' @param model_names Character. Covariance structure passed to
#'   \code{mclust::Mclust} as \code{modelNames}. Use \code{"ALL"}
#'   (default) to let mclust choose the best model for each \code{G},
#'   or a specific code such as \code{"EII"}, \code{"VII"}, \code{"EEI"},
#'   \code{"VVV"}, \code{"VEV"}, etc. See the section above.
#' @param id_col Logical or character. If \code{TRUE} (default), treat the
#'   first column as an ID. If a character string, treat the column with
#'   that name as the ID. If \code{FALSE}, no ID column is assumed.
#'
#' @return A list with:
#' \describe{
#'   \item{Cluster_Assignments}{Data frame of hard cluster memberships
#'     (ID column first, if present). Names look like \code{modVvv_3}.}
#'   \item{Cluster_Sizes}{Long data frame of cluster sizes.}
#'   \item{Soft_Probabilities}{Named list of posterior probability
#'     matrices (\code{z}), one per solution.}
#'   \item{BIC}{Named list of BIC values, one per solution.}
#'   \item{LogLikelihood}{Named list of log-likelihoods, one per solution.}
#'   \item{Model_Used}{Named list of the selected mclust model name
#'     for each solution (especially useful when \code{model_names = "ALL"}).}
#'   \item{Input_Value_Frequencies}{Frequency table of the (converted)
#'     input values.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cluster_mod(
#'   data          = for_clust,
#'   num_solutions = 3:5,
#'   model_names   = "ALL"
#' )
#' names(result$Cluster_Assignments)
#' result$Model_Used
#' result$Soft_Probabilities$modAll_3
#'
#' # Force a specific covariance structure
#' result_vvv <- cluster_mod(for_clust, num_solutions = 3:4, model_names = "VVV")
#' names(result_vvv$Cluster_Assignments)  # "resp_id" "modVvv_3" "modVvv_4"
#' }
#'
#' @importFrom mclust Mclust
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_mod <- function(data,
                        num_solutions = 3:8,
                        model_names = "ALL",
                        id_col = TRUE) {
  
  # mclust must be attached: Mclust() calls mclustBIC() without namespace
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required. Install it with install.packages(\"mclust\").")
  }
  suppressPackageStartupMessages(library(mclust))
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two analysis variables (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  # ---- value frequencies ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  inputs_mat <- as.matrix(data)
  
  # ---- solution names: modAll_3, modVvv_4, ... ----
  cap1 <- function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  }
  model_tag <- if (toupper(model_names) == "ALL") {
    "All"
  } else {
    cap1(model_names)
  }
  sol_prefix <- paste0("mod", model_tag)
  sol_names  <- paste0(sol_prefix, "_", num_solutions)
  
  # ---- fit Mclust for each G ----
  mclust_results <- purrr::map(num_solutions, function(g) {
    if (toupper(model_names) == "ALL") {
      mclust::Mclust(inputs_mat, G = g)
    } else {
      mclust::Mclust(inputs_mat, G = g, modelNames = model_names)
    }
  })
  names(mclust_results) <- sol_names
  
  failed <- vapply(mclust_results, is.null, logical(1))
  if (any(failed)) {
    stop("Mclust failed for G = ",
         paste(num_solutions[failed], collapse = ", "),
         ". Try a different model_names value or check the data.")
  }
  
  # ---- hard cluster assignments ----
  clusters <- purrr::map(mclust_results, function(obj) {
    obj$classification
  }) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- cluster sizes ----
  cluster_sizes <- lapply(sol_names, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- soft probabilities, BIC, loglik, model used ----
  z_list      <- list()
  bic_list    <- list()
  loglik_list <- list()
  model_list  <- list()
  
  for (i in seq_along(mclust_results)) {
    g   <- num_solutions[i]
    obj <- mclust_results[[i]]
    key <- sol_names[i]
    
    zmat <- as.data.frame(obj$z)
    colnames(zmat) <- paste0("Prob_", seq_len(g))
    if (!is.null(id_vector)) {
      zmat <- dplyr::bind_cols(!!id_name := id_vector, zmat)
    } else {
      zmat$Case <- seq_len(nrow(zmat))
      zmat <- zmat[, c("Case", paste0("Prob_", seq_len(g)))]
    }
    z_list[[key]] <- zmat
    
    bic_list[[key]]    <- obj$bic
    loglik_list[[key]] <- obj$loglik
    model_list[[key]]  <- obj$modelName
  }
  
  # ---- return ----
  list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Soft_Probabilities      = z_list,
    BIC                     = bic_list,
    LogLikelihood           = loglik_list,
    Model_Used              = model_list,
    Input_Value_Frequencies = value_freqs,
    Note                    = paste0(
      "Soft_Probabilities is a list named by solution ",
      "(e.g. Soft_Probabilities$", sol_names[1], "). ",
      "When model_names = \"ALL\", Model_Used shows which covariance ",
      "model mclust selected for each G. Different covariance models ",
      "can yield different assignments."
    )
  )
}

#' Canonical Correlation + K-means Clustering
#'
#' Runs canonical correlation analysis (CCA) between a primary (X) and
#' secondary (Y) variable set, keeps the Wilks-significant X canonical
#' dimensions, then clusters respondents with k-means on those scores.
#'
#' If the first column of \code{primary} is an ID variable
#' (e.g. \code{resp_id}), it is excluded from the analysis and appended
#' to \code{Cluster_Assignments}. \code{secondary} is assumed to be in
#' the same row order and should not contain an ID column (or use
#' matching IDs only for alignment checks outside this function).
#'
#' Column names follow the pattern \code{kmcc_k}, e.g. \code{kmcc_3},
#' so that downstream processes that parse \code{NAME_N} continue to work.
#'
#' @param primary A data frame or matrix of primary (X) variables.
#'   Optional ID in the first column when \code{id_col = TRUE}.
#' @param secondary A data frame or matrix of secondary (Y) variables.
#'   Must have the same number of rows as \code{primary} (after ID removal).
#' @param num_solutions Integer vector of k-means cluster counts.
#'   Default is \code{3:8}.
#' @param sig_level Numeric. Significance level for the Wilks test when
#'   choosing how many canonical dimensions to retain. Default \code{0.05}.
#' @param nstart Integer. Number of random k-means starts. Default \code{25}.
#' @param seed Integer. Random seed for reproducibility. Default \code{123}.
#' @param id_col Logical or character. If \code{TRUE} (default), treat the
#'   first column of \code{primary} as an ID. If a character string, treat
#'   that column name in \code{primary} as the ID. If \code{FALSE}, no ID.
#'
#' @return A list with:
#' \describe{
#'   \item{Cluster_Assignments}{Data frame of cluster memberships
#'     (ID first, if present). Names look like \code{kmcc_3}.}
#'   \item{Cluster_Sizes}{Long data frame of cluster sizes.}
#'   \item{Canonical_Scores_X}{Significant X canonical score matrix.}
#'   \item{Wilks_Test}{Wilks lambda test table for each dimension.}
#'   \item{Num_Significant_Dims}{Number of dimensions retained.}
#'   \item{Canonical_Correlations}{Vector of canonical correlations.}
#'   \item{Input_Value_Frequencies}{Frequency table for the primary set.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cluster_cancor(
#'   primary       = primary_df,    # resp_id + X vars
#'   secondary     = secondary_df,  # Y vars only, same row order
#'   num_solutions = 3:5,
#'   sig_level     = 0.05,
#'   seed          = 123
#' )
#' names(result$Cluster_Assignments)
#' result$Wilks_Test
#' result$Num_Significant_Dims
#' }
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_cancor <- function(primary,
                           secondary,
                           num_solutions = 3:8,
                           sig_level = 0.05,
                           nstart = 25L,
                           seed = 123L,
                           id_col = TRUE) {
  
  nstart <- as.integer(nstart)
  seed   <- as.integer(seed)
  
  # ---- Wilks test helper ----
  wilks_test <- function(rho, n, p, q) {
    k <- length(rho)
    out <- data.frame(
      Dimension      = integer(k),
      Canonical_R    = numeric(k),
      Wilks_Lambda   = numeric(k),
      Chi_Sq         = numeric(k),
      df             = numeric(k),
      p_value        = numeric(k)
    )
    for (i in seq_len(k)) {
      lam  <- prod(1 - rho[i:k]^2)
      chi  <- -(n - 1 - (p + q + 1) / 2) * log(lam)
      deg  <- (p - i + 1) * (q - i + 1)
      pval <- stats::pchisq(chi, deg, lower.tail = FALSE)
      out[i, ] <- list(i, rho[i], lam, chi, deg, pval)
    }
    out
  }
  
  to_numeric <- function(df) {
    as.data.frame(lapply(df, function(x) {
      if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
    }))
  }
  
  # ---- checks ----
  if (!is.data.frame(primary) && !is.matrix(primary)) {
    stop("`primary` must be a data frame or matrix.")
  }
  if (!is.data.frame(secondary) && !is.matrix(secondary)) {
    stop("`secondary` must be a data frame or matrix.")
  }
  primary   <- as.data.frame(primary)
  secondary <- as.data.frame(secondary)
  
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column on primary ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    if (ncol(primary) < 2L) {
      stop("`primary` needs an ID column plus at least one variable.")
    }
    id_name   <- names(primary)[1]
    id_vector <- primary[[1]]
    primary   <- primary[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(primary)) {
      stop("Column '", id_col, "' not found in `primary`.")
    }
    id_name   <- id_col
    id_vector <- primary[[id_col]]
    primary   <- primary[, setdiff(names(primary), id_col), drop = FALSE]
  }
  
  if (ncol(primary) < 2L) {
    stop("Please supply at least two primary (X) variables.")
  }
  if (ncol(secondary) < 2L) {
    stop("Please supply at least two secondary (Y) variables.")
  }
  
  primary   <- to_numeric(primary)
  secondary <- to_numeric(secondary)
  
  if (nrow(primary) != nrow(secondary)) {
    stop("`primary` and `secondary` must have the same number of rows ",
         "(after removing any ID column from primary).")
  }
  
  # ---- value frequencies (primary set) ----
  value_freqs <- lapply(names(primary), function(var) {
    tab <- as.data.frame(table(primary[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  X <- as.matrix(primary)
  Y <- as.matrix(secondary)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # ---- CCA ----
  cca <- stats::cancor(X, Y)
  wilks_tab <- wilks_test(cca$cor, n = n, p = p, q = q)
  
  num_sig <- sum(wilks_tab$p_value < sig_level, na.rm = TRUE)
  if (num_sig < 1L) {
    warning("No canonical dimensions were significant at p < ", sig_level,
            ". Using the first dimension.", call. = FALSE)
    num_sig <- 1L
  }
  
  # Significant X canonical scores (centered X %*% xcoef)
  x_scores <- as.data.frame(
    scale(X, center = TRUE, scale = FALSE) %*% cca$xcoef[, seq_len(num_sig), drop = FALSE]
  )
  colnames(x_scores) <- paste0("Xcan_", seq_len(num_sig))
  
  if (!is.null(id_vector)) {
    x_scores_out <- dplyr::bind_cols(!!id_name := id_vector, x_scores)
  } else {
    x_scores_out <- x_scores
  }
  
  # ---- solution names: kmcc_3, kmcc_4, ... ----
  # No CCA method subtypes in base stats::cancor
  sol_names <- paste0("kmcc_", num_solutions)
  
  # ---- k-means on canonical scores ----
  set.seed(seed)
  
  clusters <- purrr::map(num_solutions, function(k) {
    stats::kmeans(x_scores, centers = k, nstart = nstart)$cluster
  }) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- cluster sizes ----
  cluster_sizes <- lapply(sol_names, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- return ----
  list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Canonical_Scores_X      = x_scores_out,
    Wilks_Test              = wilks_tab,
    Num_Significant_Dims    = num_sig,
    Canonical_Correlations  = cca$cor,
    Input_Value_Frequencies = value_freqs,
    Note                    = paste0(
      "Used ", num_sig,
      " significant canonical dimension(s) from the primary (X) set ",
      "based on Wilks p < ", sig_level, "."
    )
  )
}

#' K-means and Hierarchical Clustering
#'
#' A unified interface for several common clustering procedures used in
#' survey and market-research segmentation. Multiple numbers of clusters
#' can be requested in a single call; results are returned in a consistent
#' list structure suitable for downstream profiling and evaluation.
#'
#' If the first column of \code{data} is an ID variable (e.g. \code{resp_id}),
#' it is excluded from the analysis and prepended to
#' \code{Cluster_Assignments}. Column names always end in \code{_k}
#' (e.g. \code{kmLloyd_3}) so that pipelines that parse \code{NAME_N}
#' continue to work.
#'
#' @section Available methods:
#' \describe{
#'   \item{\code{"kmeans"}}{
#'     Standard k-means with random starting centers (\code{nstart} restarts).
#'     Fast and widely used. Sensitive to scale: consider standardising
#'     variables before calling if they are on very different metrics.
#'   }
#'   \item{\code{"kmeans_seeded"}}{
#'     Hierarchical clustering (Euclidean + Ward.D2) is run first; the
#'     resulting group means are used as starting centers for k-means.
#'     Often more stable than pure random starts and reduces the chance
#'     of poor local minima.
#'   }
#'   \item{\code{"h_ward"}}{
#'     Hierarchical clustering with Euclidean distance and Ward.D2 linkage.
#'     Fully deterministic. Produces a nested sequence of partitions;
#'     \code{num_solutions} simply cuts the same tree at different heights.
#'   }
#'   \item{\code{"h_corr"}}{
#'     Hierarchical clustering using correlation distance between
#'     respondents: \(d_{ij} = 1 - r_{ij}\), where \(r\) is Pearson or
#'     Spearman correlation across variables. Useful when the \emph{shape}
#'     of the response profile matters more than absolute level (e.g. Likert
#'     batteries). Linkage is Ward.D2.
#'   }
#'   \item{\code{"h_gower"}}{
#'     Hierarchical clustering with Gower distance via
#'     \code{cluster::daisy}. After conversion, variables are treated as
#'     numeric/ordinal. Gower is a reasonable default when you may later
#'     extend the tool to mixed types; with purely numeric data it is
#'     similar in spirit to a range-normalised Manhattan distance.
#'     Linkage is Ward.D2.
#'   }
#' }
#'
#' @section K-means algorithms:
#' Passed to \code{stats::kmeans} when \code{method} is \code{"kmeans"}
#' or \code{"kmeans_seeded"}:
#' \describe{
#'   \item{\code{"Lloyd"}}{
#'     Classic alternating assignment and centroid update. Simple and
#'     predictable; a good default for seeded starts.
#'   }
#'   \item{\code{"Hartigan-Wong"}}{
#'     R's default algorithm. Often converges in fewer iterations than
#'     Lloyd and can be slightly more robust.
#'   }
#'   \item{\code{"Forgy"}}{
#'     Similar to Lloyd (batch update of centroids). Historical variant.
#'   }
#'   \item{\code{"MacQueen"}}{
#'     Online-style updates; can behave differently with noisy data.
#'   }
#' }
#' For most survey applications, \code{"Lloyd"} or \code{"Hartigan-Wong"}
#' are sufficient. Algorithm choice usually matters less than scaling,
#' \code{nstart}, and whether you seed from hierarchical clustering.
#'
#' @section Correlation distance (\code{corr_method}):
#' Only used when \code{method = "h_corr"}:
#' \describe{
#'   \item{\code{"pearson"}}{
#'     Linear correlation. Sensitive to magnitude patterns; assumes
#'     roughly interval-scaled variables.
#'   }
#'   \item{\code{"spearman"}}{
#'     Rank correlation. More robust to monotonic non-linearities and
#'     ordinal Likert data.
#'   }
#' }
#' Distance is always \(1 - r\) between respondent profile vectors
#' (rows), not between variables.
#'
#' @param data A data frame or matrix of variables to cluster.
#'   Optional respondent ID in the first column (see \code{id_col}).
#'   Factors and ordered factors are converted to numeric.
#' @param num_solutions Integer vector of cluster counts to produce
#'   (e.g. \code{3:8}). For hierarchical methods these are cuts of the
#'   same tree; for k-means each \(k\) is a separate fit.
#' @param method Clustering procedure. One of \code{"kmeans"},
#'   \code{"kmeans_seeded"}, \code{"h_ward"}, \code{"h_corr"},
#'   \code{"h_gower"}. See \strong{Available methods} above.
#' @param algorithm K-means algorithm for \code{"kmeans"} and
#'   \code{"kmeans_seeded"}: \code{"Lloyd"}, \code{"Hartigan-Wong"},
#'   \code{"Forgy"}, or \code{"MacQueen"}. Ignored for pure hierarchical
#'   methods. Default \code{"Lloyd"}.
#' @param nstart Number of random starts for \code{method = "kmeans"}.
#'   Ignored for other methods. Default \code{50}.
#' @param corr_method \code{"pearson"} or \code{"spearman"} when
#'   \code{method = "h_corr"}. Default \code{"pearson"}.
#' @param seed Integer random seed (affects k-means random starts).
#'   Default \code{123}.
#' @param id_col If \code{TRUE}, the first column of \code{data} is treated
#'   as an ID and excluded from clustering. If a string, the column with
#'   that name is the ID. If \code{FALSE}, all columns are used as inputs.
#'
#' @return A list with at least:
#' \describe{
#'   \item{Cluster_Assignments}{Data frame of memberships (ID column
#'     first if present; then one column per solution).}
#'   \item{Cluster_Sizes}{Long data frame of cluster sizes by solution.}
#'   \item{Method}{The \code{method} value used.}
#'   \item{Input_Value_Frequencies}{Frequency table of converted inputs.}
#' }
#' Additional elements depend on method (\code{Algorithm}, \code{nstart},
#' \code{Correlation_Method}, or a short \code{Note}).
#'
#' @examples
#' \dontrun{
#' # Random-start k-means
#' r1 <- cluster_km(for_clust, num_solutions = 3:5, method = "kmeans",
#'                  algorithm = "Hartigan-Wong", nstart = 50)
#'
#' # Hierarchical seed then k-means (often more stable)
#' r2 <- cluster_km(for_clust, num_solutions = 3:5, method = "kmeans_seeded")
#'
#' # Profile-shape clustering via correlation distance
#' r3 <- cluster_km(for_clust, num_solutions = 3:5, method = "h_corr",
#'                  corr_method = "spearman")
#'
#' names(r2$Cluster_Assignments)
#' r2$Cluster_Sizes
#' }
#'
#' @seealso \code{\link{cluster_fa}}, \code{\link{cluster_nmf}},
#'   \code{\link{cluster_mod}}, \code{\link[stats]{kmeans}},
#'   \code{\link[stats]{hclust}}, \code{\link[cluster]{daisy}}
#'
#' @importFrom purrr map map2
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @importFrom cluster daisy
#' @export
cluster_km <- function(data,
                       num_solutions = 3:8,
                       method = c("kmeans", "kmeans_seeded", "h_ward",
                                  "h_corr", "h_gower"),
                       algorithm = c("Lloyd", "Hartigan-Wong", "Forgy", "MacQueen"),
                       nstart = 50L,
                       corr_method = c("pearson", "spearman"),
                       seed = 123L,
                       id_col = TRUE) {
  
  method      <- match.arg(method)
  algorithm   <- match.arg(algorithm)
  corr_method <- match.arg(corr_method)
  nstart      <- as.integer(nstart)
  seed        <- as.integer(seed)
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two analysis variables (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  # ---- value frequencies ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  inputs_mat <- as.matrix(data)
  
  # ---- solution name prefix by method ----
  algo_tag <- gsub("-", "", algorithm)
  sol_prefix <- switch(
    method,
    kmeans        = paste0("km", algo_tag),
    kmeans_seeded = paste0("kmSeeded", algo_tag),
    h_ward        = "hWard",
    h_corr        = paste0("hCorr", tools::toTitleCase(corr_method)),
    h_gower       = "hGower"
  )
  sol_names <- paste0(sol_prefix, "_", num_solutions)
  
  make_sizes <- function(clusters_df, names_vec) {
    lapply(names_vec, function(nm) {
      tab <- as.data.frame(table(clusters_df[[nm]]))
      colnames(tab) <- c("Cluster", "N")
      tab$Solution <- nm
      tab
    }) |>
      dplyr::bind_rows()
  }
  
  set.seed(seed)
  
  # ============================================================
  # METHOD: K-means (random starts)
  # ============================================================
  if (method == "kmeans") {
    clusters <- purrr::map(num_solutions, function(k) {
      stats::kmeans(
        inputs_mat,
        centers   = k,
        nstart    = nstart,
        algorithm = algorithm,
        iter.max  = 1000
      )$cluster
    }) |>
      dplyr::bind_cols() |>
      stats::setNames(sol_names)
    
    # ============================================================
    # METHOD: K-means seeded by Hierarchical (Ward.D2)
    # ============================================================
  } else if (method == "kmeans_seeded") {
    hc <- stats::hclust(
      stats::dist(inputs_mat, method = "euclidean"),
      method = "ward.D2"
    )
    start_centers <- purrr::map(num_solutions, function(k) {
      membership <- stats::cutree(hc, k = k)
      as.matrix(stats::aggregate(inputs_mat, by = list(membership), FUN = mean)[, -1])
    })
    clusters <- purrr::map2(start_centers, num_solutions, function(centers, k) {
      stats::kmeans(
        inputs_mat,
        centers   = centers,
        algorithm = algorithm,
        iter.max  = 1000
      )$cluster
    }) |>
      dplyr::bind_cols() |>
      stats::setNames(sol_names)
    
    # ============================================================
    # METHOD: Hierarchical – Euclidean + Ward
    # ============================================================
  } else if (method == "h_ward") {
    hc <- stats::hclust(
      stats::dist(inputs_mat, method = "euclidean"),
      method = "ward.D2"
    )
    clusters <- purrr::map(num_solutions, function(k) {
      stats::cutree(hc, k = k)
    }) |>
      dplyr::bind_cols() |>
      stats::setNames(sol_names)
    
    # ============================================================
    # METHOD: Hierarchical – Correlation distance
    # ============================================================
  } else if (method == "h_corr") {
    dist_mat <- stats::as.dist(
      1 - stats::cor(t(inputs_mat), method = corr_method,
                     use = "pairwise.complete.obs")
    )
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    clusters <- purrr::map(num_solutions, function(k) {
      stats::cutree(hc, k = k)
    }) |>
      dplyr::bind_cols() |>
      stats::setNames(sol_names)
    
    # ============================================================
    # METHOD: Hierarchical – Gower distance
    # ============================================================
  } else if (method == "h_gower") {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package 'cluster' is required for Gower distance.")
    }
    dist_mat <- cluster::daisy(data, metric = "gower")
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    clusters <- purrr::map(num_solutions, function(k) {
      stats::cutree(hc, k = k)
    }) |>
      dplyr::bind_cols() |>
      stats::setNames(sol_names)
  }
  
  # ---- prepend ID ----
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  cluster_sizes <- make_sizes(clusters, sol_names)
  
  # ---- return ----
  out <- list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Method                  = method,
    Input_Value_Frequencies = value_freqs
  )
  
  if (method %in% c("kmeans", "kmeans_seeded")) {
    out$Algorithm <- algorithm
  }
  if (method == "kmeans") {
    out$nstart <- nstart
  }
  if (method == "kmeans_seeded") {
    out$Note <- "K-means initialized with centers from hierarchical Ward.D2 clustering."
  }
  if (method == "h_corr") {
    out$Correlation_Method <- corr_method
  }
  if (method == "h_gower") {
    out$Note <- "Gower distance used; variables treated as numeric/ordinal after conversion."
  }
  
  out
}

#' K-Centroids Clustering (flexclust)
#'
#' Runs \code{flexclust::kcca} for one or more numbers of clusters.
#' Supports several distance / centroid families (k-means, k-medians,
#' angle, Jaccard) with controllable initialisation and convergence.
#'
#' If the first column of \code{data} is an ID variable (e.g. \code{resp_id}),
#' it is excluded from the analysis and prepended to
#' \code{Cluster_Assignments}. Column names follow the pattern
#' \code{flex\{Family\}_k}, e.g. \code{flexKmeans_3}, \code{flexEjaccard_4}.
#'
#' @section Families:
#' \describe{
#'   \item{\code{"kmeans"}}{
#'     Euclidean distance with mean centroids. Same objective as classical
#'     k-means; useful when you want flexclust's control options or to
#'     compare with other families in one workflow.
#'   }
#'   \item{\code{"kmedians"}}{
#'     Manhattan distance with median centroids. More robust to outliers
#'     than k-means.
#'   }
#'   \item{\code{"angle"}}{
#'     Cosine / angular distance. Emphasises direction of the profile
#'     rather than magnitude (similar in spirit to correlation-based
#'     clustering of respondent shapes).
#'   }
#'   \item{\code{"jaccard"}, \code{"ejaccard"}}{
#'     Jaccard-type distances for \strong{binary (0/1)} data only.
#'     \code{"ejaccard"} is the extended Jaccard family in flexclust.
#'     The function stops if non-binary values are detected.
#'   }
#' }
#'
#' @section Initialisation (\code{initcent}):
#' \describe{
#'   \item{\code{"random"}}{Random starting centroids.}
#'   \item{\code{"kmeanspp"}}{k-means++ style initialisation (often more stable).}
#' }
#' Other flexclust init options may work if supported by your installed
#' version; the two above are the usual choices.
#'
#' @section Convergence:
#' \describe{
#'   \item{\code{tolerance}}{Relative change criterion for stopping.
#'     Smaller values demand tighter convergence (default \code{0.05}
#'     matches common flexclust usage).}
#'   \item{\code{iter.max}}{Maximum iterations per run (default \code{50}).}
#' }
#'
#' @param data A data frame or matrix. Optional ID in the first column.
#'   Factors / ordered factors are converted to numeric. For Jaccard
#'   families, values must be binary 0/1.
#' @param num_solutions Integer vector of cluster counts (e.g. \code{3:8}).
#' @param family Character. One of \code{"kmeans"}, \code{"kmedians"},
#'   \code{"angle"}, \code{"jaccard"}, \code{"ejaccard"}.
#' @param initcent Character. Centroid initialisation, typically
#'   \code{"random"} or \code{"kmeanspp"}. Default \code{"kmeanspp"}.
#' @param tolerance Numeric. Convergence tolerance. Default \code{0.05}.
#' @param iter.max Integer. Maximum iterations. Default \code{50}.
#' @param seed Integer random seed. Default \code{123}.
#' @param id_col If \code{TRUE}, first column is ID. If a string, that
#'   column name is the ID. If \code{FALSE}, no ID column.
#'
#' @return A list with:
#' \describe{
#'   \item{Cluster_Assignments}{Memberships (ID first if present).}
#'   \item{Cluster_Sizes}{Sizes by solution.}
#'   \item{Centroids}{Named list of centroid tables (variables x clusters),
#'     one per solution.}
#'   \item{Family, Init_Method, Tolerance, Iter_Max}{Run settings.}
#'   \item{Input_Value_Frequencies}{Frequency table of converted inputs.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cluster_flex(
#'   data          = for_clust,
#'   num_solutions = 3:5,
#'   family        = "kmeans",
#'   initcent      = "kmeanspp",
#'   seed          = 123
#' )
#' names(result$Cluster_Assignments)  # e.g. resp_id, flexKmeans_3, ...
#' result$Centroids$flexKmeans_3
#' }
#'
#' @seealso \code{\link{cluster_km}}, \code{\link[flexclust]{kcca}},
#'   \code{\link[flexclust]{kccaFamily}}
#'
#' @importFrom flexclust kcca kccaFamily parameters
#' @importFrom purrr map
#' @importFrom dplyr bind_cols bind_rows select arrange
#' @export
cluster_flex <- function(data,
                         num_solutions = 3:8,
                         family = c("kmeans", "kmedians", "angle",
                                    "jaccard", "ejaccard"),
                         initcent = "kmeanspp",
                         tolerance = 0.05,
                         iter.max = 50L,
                         seed = 123L,
                         id_col = TRUE) {
  
  family   <- match.arg(family)
  iter.max <- as.integer(iter.max)
  seed     <- as.integer(seed)
  
  if (!requireNamespace("flexclust", quietly = TRUE)) {
    stop("Package 'flexclust' is required. Install it with install.packages(\"flexclust\").")
  }
  # Attach so internal class/methods resolve cleanly
  suppressPackageStartupMessages(library(flexclust))
  
  # ---- checks ----
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.")
  }
  data <- as.data.frame(data)
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two columns (ID + variables, or just variables).")
  }
  if (length(num_solutions) < 1L || any(num_solutions < 1L)) {
    stop("`num_solutions` must be a vector of positive integers.")
  }
  num_solutions <- sort(unique(as.integer(num_solutions)))
  
  # ---- handle ID column ----
  id_vector <- NULL
  id_name   <- NULL
  
  if (isTRUE(id_col)) {
    id_name   <- names(data)[1]
    id_vector <- data[[1]]
    data      <- data[, -1, drop = FALSE]
  } else if (is.character(id_col) && length(id_col) == 1L) {
    if (!id_col %in% names(data)) {
      stop("Column '", id_col, "' not found in `data`.")
    }
    id_name   <- id_col
    id_vector <- data[[id_col]]
    data      <- data[, setdiff(names(data), id_col), drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("Please supply at least two analysis variables (after removing any ID column).")
  }
  
  # ---- convert to numeric ----
  data <- as.data.frame(lapply(data, function(x) {
    if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
  }))
  
  # ---- value frequencies ----
  value_freqs <- lapply(names(data), function(var) {
    tab <- as.data.frame(table(data[[var]], useNA = "ifany"))
    colnames(tab) <- c("Value", "N")
    tab$Variable <- var
    tab
  }) |>
    dplyr::bind_rows() |>
    dplyr::select("Variable", "Value", "N") |>
    dplyr::arrange(.data$Variable, .data$Value)
  
  inputs_mat <- as.matrix(data)
  
  # ---- binary check for Jaccard families ----
  if (family %in% c("jaccard", "ejaccard")) {
    unique_vals <- unique(as.vector(inputs_mat))
    unique_vals <- unique_vals[!is.na(unique_vals)]
    if (!all(unique_vals %in% c(0, 1))) {
      stop("Jaccard and ejaccard families require binary (0/1) data.\n",
           "Use family = \"kmeans\", \"kmedians\", or \"angle\" for continuous data.")
    }
  }
  
  # ---- solution names: flexKmeans_3, flexEjaccard_4, ... ----
  cap1 <- function(x) {
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }
  sol_prefix <- paste0("flex", cap1(family))
  sol_names  <- paste0(sol_prefix, "_", num_solutions)
  
  # ---- flexclust control ----
  fc_cont <- new("flexclustControl")
  fc_cont@tolerance <- tolerance
  fc_cont@iter.max  <- iter.max
  fc_cont@verbose   <- 0
  fc_cont@initcent  <- initcent
  
  set.seed(seed)
  
  # ---- fit kcca for each k ----
  kcca_results <- purrr::map(num_solutions, function(k) {
    flexclust::kcca(
      x         = inputs_mat,
      k         = k,
      family    = flexclust::kccaFamily(family),
      control   = fc_cont,
      save.data = FALSE
    )
  })
  names(kcca_results) <- sol_names
  
  # ---- cluster assignments ----
  clusters <- purrr::map(kcca_results, function(obj) obj@cluster) |>
    dplyr::bind_cols() |>
    stats::setNames(sol_names)
  
  if (!is.null(id_vector)) {
    clusters <- dplyr::bind_cols(
      !!id_name := id_vector,
      clusters
    )
  }
  
  # ---- cluster sizes ----
  cluster_sizes <- lapply(sol_names, function(nm) {
    tab <- as.data.frame(table(clusters[[nm]]))
    colnames(tab) <- c("Cluster", "N")
    tab$Solution <- nm
    tab
  }) |>
    dplyr::bind_rows()
  
  # ---- centroids (variables x clusters) for each k ----
  centroids_list <- list()
  for (i in seq_along(kcca_results)) {
    k        <- num_solutions[i]
    key      <- sol_names[i]
    cent_mat <- flexclust::parameters(kcca_results[[i]])
    # parameters(): rows = clusters, columns = variables → transpose
    cent <- as.data.frame(t(cent_mat))
    colnames(cent) <- paste0("Centroid_", seq_len(k))
    cent$Variable <- colnames(inputs_mat)
    cent <- cent[, c("Variable", paste0("Centroid_", seq_len(k)))]
    centroids_list[[key]] <- cent
  }
  
  # ---- return ----
  list(
    Cluster_Assignments     = clusters,
    Cluster_Sizes           = cluster_sizes,
    Centroids               = centroids_list,
    Family                  = family,
    Init_Method             = initcent,
    Tolerance               = tolerance,
    Iter_Max                = iter.max,
    Input_Value_Frequencies = value_freqs,
    Note                    = paste0(
      "Centroids is a list named by solution ",
      "(e.g. Centroids$", sol_names[1], ")."
    )
  )
}

