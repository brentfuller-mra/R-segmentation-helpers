#' Evaluate and Score Multiple Segmentation Solutions
#'
#' Compares multiple cluster solutions across stability, sample balance,
#' attribute differentiation, and predictive discriminant capacity. Computes a
#' unified 0-100 evaluation score and categorization flag per solution.
#'
#' @param assignments A data frame where each column represents a cluster solution, 
#'   OR a list output directly from a \code{cluster_} generation function.
#' @param cluster_vars A data frame or matrix of attitudinal/core driver variables 
#'   used to calculate ANOVA F-statistics and separation profiling.
#' @param disc_vars Optional data frame of secondary demographic or behavioral 
#'   discriminant variables used to run cross-validated classification models.
#' @param min_size_pct Numeric. Minimum acceptable percentage size of a segment 
#'   before triggering a size warning flag. Default \code{5.0}.
#' @param do_lda Logical. If \code{TRUE}, runs a Linear Discriminant Analysis to 
#'   verify cross-validated predictability. Default \code{TRUE}.
#' @param do_migration Logical. If \code{TRUE}, maps transitions across 
#'   consecutive cluster runs (e.g., K=3 to K=4 migrations). Default \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{Solution_Evaluation}{Master tracking frame with scores and Go/No-Go classifications.}
#'   \item{ARI_Matrix}{Adjusted Rand Index consistency matrix between all solutions.}
#'   \item{Size_and_Balance}{Detailed entropy and min/max cluster sizes.}
#'   \item{Differentiation_Metrics}{ANOVA summary variables across cluster inputs.}
#'   \item{Discriminant_Impact}{Eta-squared feature rankings from the LDA fit.}
#'   \item{Migration_Tables}{Transition matrices between sequential cluster cuts.}
#' }
#'
#' @importFrom purrr map_dfr
#' @importFrom dplyr bind_cols bind_rows select arrange mutate left_join case_when desc
#' @importFrom mclust adjustedRandIndex
#' @importFrom MASS lda
#' @export
evaluate_segmentations <- function(assignments,
                                   cluster_vars,
                                   disc_vars = NULL,
                                   min_size_pct = 5.0,
                                   do_lda = TRUE,
                                   do_migration = TRUE) {
  
  # ---- 1. Safe alignment & assignment extraction ----
  if (is.list(assignments) && !is.data.frame(assignments) &&
      !is.null(assignments$Cluster_Assignments)) {
    assignments <- assignments$Cluster_Assignments
  }
  
  if (is.matrix(assignments)) assignments <- as.data.frame(assignments)
  if (is.vector(assignments)) assignments <- data.frame(Solution = assignments)
  
  assignments <- as.data.frame(assignments)
  
  id_candidates <- c("id", "resp_id", "case", "ID", "respondent_id")
  if (length(names(assignments)) > 0 &&
      tolower(names(assignments)[1]) %in% tolower(id_candidates)) {
    assignments <- assignments[, -1, drop = FALSE]
  }
  
  sol_names   <- names(assignments)
  n_solutions <- ncol(assignments)
  
  if (n_solutions < 1) stop("Please supply at least one cluster solution column.")
  
  to_numeric <- function(df) {
    as.data.frame(lapply(df, function(x) {
      if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
    }))
  }
  
  cluster_vars <- to_numeric(as.data.frame(cluster_vars))
  if (!is.null(disc_vars)) {
    disc_vars <- to_numeric(as.data.frame(disc_vars))
  } else {
    disc_vars <- cluster_vars
  }
  
  if (nrow(cluster_vars) != nrow(assignments)) {
    stop(
      "`cluster_vars` has ", nrow(cluster_vars),
      " rows but `assignments` has ", nrow(assignments),
      " rows. Align respondents before calling."
    )
  }
  if (nrow(disc_vars) != nrow(assignments)) {
    stop(
      "`disc_vars` has ", nrow(disc_vars),
      " rows but `assignments` has ", nrow(assignments),
      " rows. Align respondents before calling."
    )
  }
  
  # ---- 2. Adjusted Rand Index matrix ----
  ari_matrix <- matrix(NA_real_, n_solutions, n_solutions)
  rownames(ari_matrix) <- colnames(ari_matrix) <- sol_names
  
  for (i in seq_len(n_solutions)) {
    for (j in i:n_solutions) {
      ari <- mclust::adjustedRandIndex(assignments[[i]], assignments[[j]])
      ari_matrix[i, j] <- ari_matrix[j, i] <- ari
    }
  }
  ari_df <- as.data.frame(ari_matrix)
  
  # ---- 3. Size & balance ----
  size_summary <- purrr::map_dfr(sol_names, function(nm) {
    tab   <- table(assignments[[nm]])
    props <- prop.table(tab)
    entropy <- -sum(props * log(props, base = 2), na.rm = TRUE)
    
    data.frame(
      Solution        = nm,
      N_Clusters      = length(tab),
      Min_N           = as.numeric(min(tab)),
      Max_N           = as.numeric(max(tab)),
      Min_Pct         = round(min(props) * 100, 1),
      Max_Pct         = round(max(props) * 100, 1),
      Entropy         = round(entropy, 3),
      Flag_Small      = min(props) * 100 < min_size_pct,
      Flag_Unbalanced = max(props) > 0.60,
      stringsAsFactors = FALSE
    )
  })
  
  # ---- 4. Differentiation (ANOVA) ----
  diff_metrics <- purrr::map_dfr(sol_names, function(nm) {
    cl <- as.factor(assignments[[nm]])
    
    f_stats <- vapply(cluster_vars, function(v) {
      tryCatch({
        a <- stats::aov(v ~ cl)
        s <- summary(a)[[1]]
        as.numeric(s[["F value"]][1])
      }, error = function(e) NA_real_)
    }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
    
    p_vals <- vapply(cluster_vars, function(v) {
      tryCatch({
        a <- stats::aov(v ~ cl)
        s <- summary(a)[[1]]
        as.numeric(s[["Pr(>F)"]][1])
      }, error = function(e) NA_real_)
    }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
    
    data.frame(
      Solution           = nm,
      Mean_F             = round(mean(f_stats, na.rm = TRUE), 2),
      Median_F           = round(stats::median(f_stats, na.rm = TRUE), 2),
      Prop_Significant   = round(mean(p_vals < 0.05, na.rm = TRUE), 3),
      N_Significant_Vars = sum(p_vals < 0.05, na.rm = TRUE),
      stringsAsFactors   = FALSE
    )
  })
  
  # ---- 5. LDA ----
  hit_resub_df <- data.frame(
    Solution = sol_names,
    Hit_Rate_Resubstitution = NA_real_,
    stringsAsFactors = FALSE
  )
  hit_cv_df <- data.frame(
    Solution = sol_names,
    Hit_Rate_CV = NA_real_,
    stringsAsFactors = FALSE
  )
  confusion_list       <- list()
  variable_impact_list <- list()
  
  if (isTRUE(do_lda)) {
    for (i in seq_along(sol_names)) {
      nm <- sol_names[i]
      cl <- as.factor(assignments[[nm]])
      
      if (length(unique(stats::na.omit(cl))) < 2) {
        confusion_list[[nm]] <- "Fewer than 2 clusters"
        next
      }
      
      dat <- data.frame(cl = cl, disc_vars, stringsAsFactors = FALSE)
      dat <- stats::na.omit(dat)
      
      if (nrow(dat) < 2 || length(unique(dat$cl)) < 2) {
        confusion_list[[nm]] <- "Insufficient complete cases for LDA"
        next
      }
      
      lda_fit <- tryCatch(
        MASS::lda(cl ~ ., data = dat, CV = FALSE),
        error = function(e) NULL
      )
      lda_cv <- tryCatch(
        MASS::lda(cl ~ ., data = dat, CV = TRUE),
        error = function(e) NULL
      )
      
      if (is.null(lda_fit)) {
        confusion_list[[nm]] <- "LDA failed"
        next
      }
      
      pred_resub <- predict(lda_fit)$class
      hit_resub_df$Hit_Rate_Resubstitution[i] <- round(mean(pred_resub == dat$cl), 3)
      hit_cv_df$Hit_Rate_CV[i] <- if (!is.null(lda_cv)) {
        round(mean(lda_cv$class == dat$cl), 3)
      } else {
        NA_real_
      }
      confusion_list[[nm]] <- table(Actual = dat$cl, Predicted = pred_resub)
      
      impact <- vapply(disc_vars, function(v) {
        tryCatch({
          ss <- summary(stats::aov(v ~ cl))[[1]]
          as.numeric(ss["cl", "Sum Sq"] / sum(ss[["Sum Sq"]], na.rm = TRUE))
        }, error = function(e) NA_real_)
      }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
      
      variable_impact_list[[nm]] <- data.frame(
        Variable    = names(disc_vars),
        Eta_Squared = as.numeric(impact),
        stringsAsFactors = FALSE
      ) |>
        dplyr::arrange(dplyr::desc(.data$Eta_Squared))
    }
  }
  
  # ---- 6. Overall evaluation score ----
  eval_df <- size_summary |>
    dplyr::left_join(diff_metrics, by = "Solution") |>
    dplyr::left_join(hit_cv_df, by = "Solution") |>
    dplyr::left_join(hit_resub_df, by = "Solution")
  
  eval_df <- eval_df |>
    dplyr::mutate(
      Score_CV   = ifelse(is.na(.data$Hit_Rate_CV), 0, .data$Hit_Rate_CV * 100),
      Score_Diff = ifelse(is.na(.data$Prop_Significant), 0, .data$Prop_Significant * 100),
      Score_Size = 100 - pmax(0, min_size_pct - .data$Min_Pct) * 4 -
        pmax(0, .data$Max_Pct - 55) * 1.5,
      Score_Size = pmax(0, pmin(100, .data$Score_Size)),
      Score_F    = pmin(100, ifelse(is.na(.data$Mean_F), 0, .data$Mean_F * 3)),
      Overall_Score = round(
        0.40 * .data$Score_CV +
          0.30 * .data$Score_Diff +
          0.20 * .data$Score_Size +
          0.10 * .data$Score_F,
        1
      ),
      Flag_NoGo_Size = .data$Flag_Small | .data$Flag_Unbalanced,
      Flag_NoGo_Sep  = !is.na(.data$Hit_Rate_CV) & .data$Hit_Rate_CV < 0.60,
      Flag_NoGo_Diff = !is.na(.data$Prop_Significant) & .data$Prop_Significant < 0.40,
      Overall_Flag = dplyr::case_when(
        .data$Flag_NoGo_Size | .data$Flag_NoGo_Sep | .data$Flag_NoGo_Diff ~ "No-Go",
        .data$Overall_Score >= 75 ~ "Strong",
        .data$Overall_Score >= 60 ~ "Acceptable",
        TRUE ~ "Weak"
      )
    ) |>
    dplyr::select(dplyr::any_of(c(
      "Solution", "N_Clusters", "Min_Pct", "Max_Pct", "Entropy",
      "Mean_F", "Prop_Significant",
      "Hit_Rate_CV", "Hit_Rate_Resubstitution",
      "Overall_Score", "Overall_Flag"
    )))
  
  # ---- 7. Migration tables ----
  migration_tables <- list()
  if (isTRUE(do_migration)) {
    method_prefixes <- unique(gsub("[_]?[0-9]+$", "", sol_names))
    for (pref in method_prefixes) {
      cols <- grep(paste0("^", pref), sol_names, value = TRUE)
      if (length(cols) >= 2) {
        cols <- cols[order(as.numeric(gsub("[^0-9]", "", cols)))]
        tabs <- list()
        for (i in 2:length(cols)) {
          tabs[[paste0(cols[i - 1], " to ", cols[i])]] <-
            table(assignments[[cols[i - 1]]], assignments[[cols[i]]])
        }
        migration_tables[[pref]] <- tabs
      }
    }
  }
  
  # ---- return ----
  list(
    Solution_Evaluation   = eval_df,
    ARI_Matrix            = ari_df,
    Size_and_Balance      = size_summary,
    Differentiation_Metric = diff_metrics,
    Discriminant_Impact   = variable_impact_list,
    Confusion_Matrices    = confusion_list,
    Migration_Tables      = if (isTRUE(do_migration)) migration_tables else "Not requested"
  )
}
