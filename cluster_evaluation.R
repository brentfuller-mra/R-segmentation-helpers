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
  
  # ---- 1. Safe Alignment & Assignment Extraction ----
  # If passed directly as a raw output list object from generation engine
  if (is.list(assignments) && !is.data.frame(assignments) && !is.null(assignments$Cluster_Assignments)) {
    assignments <- assignments$Cluster_Assignments
  }
  
  if (is.matrix(assignments)) assignments <- as.data.frame(assignments)
  if (is.vector(assignments)) assignments <- data.frame(Solution = assignments)
  
  # Strip away unique ID column if present to isolate solutions
  if (names(assignments)[1] %in% c("id", "resp_id", "Case")) {
    assignments <- assignments[, -1, drop = FALSE]
  }
  
  sol_names   <- names(assignments)
  n_solutions <- ncol(assignments)
  
  if (n_solutions < 1) stop("Please supply at least one cluster solution column.")
  
  # Ensure parsing structures are clean numeric data frames
  to_numeric <- function(df) {
    as.data.frame(lapply(df, function(x) {
      if (is.factor(x) || is.ordered(x)) as.numeric(x) else as.numeric(x)
    }))
  }
  
  cluster_vars <- to_numeric(as.data.frame(cluster_vars))
  if (!is.null(disc_vars)) {
    disc_vars <- to_numeric(as.data.frame(disc_vars))
  } else {
    disc_vars <- cluster_vars # Fallback to cluster inputs if no demographics provided
  }
  
  # ---- 2. Adjusted Rand Index Stability Matrix ----
  ari_matrix <- matrix(NA, n_solutions, n_solutions)
  rownames(ari_matrix) <- colnames(ari_matrix) <- sol_names
  
  for (i in 1:n_solutions) {
    for (j in i:n_solutions) {
      ari <- mclust::adjustedRandIndex(assignments[[i]], assignments[[j]])
      ari_matrix[i, j] <- ari_matrix[j, i] <- ari
    }
  }
  ari_df <- as.data.frame(ari_matrix)
  
  # ---- 3. Size & Balance Topology Metrics ----
  size_summary <- purrr::map_dfr(sol_names, function(nm) {
    tab   <- table(assignments[[nm]])
    props <- prop.table(tab)
    entropy <- -sum(props * log(props, base = 2), na.rm = TRUE)
    
    data.frame(
      Solution        = nm,
      N_Clusters      = length(tab),
      Min_N           = min(tab),
      Max_N           = max(tab),
      Min_Pct         = round(min(props) * 100, 1),
      Max_Pct         = round(max(props) * 100, 1),
      Entropy         = round(entropy, 3),
      Flag_Small      = min(props) * 100 < min_size_pct,
      Flag_Unbalanced = max(props) > 0.60
    )
  })
  
  # ---- 4. Differentiation (ANOVA Space Mappings) ----
  diff_metrics <- purrr::map_dfr(sol_names, function(nm) {
    cl <- as.factor(assignments[[nm]])
    
    f_stats <- sapply(cluster_vars, function(v) {
      tryCatch(summary(aov(v ~ cl))[[1]][["F value"]][1], error = function(e) NA)
    })
    p_vals <- sapply(cluster_vars, function(v) {
      tryCatch(summary(aov(v ~ cl))[[1]][["Pr(>F)"]][1], error = function(e) NA)
    })
    
    data.frame(
      Solution           = nm,
      Mean_F             = round(mean(f_stats, na.rm = TRUE), 2),
      Median_F           = round(median(f_stats, na.rm = TRUE), 2),
      Prop_Significant   = round(mean(p_vals < 0.05, na.rm = TRUE), 3),
      N_Significant_Vars = sum(p_vals < 0.05, na.rm = TRUE)
    )
  })
  
  # ---- 5. Linear Discriminant Classification Analysis ----
  hit_resub_df <- data.frame(Solution = sol_names, Hit_Rate_Resubstitution = NA_real_)
  hit_cv_df    <- data.frame(Solution = sol_names, Hit_Rate_CV = NA_real_)
  confusion_list       <- list()
  variable_impact_list <- list()
  
  if (do_lda) {
    for (i in seq_along(sol_names)) {
      nm <- sol_names[i]
      cl <- as.factor(assignments[[nm]])
      
      if (length(unique(cl)) < 2) {
        confusion_list[[nm]] <- "Fewer than 2 clusters"
        next
      }
      
      dat <- data.frame(cl = cl, disc_vars)
      dat <- na.omit(dat)
      
      lda_fit <- tryCatch(MASS::lda(cl ~ ., data = dat, CV = FALSE), error = function(e) NULL)
      lda_cv  <- tryCatch(MASS::lda(cl ~ ., data = dat, CV = TRUE),  error = function(e) NULL)
      
      if (is.null(lda_fit)) {
        confusion_list[[nm]] <- "LDA failed"
        next
      }
      
      pred_resub <- predict(lda_fit)$class
      hit_resub_df$Hit_Rate_Resubstitution[i] <- round(mean(pred_resub == dat$cl), 3)
      hit_cv_df$Hit_Rate_CV[i] <- if (!is.null(lda_cv)) round(mean(lda_cv$class == dat$cl), 3) else NA
      confusion_list[[nm]] <- table(Actual = dat$cl, Predicted = pred_resub)
      
      # Variable impact rankings (Eta-Squared calculations)
      impact <- sapply(disc_vars, function(v) {
        tryCatch({
          ss <- summary(aov(v ~ cl))[[1]]
          round(ss["cl", "Sum Sq"] / sum(ss$"Sum Sq", na.rm = TRUE), 3)
        }, error = function(e) NA)
      })
      
      variable_impact_list[[nm]] <- data.frame(
        Variable    = names(disc_vars),
        Eta_Squared = as.numeric(impact)
      ) |> dplyr::arrange(dplyr::desc(Eta_Squared))
    }
  }
  
  # ---- 6. Synthesis Master Evaluation Score Modeling ----
  eval_df <- size_summary |>
    dplyr::left_join(diff_metrics, by = "Solution") |>
    dplyr::left_join(hit_cv_df, by = "Solution") |>
    dplyr::left_join(hit_resub_df, by = "Solution")
  
  eval_df <- eval_df |>
    dplyr::mutate(
      Score_CV   = ifelse(is.na(Hit_Rate_CV), 0, Hit_Rate_CV * 100),
      Score_Diff = ifelse(is.na(Prop_Significant), 0, Prop_Significant * 100),
      Score_Size = 100 - pmax(0, min_size_pct - Min_Pct) * 4 - pmax(0, Max_Pct - 55) * 1.5,
      Score_Size = pmax(0, pmin(100, Score_Size)),
      Score_F    = pmin(100, ifelse(is.na(Mean_F), 0, Mean_F * 3)),
      Overall_Score = round(0.40 * Score_CV + 0.30 * Score_Diff + 0.20 * Score_Size + 0.10 * Score_F, 1),
      Flag_NoGo_Size = Flag_Small | Flag_Unbalanced,
      Flag_NoGo_Sep  = !is.na(Hit_Rate_CV) & Hit_Rate_CV < 0.60,
      Flag_NoGo_Diff = !is.na(Prop_Significant) & Prop_Significant < 0.40,
      Overall_Flag = dplyr::case_when(
        Flag_NoGo_Size | Flag_NoGo_Sep | Flag_NoGo_Diff ~ "No-Go",
        Overall_Score >= 75 ~ "Strong",
        Overall_Score >= 60 ~ "Acceptable",
        TRUE ~ "Weak")) |> 
    dplyr::select(
      dplyr::any_of(c("Solution", "N_Clusters", "Min_Pct", "Max_Pct", "Entropy","Mean_F", "Prop_Significant", "Hit_Rate_CV", "Hit_Rate_Resubstitution","Overall_Score", "Overall_Flag"))
    )
  
  #---- 7. Cluster Cut Migration Pathways ----
  migration_tables <- list()
  if (do_migration) {
    method_prefixes <- unique(gsub("[_]?[0-9]+$", "", sol_names))
    for (pref in method_prefixes) {
      cols <- grep(paste0("^", pref), sol_names, value = TRUE)
      if (length(cols) >= 2) {
        # Sort sequentially by extracted numeric cluster numbers
        #cols <- cols[order(as.numeric(gsub("\\D", "", cols)))]
        cols <- cols[order(as.numeric(gsub("[^0-9]", "", cols)))]
        tabs <- list()
        for (i in 2:length(cols)) {
          tabs[[paste0(cols[i-1], " to ", cols[i])]] <- table(assignments[[cols[i-1]]], assignments[[cols[i]]])
        }
        migration_tables[[pref]] <- tabs
      }
    }
  }
  #---- Return Unified Structure ----
  list(Solution_Evaluation=eval_df,
       ARI_Matrix=ari_df,
       Size_and_Balance=size_summary,
       Differentiation_Metric=diff_metrics,
       Discriminant_Impact=variable_impact_list,
       Confusion_Matrices=confusion_list,
       Migration_Tables=if (do_migration) migration_tables else "Not requested"
  )
  }