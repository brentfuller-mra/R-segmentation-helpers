# Test segment generation and evaluation workflow functions
# Not showcasing an actual workflow, just showing how the functions work

seg_packages <- c(
  "cluster", "psych", "NMF", "archetypes", "mclust", "flexclust", 
  "poLCA", "poLCAParallel", "kohonen", "randomForest", "uwot", "dbscan", "MASS",
  "haven", "purrr", "dplyr", "tidyverse" # Tidyverse loaded DEAD LAST to own master namespace
)
invisible(lapply(seg_packages, library, character.only = TRUE))

source("https://raw.githubusercontent.com/brentfuller-mra/R-segmentation-helpers/main/cluster_generation.R")
source("https://raw.githubusercontent.com/brentfuller-mra/R-segmentation-helpers/main/cluster_evaluation.R")
survey_data <- haven::read_sav("https://raw.githubusercontent.com/brentfuller-mra/R-segmentation-helpers/main/data_example.sav")

# input_data
# First column = respondent ID string, remaining columns are inputs
input_data <- survey_data |> 
  dplyr::select(ID, quality, value_for_money, innovation, trustworthiness, 
                sustainability, ease_of_use, customer_service, modern)

# Discrete categorical copy specifically isolated for testing LCA
categorical_data <- survey_data |>
  dplyr::select(ID, quality, value_for_money, innovation) |>
  dplyr::mutate(across(-ID, ~ as.integer(cut(.x, breaks = 3, labels = c(1, 2, 3)))))


km_results <- cluster_km(data = input_data, num_solutions = 3:4, method = "kmeans_seeded", algorithm = "Lloyd", id_col = TRUE)
flex_results <- cluster_flex(data = input_data, num_solutions = 3:4, family = "kmeans", initcent = "kmeanspp", id_col = TRUE)
clara_results <- cluster_clara(data = input_data, num_solutions = 3:4, samples = 5L, sampsize = 100L, id_col = TRUE)
  primary_vars   <- input_data |> dplyr::select(ID,quality, value_for_money, innovation, trustworthiness)
  secondary_vars <- input_data |> dplyr::select(ID,sustainability, ease_of_use, customer_service, modern)
cancor_results <- cluster_cancor(primary = primary_vars, secondary = secondary_vars, num_solutions = 3:4, seed = 123L, id_col = TRUE)
arch_results <- cluster_arch(data = input_data, num_solutions = 3:4, id_col = TRUE)
fa_results <- cluster_fa(data = input_data, cor_type= "cor", num_solutions = 3:4, fm = "pa", scoring = "max", id_col = TRUE)
nmf_results <- cluster_nmf(data = input_data, num_solutions = 3:4, method = "lee", nrun = 5L, seed = 123L, id_col = TRUE)
mod_results <- cluster_mod(data = input_data, num_solutions = 3:4, model_names = "EII", id_col = TRUE)
twostep_results <- cluster_twostep(data = input_data, num_solutions = 3:4, id_col = TRUE)
lca_results <- cluster_lca(data = categorical_data, num_solutions = 3:4, nrep = 2L, maxiter = 500L, id_col = TRUE)
lca2_results <- cluster_lca2(data = categorical_data, num_solutions = 3:4, nrep = 2L, maxiter = 500L, n_thread = 2L, id_col = TRUE)
som_results <- cluster_som(data = input_data, num_solutions = 3:4, topo = "hexagonal", rlen = 50L, id_col = TRUE)
rf_results <- cluster_rf(data = input_data, num_solutions = 3:4, ntree = 200L, seed = 123L, id_col = TRUE)
umap_hdb_results <- cluster_umap_hdbscan(data = input_data, num_solutions = 3:4, n_neighbors = 10L, id_col = TRUE)
autoenc_results <- cluster_autoencoder(data = input_data, num_solutions = 3:4, layers = c(6, 3), seed = 123L, id_col = TRUE)

all_assignments <- km_results$Cluster_Assignments |>
  dplyr::left_join(flex_results$Cluster_Assignments,      by = "ID") |>
  dplyr::left_join(clara_results$Cluster_Assignments,     by = "ID") |> 
  dplyr::left_join(cancor_results$Cluster_Assignments,    by = "ID") |>
  dplyr::left_join(arch_results$Cluster_Assignments,      by = "ID") |>
  dplyr::left_join(fa_results$Cluster_Assignments,        by = "ID") |>
  dplyr::left_join(nmf_results$Cluster_Assignments,       by = "ID") |>
  dplyr::left_join(mod_results$Cluster_Assignments,       by = "ID") |>
  dplyr::left_join(twostep_results$Cluster_Assignments,   by = "ID") |>
  dplyr::left_join(lca_results$Cluster_Assignments,       by = "ID") |>
  dplyr::left_join(lca2_results$Cluster_Assignments,      by = "ID") |>
  dplyr::left_join(som_results$Cluster_Assignments,       by = "ID") |>
  dplyr::left_join(rf_results$Cluster_Assignments,        by = "ID") |>
  dplyr::left_join(umap_hdb_results$Cluster_Assignments,  by = "ID") |>
  dplyr::left_join(autoenc_results$Cluster_Assignments,   by = "ID")

dim(all_assignments)

ensemble_results <- cluster_ensemble(
  data          = all_assignments, 
  num_solutions = 3:4, 
  distance      = "jaccard", 
  method        = "kmedians", 
  id_col        = TRUE
)

# Append final ensemble segs
all_assignments_out <- all_assignments |>
  dplyr::left_join(ensemble_results$Cluster_Assignments, by = "ID")

eval_output <- evaluate_segmentations(
  assignments  = all_assignments_out,
  cluster_vars = input_data |> dplyr::select(-ID),
  disc_vars    = NULL,
  min_size_pct = 5.0,
  do_lda       = T,
  do_migration = T
)

eval_output$Size_and_Balance
eval_output$Solution_Evaluation
eval_output$ARI_Matrix[1:6, 1:6]
eval_output$ARI_Matrix
eval_output$Differentiation_Metric
eval_output$Migration_Tables
eval_output$Confusion_Matrices
eval_output$Discriminant_Impact
