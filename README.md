# R-segmentation-helpers
A collection of custom R wrappers and functions designed to aid in market research segmentation generation, exploration, and evaluation workflows. 

## 🛠️ Cluster Generation Engine (`cluster_generation.R`)

All functions accept an optional ID column (e.g., `resp_id`) as the first column, cycle through a vector of cluster solutions (e.g., `3:8`), and return a unified list structure.

### Available Algorithms & Functions

| Function | Method | Used For... | Output Example |
| :--- | :--- | :--- | :--- |
| `cluster_km()` | K-Means & Hierarchical | Fast baseline solutions, profile-shape clustering via correlation distance, or Ward's linkage. | `kmLloyd_k`, `hWard_k`, `hCorr_k` |
| `cluster_flex()` | K-Centroids (`flexclust`) | K-medians for outlier control, or Jaccard/Angular distances for binary data grids. | `flexKmeans_k` |
| `cluster_clara()` | Large Application Medoids | Scaling up segment definitions to exceptionally massive customer CRM databases ($N > 100,000$). | `claraEucl_k` |
| `cluster_cancor()` | Canonical Correlation K-Means | Linking primary segment drivers (e.g., needs) directly with secondary metrics (e.g., behaviors). | `kmcc_k` |
| `cluster_arch()` | Archetypal Analysis | Finding "extreme" prototypical respondents ("standard bearer") and treating others as blends of these archetypes. | `arch_k` |
| `cluster_fa()` | Factor Analysis Clustering | Condensing massive attribute batteries into dimensions via `psych::fa` before grouping. | `factsCorPa_k` |
| `cluster_nmf()` | Non-Negative Matrix Factorization | MaxDiff, pick-any binary grids, or strictly positive attitudinal datasets. | `nmfLee_k` |
| `cluster_mod()` | Model-Based Mixture Models | Gaussian Mixture Models (`mclust`) when segments vary drastically in size, shape, or orientation. | `modAll_k` |
| `cluster_twostep()` | SPSS Two-Step Emulation | Replicating the dual data-compression and log-likelihood cluster generation framework native to SPSS. | `twostep_k` |
| `cluster_lca()` | Latent Class Analysis | Categorical survey questions, pick-any matrices, or ordinal scales where categorical distributions determine groups. | `lca_k` |
| `cluster_lca2()` | High-Performance LCA | C++ accelerated multi-threaded Latent Class Analysis (via `poLCAParallel`). Dramatically cuts computation runtimes on multi-core CPUs. | `lcaParallel_k` |
| `cluster_som()` | Self-Organizing Maps | Uncovering multi-dimensional topologies where cluster positioning on a grid displays structural likeness. | `somHex_k` |
| `cluster_rf()` | Unsupervised Random Forest | Finding organic groupings via tree-ensemble co-occurrence proximities. | `rfWard_k` |
| `cluster_umap_hdbscan()` | UMAP + HDBSCAN | Modern manifold learning combined with density clustering. Identifies organic cluster boundaries and automatically separates noisy outliers (Cluster 0). | `umaphdb_k` |
| `cluster_autoencoder()` | Deep Representation Space | Utilizing non-linear neural compression network layers to build a deep feature bottleneck space before segment isolation. | `deepkm_k` |

### Output Structure
Every function returns a `list` containing:
1. **`Cluster_Assignments`**: Data frame with your original ID column and hard cluster flags for every solution count.
2. **`Cluster_Sizes`**: Tidy long-form data frame showing counts ($N$) per cluster variant (perfect for quick `ggplot` checks).
3. **Algorithm Specific Entries**: Extended outputs like `Centroids`, `Soft_Probabilities`, `Alphas`, etc.

---

## 🛠️ Cluster Evaluation Engine (`cluster_evaluation.R`)
The `evaluate_segmentations()` function processes multiple cluster assignment outputs, scoring them on cross-validation classification (`MASS::lda`), sample balancing, and stability.

### Key Evaluation Dimensions Evaluated
1. **Adjusted Rand Index (ARI)**: Verifies crossover consistency to see how solutions overlap.
2. **Information Entropy**: Tracks allocation variance to flag solutions dominated by a massive single group (>60%) or split by micro-segments (<5%).
3. **ANOVA Driver Metrics**: Compares input factor separation via aggregate Mean F-statistics.
4. **Predictive Classification (LDA Cross-Validation)**: Calculates hit rates.

## 🛠️ Cluster Exploration Engine (`cluster_exploration.R`)
a.k.a. profiling, naming, description, & viz
T.K.

## ✨ Usage
```r
source("https://raw.githubusercontent.com/brentfuller-mra/R-segmentation-helpers/main/cluster_generation.R")

km_results <- cluster_km(
  data          = input_data,  # First column = respondent ID string, remaining columns are inputs
  num_solutions = 3:6,
  method        = "kmeans_seeded",
  algorithm     = "Lloyd",
  id_col        = TRUE
)

nmf_results <- cluster_nmf(
  data          = input_data,  # First column = respondent ID string, remaining columns are inputs
  num_solutions = 3:6,
  method        = "lee",
  nrun          = 30L,
  seed          = 123L,
  scoring       = "max",
  check_nonnegative = TRUE,
  id_col = TRUE
)

km_results$Cluster_Assignments
nmf_results$Cluster_Assignments

all_assngments <- km_results$Cluster_Assignments |> 
  left_join(nmf_results$Cluster_Assignments)

eval_output <- evaluate_segmentations(assignments=all_assngments,
                                   cluster_vars=input_data[,-1],
                                   disc_vars = NULL,
                                   min_size_pct = 5.0,
                                   do_lda = TRUE,
                                   do_migration = TRUE)
eval_output$ARI_Matrix
eval_output$Solution_Evaluation
eval_output$Size_and_Balance
eval_output$Differentiation_Metric
eval_output$Migration_Tables
eval_output$Confusion_Matrices
eval_output$Discriminant_Impact

```
