library(glmnet)
library(caret)
library(pROC)
library(randomForest)
library(nnet)
library(ggplot2)

# ==============================================================================
# 0. SETUP
# ==============================================================================
BASE_PATH <- "C:/UPC/Thesis/Rcode/"

# Load all data sources once
demo       <- read.csv(paste0(BASE_PATH, "demographics_cleaned_filtered_259.csv"))
afc_matrix <- read.csv(paste0(BASE_PATH, "AFC_vectors.csv"))
afc_pca    <- read.csv(paste0(BASE_PATH, "AFC_PCA_scores_120.csv"))
afc_score  <- read.csv(paste0(BASE_PATH, "afc_scores.csv"))
irr_df     <- read.csv(paste0(BASE_PATH, "irrli_259.csv"))
features_simple <- read.csv(paste0(BASE_PATH, "old_features/features_simple.csv"))
graph_features  <- read.csv(paste0(BASE_PATH, "graph_features.csv"))

# Standardize ID columns
for (df_name in c("afc_matrix", "afc_pca", "irr_df")) {
  df_obj <- get(df_name)
  if ("subject_id" %in% names(df_obj))
    names(df_obj)[names(df_obj) == "subject_id"] <- "ID"
  assign(df_name, df_obj)
}
if ("subject_id" %in% names(afc_score))
  names(afc_score)[names(afc_score) == "subject_id"] <- "ID"
afc_score_data_cols <- setdiff(names(afc_score), "ID")
names(afc_score)[names(afc_score) == afc_score_data_cols] <- "afc_score"

names(features_simple)[names(features_simple) == "subject_id"] <- "ID"
features_simple$ID <- sub("^sub-", "", features_simple$ID)
names(graph_features)[names(graph_features) == "subject_id"] <- "ID"
graph_features$ID  <- sub("^sub-", "", graph_features$ID)

# ==============================================================================
# 1. BUILD df_full - merge ALL sources at once
# ==============================================================================
pc_cols_all <- paste0("PC", 1:120)
pc_cols     <- paste0("PC", 1:10)

df_full <- merge(
  demo[, c("ID", "HEMISPHERIC", "HEMISPHERIC_2",
           "Laterality_Quotient", "LI_GEN", "cohort", "Edad", "Sexo")],
  afc_matrix, by = "ID")
df_full <- merge(df_full, afc_pca[, c("ID", pc_cols_all)], by = "ID")
df_full <- merge(df_full, afc_score[, c("ID", "afc_score")], by = "ID")
df_full <- merge(df_full, irr_df[, c("ID", "IrrLI")], by = "ID")
df_full <- merge(df_full,
                 features_simple[, !names(features_simple) %in% "LI_GEN"],
                 by = "ID")
df_full <- merge(df_full, graph_features, by = "ID")
df_full <- df_full[!is.na(df_full$HEMISPHERIC) &
                     !is.na(df_full$Laterality_Quotient), ]
df_full$sex_bin <- ifelse(df_full$Sexo == "Mujer", 1, 0)

cat("Final n:", nrow(df_full), "\n")
print(table(df_full$HEMISPHERIC))
print(table(df_full$HEMISPHERIC_2))

# ==============================================================================
# 2. FEATURE DEFINITIONS
# ==============================================================================
na_counts  <- sapply(pc_cols_all, function(col) sum(is.na(df_full[[col]])))
max_safe_k <- ifelse(any(na_counts > 0), which(na_counts > 0)[1] - 1, 120)
cat("Max usable PCs:", max_safe_k, "\n")

y        <- ifelse(df_full$HEMISPHERIC == "Typical", 1, 0)
y_3class <- factor(df_full$HEMISPHERIC_2)

# HEMISPHERIC_3: three-class labels using -20/+20 LI_GEN threshold
df_full$HEMISPHERIC_3 <- cut(
  df_full$LI_GEN,
  breaks = c(-Inf, -20, 20, Inf),
  labels = c("Right", "Bilateral", "Left")
)

# HEMISPHERIC_4: binary labels derived from HEMISPHERIC_3
df_full$HEMISPHERIC_4 <- ifelse(df_full$HEMISPHERIC_3 == "Left",
                                "Typical", "Atypical")
y_h3 <- ifelse(df_full$HEMISPHERIC_3 == "Left", 1, 0)

cat("\nHEMISPHERIC_3 distribution:\n")
print(table(df_full$HEMISPHERIC_3))
cat("\nHEMISPHERIC_4 distribution:\n")
print(table(df_full$HEMISPHERIC_4))
cat("\nCrosstab original vs new binary label:\n")
print(table(Original = df_full$HEMISPHERIC, New = df_full$HEMISPHERIC_4))
afc_cols <- names(afc_matrix)[names(afc_matrix) != "ID"]
ezm2_idx <- which(df_full$cohort == 2)
ezm3_idx <- which(df_full$cohort == 3)

# Strategy 1 feature matrices
feat_pca_lq_age     <- c(pc_cols, "Laterality_Quotient", "Edad")
feat_pca_lq_age_irr <- c(pc_cols, "Laterality_Quotient", "Edad", "IrrLI")

X_afc            <- as.matrix(df_full[, afc_cols])
X_pca            <- as.matrix(df_full[, pc_cols])
X_pca_lq         <- as.matrix(df_full[, c(pc_cols, "Laterality_Quotient")])
X_pca_lq_age     <- as.matrix(df_full[, feat_pca_lq_age])
X_pca_lq_age_irr <- as.matrix(df_full[, feat_pca_lq_age_irr])

# Strategy 2 feature matrices
roi_idx_simple <- c(1:46, 48:101)
nodeai_cols    <- paste0("nodeAI_R",    roi_idx_simple)
homotopic_cols <- paste0("homotopic_R", roi_idx_simple)
graph_cols     <- c(paste0("asym_deg_R",   1:100),
                    paste0("asym_clust_R", 1:100),
                    paste0("asym_btwn_R",  1:100),
                    paste0("asym_leff_R",  1:100))

stopifnot(all(nodeai_cols    %in% names(df_full)))
stopifnot(all(homotopic_cols %in% names(df_full)))
stopifnot(all(graph_cols     %in% names(df_full)))

feat_nodeai_lq     <- c(nodeai_cols, "Laterality_Quotient", "Edad")
feat_conn_lq       <- c(nodeai_cols, homotopic_cols,
                        "Laterality_Quotient", "Edad")
feat_graph_lq      <- c(graph_cols, "Laterality_Quotient", "Edad")
feat_all_interp_lq <- c(nodeai_cols, homotopic_cols, graph_cols,
                        "Laterality_Quotient", "Edad")

X_nodeai        <- as.matrix(df_full[, nodeai_cols])
X_nodeai_lq     <- as.matrix(df_full[, feat_nodeai_lq])
X_conn_lq       <- as.matrix(df_full[, feat_conn_lq])
X_graph_lq      <- as.matrix(df_full[, feat_graph_lq])
X_all_interp_lq <- as.matrix(df_full[, feat_all_interp_lq])

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================
run_ridge_cv <- function(X_mat, y_vec, folds) {
  preds <- numeric(length(y_vec))
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_along(y_vec), train_idx)
    tr_means  <- colMeans(X_mat[train_idx, , drop = FALSE])
    tr_sds    <- apply(X_mat[train_idx, , drop = FALSE], 2, sd)
    tr_sds[tr_sds == 0] <- 1
    X_tr <- scale(X_mat[train_idx, , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    X_te <- scale(X_mat[test_idx,  , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    fit  <- cv.glmnet(X_tr, y_vec[train_idx],
                      alpha = 0, family = "binomial", nfolds = 5)
    preds[test_idx] <- predict(fit, X_te, s = "lambda.min",
                               type = "response")[, 1]
  }
  preds
}

library(glmnet)
library(caret)
library(pROC)
library(randomForest)
library(nnet)
library(ggplot2)

# ==============================================================================
# 0. SETUP
# ==============================================================================
BASE_PATH <- "C:/UPC/Thesis/Rcode/"

# Load all data sources once
demo       <- read.csv(paste0(BASE_PATH, "demographics_cleaned_filtered_259.csv"))
afc_matrix <- read.csv(paste0(BASE_PATH, "AFC_vectors.csv"))
afc_pca    <- read.csv(paste0(BASE_PATH, "AFC_PCA_scores_120.csv"))
afc_score  <- read.csv(paste0(BASE_PATH, "afc_scores.csv"))
irr_df     <- read.csv(paste0(BASE_PATH, "irrli_259.csv"))
features_simple <- read.csv(paste0(BASE_PATH, "old_features/features_simple.csv"))
graph_features  <- read.csv(paste0(BASE_PATH, "graph_features.csv"))

# Standardize ID columns
for (df_name in c("afc_matrix", "afc_pca", "irr_df")) {
  df_obj <- get(df_name)
  if ("subject_id" %in% names(df_obj))
    names(df_obj)[names(df_obj) == "subject_id"] <- "ID"
  assign(df_name, df_obj)
}
if ("subject_id" %in% names(afc_score))
  names(afc_score)[names(afc_score) == "subject_id"] <- "ID"
afc_score_data_cols <- setdiff(names(afc_score), "ID")
names(afc_score)[names(afc_score) == afc_score_data_cols] <- "afc_score"

names(features_simple)[names(features_simple) == "subject_id"] <- "ID"
features_simple$ID <- sub("^sub-", "", features_simple$ID)
names(graph_features)[names(graph_features) == "subject_id"] <- "ID"
graph_features$ID  <- sub("^sub-", "", graph_features$ID)

# ==============================================================================
# 1. BUILD df_full - merge ALL sources at once
# ==============================================================================
pc_cols_all <- paste0("PC", 1:120)
pc_cols     <- paste0("PC", 1:10)

df_full <- merge(
  demo[, c("ID", "HEMISPHERIC", "HEMISPHERIC_2",
           "Laterality_Quotient", "LI_GEN", "cohort", "Edad", "Sexo")],
  afc_matrix, by = "ID")
df_full <- merge(df_full, afc_pca[, c("ID", pc_cols_all)], by = "ID")
df_full <- merge(df_full, afc_score[, c("ID", "afc_score")], by = "ID")
df_full <- merge(df_full, irr_df[, c("ID", "IrrLI")], by = "ID")
df_full <- merge(df_full,
                 features_simple[, !names(features_simple) %in% "LI_GEN"],
                 by = "ID")
df_full <- merge(df_full, graph_features, by = "ID")
df_full <- df_full[!is.na(df_full$HEMISPHERIC) &
                     !is.na(df_full$Laterality_Quotient), ]
df_full$sex_bin <- ifelse(df_full$Sexo == "Mujer", 1, 0)

cat("Final n:", nrow(df_full), "\n")
print(table(df_full$HEMISPHERIC))
print(table(df_full$HEMISPHERIC_2))

# ==============================================================================
# 2. FEATURE DEFINITIONS
# ==============================================================================
na_counts  <- sapply(pc_cols_all, function(col) sum(is.na(df_full[[col]])))
max_safe_k <- ifelse(any(na_counts > 0), which(na_counts > 0)[1] - 1, 120)
cat("Max usable PCs:", max_safe_k, "\n")

y        <- ifelse(df_full$HEMISPHERIC == "Typical", 1, 0)
y_3class <- factor(df_full$HEMISPHERIC_2)
afc_cols <- names(afc_matrix)[names(afc_matrix) != "ID"]
ezm2_idx <- which(df_full$cohort == 2)
ezm3_idx <- which(df_full$cohort == 3)

# Strategy 1 feature matrices
feat_pca_lq_age     <- c(pc_cols, "Laterality_Quotient", "Edad")
feat_pca_lq_age_irr <- c(pc_cols, "Laterality_Quotient", "Edad", "IrrLI")

X_afc            <- as.matrix(df_full[, afc_cols])
X_pca            <- as.matrix(df_full[, pc_cols])
X_pca_lq         <- as.matrix(df_full[, c(pc_cols, "Laterality_Quotient")])
X_pca_lq_age     <- as.matrix(df_full[, feat_pca_lq_age])
X_pca_lq_age_irr <- as.matrix(df_full[, feat_pca_lq_age_irr])

# Strategy 2 feature matrices
roi_idx_simple <- c(1:46, 48:101)
nodeai_cols    <- paste0("nodeAI_R",    roi_idx_simple)
homotopic_cols <- paste0("homotopic_R", roi_idx_simple)
graph_cols     <- c(paste0("asym_deg_R",   1:100),
                    paste0("asym_clust_R", 1:100),
                    paste0("asym_btwn_R",  1:100),
                    paste0("asym_leff_R",  1:100))

stopifnot(all(nodeai_cols    %in% names(df_full)))
stopifnot(all(homotopic_cols %in% names(df_full)))
stopifnot(all(graph_cols     %in% names(df_full)))

feat_nodeai_lq     <- c(nodeai_cols, "Laterality_Quotient", "Edad")
feat_conn_lq       <- c(nodeai_cols, homotopic_cols,
                        "Laterality_Quotient", "Edad")
feat_graph_lq      <- c(graph_cols, "Laterality_Quotient", "Edad")
feat_all_interp_lq <- c(nodeai_cols, homotopic_cols, graph_cols,
                        "Laterality_Quotient", "Edad")

X_nodeai        <- as.matrix(df_full[, nodeai_cols])
X_nodeai_lq     <- as.matrix(df_full[, feat_nodeai_lq])
X_conn_lq       <- as.matrix(df_full[, feat_conn_lq])
X_graph_lq      <- as.matrix(df_full[, feat_graph_lq])
X_all_interp_lq <- as.matrix(df_full[, feat_all_interp_lq])

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================
run_ridge_cv <- function(X_mat, y_vec, folds) {
  preds <- numeric(length(y_vec))
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_along(y_vec), train_idx)
    tr_means  <- colMeans(X_mat[train_idx, , drop = FALSE])
    tr_sds    <- apply(X_mat[train_idx, , drop = FALSE], 2, sd)
    tr_sds[tr_sds == 0] <- 1
    X_tr <- scale(X_mat[train_idx, , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    X_te <- scale(X_mat[test_idx,  , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    fit  <- cv.glmnet(X_tr, y_vec[train_idx],
                      alpha = 0, family = "binomial", nfolds = 5)
    preds[test_idx] <- predict(fit, X_te, s = "lambda.min",
                               type = "response")[, 1]
  }
  preds
}
run_ridge_cv_varfilter <- function(X_mat, y_vec, folds, var_threshold = 0.01) {
  preds <- numeric(length(y_vec))
  kept_counts <- numeric(length(folds))
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_along(y_vec), train_idx)
    
    # Variance filter on training fold only
    train_vars <- apply(X_mat[train_idx, , drop = FALSE], 2, var)
    keep <- train_vars > var_threshold
    kept_counts[i] <- sum(keep)
    
    X_tr_f <- X_mat[train_idx, keep, drop = FALSE]
    X_te_f <- X_mat[test_idx,  keep, drop = FALSE]
    
    tr_means <- colMeans(X_tr_f)
    tr_sds   <- apply(X_tr_f, 2, sd)
    tr_sds[tr_sds == 0] <- 1
    X_tr <- scale(X_tr_f, center = tr_means, scale = tr_sds)
    X_te <- scale(X_te_f, center = tr_means, scale = tr_sds)
    
    fit  <- cv.glmnet(X_tr, y_vec[train_idx],
                      alpha = 0, family = "binomial", nfolds = 5)
    preds[test_idx] <- predict(fit, X_te, s = "lambda.min",
                               type = "response")[, 1]
  }
  cat(sprintf("Features kept per fold: min=%d, max=%d, mean=%.1f\n",
              min(kept_counts), max(kept_counts), mean(kept_counts)))
  preds
}

run_enet_cv <- function(X_mat, y_vec, folds,
                        alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  preds <- numeric(length(y_vec))
  best_alphas <- numeric(length(folds))
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_along(y_vec), train_idx)
    tr_means  <- colMeans(X_mat[train_idx, , drop = FALSE])
    tr_sds    <- apply(X_mat[train_idx, , drop = FALSE], 2, sd)
    tr_sds[tr_sds == 0] <- 1
    X_tr <- scale(X_mat[train_idx, , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    X_te <- scale(X_mat[test_idx,  , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    best_alpha <- 0.5; best_cvm <- Inf
    for (a in alpha_grid) {
      set.seed(99)
      fit_inner <- cv.glmnet(X_tr, y_vec[train_idx],
                             alpha = a, family = "binomial", nfolds = 5)
      if (min(fit_inner$cvm) < best_cvm) {
        best_cvm <- min(fit_inner$cvm); best_alpha <- a
      }
    }
    best_alphas[i] <- best_alpha
    fit <- cv.glmnet(X_tr, y_vec[train_idx],
                     alpha = best_alpha, family = "binomial", nfolds = 5)
    preds[test_idx] <- predict(fit, X_te, s = "lambda.min",
                               type = "response")[, 1]
  }
  cat("Best alpha per fold:", round(best_alphas, 2), "\n")
  cat(sprintf("Mean alpha: %.3f (SD: %.3f)\n",
              mean(best_alphas), sd(best_alphas)))
  preds
}

evaluate <- function(preds, y_vec, label) {
  roc_obj    <- roc(y_vec, preds, quiet = TRUE)
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = c("threshold", "sensitivity", "specificity"))
  y_class <- ifelse(preds >= coords_obj$threshold, 1, 0)
  cm      <- confusionMatrix(factor(y_class), factor(y_vec), positive = "1")
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("AUC:               %.3f (95%% CI: %.3f-%.3f)\n",
              auc(roc_obj), ci(roc_obj)[1], ci(roc_obj)[3]))
  cat(sprintf("Balanced Accuracy: %.3f\n", cm$byClass["Balanced Accuracy"]))
  cat(sprintf("Sensitivity:       %.3f\n", cm$byClass["Sensitivity"]))
  cat(sprintf("Specificity:       %.3f\n", cm$byClass["Specificity"]))
  cat(sprintf("Kappa:             %.3f\n", cm$overall["Kappa"]))
  print(cm$table)
  invisible(roc_obj)
}

# ==============================================================================
# 4. CREATE FOLDS ONCE - used by ALL models
# ==============================================================================
set.seed(42)
folds <- createFolds(df_full$HEMISPHERIC, k = 10, returnTrain = TRUE)

# ==============================================================================
# 5. PCA COMPONENT SWEEP
# ==============================================================================
cat("\n====== PCA COMPONENT SWEEP ======\n")
sweep_k   <- c(5, 10, 15, 20, 25, 50, 75, 100)
sweep_k   <- sweep_k[sweep_k <= max_safe_k]
sweep_auc <- numeric(length(sweep_k))
sweep_bil <- numeric(length(sweep_k))

for (j in seq_along(sweep_k)) {
  n       <- sweep_k[j]
  pc_feat <- c(paste0("PC", 1:n), "Laterality_Quotient", "Edad")
  X_n     <- as.matrix(df_full[, pc_feat])
  preds   <- run_ridge_cv(X_n, y, folds)
  roc_obj <- roc(y, preds, quiet = TRUE)
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = "threshold")
  y_class  <- ifelse(preds >= coords_obj$threshold, 1, 0)
  bil_idx  <- which(df_full$HEMISPHERIC_2 == "Bilateral")
  bil_sens <- sum(y_class[bil_idx] == 0) / length(bil_idx)
  sweep_auc[j] <- round(as.numeric(auc(roc_obj)), 3)
  sweep_bil[j] <- round(bil_sens, 3)
  cat(sprintf("k = %3d: AUC = %.3f | Bilateral sens = %.3f\n",
              n, sweep_auc[j], sweep_bil[j]))
}

sweep_results <- data.frame(k = sweep_k, AUC = sweep_auc,
                            Bilateral_sensitivity = sweep_bil)
print(sweep_results)

ggplot(sweep_results, aes(x = k, y = AUC)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_point(data = sweep_results[sweep_results$k == 10, ],
             color = "red", size = 4) +
  scale_x_continuous(breaks = sweep_k) +
  labs(x = "Number of PCA components (k)", y = "Cross-validated AUC") +
  theme_minimal()
ggsave(paste0(BASE_PATH, "figures/pca_sweep.png"),
       width = 7, height = 4, dpi = 300)

# PC correlations with LI_GEN
cat("\nSpearman correlations: PC scores vs LI_GEN\n")
for (k in 1:10) {
  r <- cor.test(df_full[[paste0("PC", k)]], df_full$LI_GEN,
                method = "spearman")
  cat(sprintf("PC%d: rho = %.3f, p = %.4f\n", k, r$estimate, r$p.value))
}

# ==============================================================================
# 6. NESTED CV VALIDATION
# ==============================================================================
cat("\n====== NESTED CV VALIDATION ======\n")
set.seed(42)
outer_folds  <- createFolds(df_full$HEMISPHERIC, k = 10, returnTrain = TRUE)
k_grid       <- sweep_k
nested_preds <- numeric(length(y))
selected_k   <- numeric(length(outer_folds))

for (i in seq_along(outer_folds)) {
  outer_train <- outer_folds[[i]]
  outer_test  <- setdiff(seq_along(y), outer_train)
  inner_folds <- createFolds(y[outer_train], k = 5, returnTrain = TRUE)
  best_k <- 10; best_auc <- 0
  for (k in k_grid) {
    pc_feat <- c(paste0("PC", 1:k), "Laterality_Quotient", "Edad")
    X_k     <- as.matrix(df_full[outer_train, pc_feat])
    inner_preds <- numeric(length(outer_train))
    for (j in seq_along(inner_folds)) {
      inner_train <- inner_folds[[j]]
      inner_test  <- setdiff(seq_along(outer_train), inner_train)
      tr_means <- colMeans(X_k[inner_train, ])
      tr_sds   <- apply(X_k[inner_train, ], 2, sd)
      tr_sds[tr_sds == 0] <- 1
      X_tr <- scale(X_k[inner_train, ], center = tr_means, scale = tr_sds)
      X_te <- scale(X_k[inner_test,  ], center = tr_means, scale = tr_sds)
      fit  <- cv.glmnet(X_tr, y[outer_train][inner_train],
                        alpha = 0, family = "binomial", nfolds = 5)
      inner_preds[inner_test] <- predict(fit, X_te, s = "lambda.min",
                                         type = "response")[, 1]
    }
    k_auc <- as.numeric(auc(roc(y[outer_train], inner_preds, quiet = TRUE)))
    if (k_auc > best_auc) { best_auc <- k_auc; best_k <- k }
  }
  selected_k[i] <- best_k
  pc_feat  <- c(paste0("PC", 1:best_k), "Laterality_Quotient", "Edad")
  X_outer  <- as.matrix(df_full[, pc_feat])
  tr_means <- colMeans(X_outer[outer_train, ])
  tr_sds   <- apply(X_outer[outer_train, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- scale(X_outer[outer_train, ], center = tr_means, scale = tr_sds)
  X_te <- scale(X_outer[outer_test,  ], center = tr_means, scale = tr_sds)
  fit  <- cv.glmnet(X_tr, y[outer_train],
                    alpha = 0, family = "binomial", nfolds = 5)
  nested_preds[outer_test] <- predict(fit, X_te, s = "lambda.min",
                                      type = "response")[, 1]
}
roc_nested <- roc(y, nested_preds, quiet = TRUE)
cat(sprintf("Nested CV AUC: %.3f\n", auc(roc_nested)))
cat("\nSelected k per outer fold:\n")
print(table(selected_k))

# ==============================================================================
# 7. STRATEGY 1 MODELS
# ==============================================================================
cat("\n--- Config 1: Full AFC + Ridge ---\n")
preds_s1 <- run_ridge_cv(X_afc, y, folds)
roc_s1   <- evaluate(preds_s1, y, "Config 1: Full AFC + Ridge")

cat("\n--- Config 2: AFC PCA + Ridge ---\n")
preds_s2 <- run_ridge_cv(X_pca, y, folds)
roc_s2   <- evaluate(preds_s2, y, "Config 2: AFC PCA + Ridge")

cat("\n--- Config 3: AFC PCA + LQ + Ridge ---\n")
preds_s3 <- run_ridge_cv(X_pca_lq, y, folds)
roc_s3   <- evaluate(preds_s3, y, "Config 3: AFC PCA + LQ + Ridge")

cat("\n--- Config 4: AFC PCA + LQ + Age + Ridge (PRIMARY) ---\n")
preds_s4 <- run_ridge_cv(X_pca_lq_age, y, folds)
roc_s4   <- evaluate(preds_s4, y, "Config 4: AFC PCA + LQ + Age + Ridge")

cat("\n--- Config 5: Elastic Net PCA + LQ + Age ---\n")
set.seed(42)
preds_s5 <- run_enet_cv(X_pca_lq_age, y, folds)
roc_s5   <- evaluate(preds_s5, y, "Config 5: Elastic Net PCA + LQ + Age")

cat("\n--- Config 6: Random Forest PCA + LQ + Age ---\n")
preds_s6 <- numeric(length(y))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y), train_idx)
  train_df  <- data.frame(df_full[train_idx, feat_pca_lq_age],
                          y = factor(y[train_idx]))
  test_df   <- data.frame(df_full[test_idx,  feat_pca_lq_age])
  rf_fit    <- randomForest(y ~ ., data = train_df, ntree = 500,
                            importance = FALSE)
  preds_s6[test_idx] <- predict(rf_fit, test_df, type = "prob")[, "1"]
}
roc_s6 <- evaluate(preds_s6, y, "Config 6: Random Forest PCA + LQ + Age")

cat("\n--- Config 7: AFC PCA + LQ + Age + IrrLI + Ridge ---\n")
preds_s7 <- run_ridge_cv(X_pca_lq_age_irr, y, folds)
roc_s7   <- evaluate(preds_s7, y, "Config 7: AFC PCA + LQ + Age + IrrLI")

cat(sprintf("\nNested vs fixed-k difference: %.3f\n",
            auc(roc_s4) - auc(roc_nested)))

# Plain LR reference
cat("\n--- Plain LR (Config 4 features, reference) ---\n")
preds_glm <- numeric(length(y))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- as.data.frame(scale(X_pca_lq_age[train_idx, ],
                              center = tr_means, scale = tr_sds))
  X_te <- as.data.frame(scale(X_pca_lq_age[test_idx, ],
                              center = tr_means, scale = tr_sds))
  fit  <- glm(y_train ~ ., data = cbind(X_tr, y_train = y[train_idx]),
              family = binomial)
  preds_glm[test_idx] <- predict(fit, X_te, type = "response")
}
roc_glm <- roc(y, preds_glm, quiet = TRUE)
cat(sprintf("Plain LR AUC: %.3f\n", auc(roc_glm)))

# ==============================================================================
# 8. STRATEGY 1 DELONG TESTS
# ==============================================================================
cat("\n====== STRATEGY 1 DELONG TESTS ======\n")
cat("Config 4 vs Config 1:\n"); print(roc.test(roc_s4, roc_s1, method = "delong"))
cat("Config 4 vs Config 2:\n"); print(roc.test(roc_s4, roc_s2, method = "delong"))
cat("Config 4 vs Config 3:\n"); print(roc.test(roc_s4, roc_s3, method = "delong"))
cat("Config 4 vs Config 5:\n"); print(roc.test(roc_s4, roc_s5, method = "delong"))
cat("Config 4 vs Config 6:\n"); print(roc.test(roc_s4, roc_s6, method = "delong"))
cat("Config 4 vs Config 7:\n"); print(roc.test(roc_s4, roc_s7, method = "delong"))
cat("Config 4 vs Plain LR:\n"); print(roc.test(roc_s4, roc_glm, method = "delong"))

# ==============================================================================
# 9. PERMUTATION TEST
# ==============================================================================
cat("\n====== PERMUTATION TEST ======\n")
set.seed(42)
n_perm    <- 1000
perm_aucs <- numeric(n_perm)
for (p in 1:n_perm) {
  y_perm       <- sample(y)
  preds_perm   <- run_ridge_cv(X_pca_lq_age, y_perm, folds)
  perm_aucs[p] <- as.numeric(auc(roc(y_perm, preds_perm, quiet = TRUE)))
}
obs_auc <- as.numeric(auc(roc_s4))
cat(sprintf("Observed AUC:        %.3f\n", obs_auc))
cat(sprintf("Permutation p-value: %.4f\n", mean(perm_aucs >= obs_auc)))
cat(sprintf("Null AUC mean:       %.3f (SD: %.3f)\n",
            mean(perm_aucs), sd(perm_aucs)))

ggplot(data.frame(AUC = perm_aucs), aes(x = AUC)) +
  geom_histogram(bins = 40, fill = "#5369A0", color = "white", alpha = 0.8) +
  geom_vline(xintercept = obs_auc, color = "#D64E4E",
             linewidth = 1.2, linetype = "solid") +
  annotate("text", x = obs_auc + 0.005, y = Inf,
           label = sprintf("Observed AUC = %.3f\np < 0.001", obs_auc),
           hjust = 0, vjust = 1.4, size = 3.8, color = "#D64E4E") +
  labs(x = "Permuted AUC", y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
ggsave(paste0(BASE_PATH, "figures/permutation_null.png"),
       width = 6, height = 4, dpi = 300)

# ==============================================================================
# 10. CROSS-COHORT GENERALIZATION
# ==============================================================================
cat("\n====== CROSS-COHORT GENERALIZATION ======\n")
tr_means <- colMeans(X_pca_lq_age[ezm2_idx, ])
tr_sds   <- apply(X_pca_lq_age[ezm2_idx, ], 2, sd)
X_tr <- scale(X_pca_lq_age[ezm2_idx, ], center = tr_means, scale = tr_sds)
X_te <- scale(X_pca_lq_age[ezm3_idx, ], center = tr_means, scale = tr_sds)
fit_c1       <- cv.glmnet(X_tr, y[ezm2_idx], alpha = 0,
                          family = "binomial", nfolds = 5)
preds_cross1 <- predict(fit_c1, X_te, s = "lambda.min",
                        type = "response")[, 1]
roc_cross1   <- roc(y[ezm3_idx], preds_cross1, quiet = TRUE)
cat(sprintf("Train EZM2 -> Test EZM3: AUC = %.3f\n", auc(roc_cross1)))

tr_means2 <- colMeans(X_pca_lq_age[ezm3_idx, ])
tr_sds2   <- apply(X_pca_lq_age[ezm3_idx, ], 2, sd)
X_tr2 <- scale(X_pca_lq_age[ezm3_idx, ], center = tr_means2, scale = tr_sds2)
X_te2 <- scale(X_pca_lq_age[ezm2_idx, ], center = tr_means2, scale = tr_sds2)
fit_c2       <- cv.glmnet(X_tr2, y[ezm3_idx], alpha = 0,
                          family = "binomial", nfolds = 5)
preds_cross2 <- predict(fit_c2, X_te2, s = "lambda.min",
                        type = "response")[, 1]
roc_cross2   <- roc(y[ezm2_idx], preds_cross2, quiet = TRUE)
cat(sprintf("Train EZM3 -> Test EZM2: AUC = %.3f\n", auc(roc_cross2)))

# ==============================================================================
# 11. STRATEGY 2 MODELS - same folds object
# ==============================================================================
cat("\n--- Config 8: NodeAI + Ridge ---\n")
preds_s8 <- run_ridge_cv(X_nodeai, y, folds)
roc_s8   <- evaluate(preds_s8, y, "Config 8: NodeAI + Ridge")

cat("\n--- Config 9: NodeAI + LQ + Age + Ridge ---\n")
preds_s9 <- run_ridge_cv(X_nodeai_lq, y, folds)
roc_s9   <- evaluate(preds_s9, y, "Config 9: NodeAI + LQ + Age + Ridge")

cat("\n--- Config 10: NodeAI + Homotopic + LQ + Age + Ridge ---\n")
preds_s10 <- run_ridge_cv(X_conn_lq, y, folds)
roc_s10   <- evaluate(preds_s10, y,
                      "Config 10: NodeAI + Homotopic + LQ + Age + Ridge")

cat("\n--- Config 11: Graph Features + LQ + Age + Ridge ---\n")
preds_s11 <- run_ridge_cv(X_graph_lq, y, folds)
roc_s11   <- evaluate(preds_s11, y,
                      "Config 11: Graph Features + LQ + Age + Ridge")

cat("\n--- Config 12: All Interpretable + LQ + Age + Ridge ---\n")
preds_s12 <- run_ridge_cv(X_all_interp_lq, y, folds)
roc_s12   <- evaluate(preds_s12, y,
                      "Config 12: All Interpretable + LQ + Age + Ridge")

# ==============================================================================
# 12. STRATEGY 2 DELONG TESTS vs Config 4
# ==============================================================================
cat("\n====== STRATEGY 2 DELONG TESTS vs Config 4 ======\n")
cat("Config 9  vs Config 4:\n"); print(roc.test(roc_s9,  roc_s4, method = "delong"))
cat("Config 10 vs Config 4:\n"); print(roc.test(roc_s10, roc_s4, method = "delong"))
cat("Config 11 vs Config 4:\n"); print(roc.test(roc_s11, roc_s4, method = "delong"))
cat("Config 12 vs Config 4:\n"); print(roc.test(roc_s12, roc_s4, method = "delong"))
cat("Config 10 vs Config 9:\n"); print(roc.test(roc_s10, roc_s9, method = "delong"))

# ==============================================================================
# 13. SENSITIVITY ANALYSES
# ==============================================================================
cat("\n====== BILATERAL REMOVED ======\n")
no_bil_idx <- which(df_full$HEMISPHERIC_2 != "Bilateral")
df_nobi    <- df_full[no_bil_idx, ]
y_nobi     <- ifelse(df_nobi$HEMISPHERIC == "Typical", 1, 0)
X_nobi     <- as.matrix(df_nobi[, feat_pca_lq_age])
set.seed(42)
folds_nobi <- createFolds(df_nobi$HEMISPHERIC, k = 10, returnTrain = TRUE)
preds_nobi <- run_ridge_cv(X_nobi, y_nobi, folds_nobi)
roc_nobi   <- evaluate(preds_nobi, y_nobi,
                       "Bilateral removed: AFC PCA + LQ + Age + Ridge")

# Misclassification breakdown
coords_s4 <- coords(roc_s4, "best", best.method = "closest.topleft",
                    ret = "threshold")
y_class_s4 <- ifelse(preds_s4 >= coords_s4$threshold, 1, 0)
misclassified_atyp <- which(y == 0 & preds_s4 >= coords_s4$threshold)
cat("Misclassified atypical subjects:", length(misclassified_atyp), "\n")
print(table(df_full$HEMISPHERIC_2[misclassified_atyp]))

# Three-class classification
cat("\n====== THREE-CLASS CLASSIFICATION ======\n")
set.seed(42)
folds_3      <- createFolds(df_full$HEMISPHERIC_2, k = 10, returnTrain = TRUE)
preds_3class <- character(length(y_3class))
for (i in seq_along(folds_3)) {
  train_idx <- folds_3[[i]]
  test_idx  <- setdiff(seq_along(y_3class), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- as.data.frame(scale(X_pca_lq_age[train_idx, ],
                              center = tr_means, scale = tr_sds))
  X_te <- as.data.frame(scale(X_pca_lq_age[test_idx, ],
                              center = tr_means, scale = tr_sds))
  train_df <- cbind(X_tr, y = y_3class[train_idx])
  fit_3    <- multinom(y ~ ., data = train_df, trace = FALSE)
  preds_3class[test_idx] <- as.character(predict(fit_3, X_te))
}
cm_3 <- confusionMatrix(factor(preds_3class), factor(y_3class),
                        mode = "everything")
print(cm_3$table)
print(cm_3$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")])
cat(sprintf("Overall Accuracy: %.3f\n", cm_3$overall["Accuracy"]))
cat(sprintf("Kappa:            %.3f\n",   cm_3$overall["Kappa"]))

# ==============================================================================
# 14. CONTINUOUS LI PREDICTION
# ==============================================================================
cat("\n====== CONTINUOUS LI PREDICTION ======\n")
y_cont     <- df_full$LI_GEN
preds_cont <- numeric(length(y_cont))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y_cont), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- scale(X_pca_lq_age[train_idx, ], center = tr_means, scale = tr_sds)
  X_te <- scale(X_pca_lq_age[test_idx,  ], center = tr_means, scale = tr_sds)
  fit  <- cv.glmnet(X_tr, y_cont[train_idx],
                    alpha = 0, family = "gaussian", nfolds = 5)
  preds_cont[test_idx] <- predict(fit, X_te, s = "lambda.min")[, 1]
}
r_cont <- cor.test(preds_cont, y_cont, method = "spearman")
cat(sprintf("Spearman rho = %.3f, p = %.2e\n",
            r_cont$estimate, r_cont$p.value))
print(tapply(preds_cont, df_full$HEMISPHERIC_2, median))
kw <- kruskal.test(preds_cont ~ df_full$HEMISPHERIC_2)
cat(sprintf("Kruskal-Wallis: chi2 = %.3f, p = %.2e\n",
            kw$statistic, kw$p.value))
print(pairwise.wilcox.test(preds_cont, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))

plot_df <- data.frame(Observed  = y_cont, Predicted = preds_cont,
                      Group = factor(df_full$HEMISPHERIC_2,
                                     levels = c("Left", "Right", "Bilateral")))
group_colors <- c("Left" = "#D64E4E", "Right" = "#2B81AC",
                  "Bilateral" = "#7A51A0")
ggplot(plot_df, aes(x = Observed, y = Predicted, color = Group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40",
              linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Spearman \u03c1 = %.3f, p < 0.001",
                           r_cont$estimate),
           hjust = 1.05, vjust = -0.8, size = 3.8, color = "grey30") +
  labs(x = "Observed LI (task-based)",
       y = "Predicted LI (AFC PCA + LQ + Age)",
       color = "Lateralization group") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
ggsave(paste0(BASE_PATH, "figures/li_prediction_scatter.png"),
       width = 6, height = 5, dpi = 300)

# ==============================================================================
# 15. GROUP DIFFERENCE ANALYSES
# ==============================================================================
cat("\n====== GROUP DIFFERENCE ANALYSES ======\n")

demo_df <- merge(demo[, c("ID", "HEMISPHERIC", "Edad", "Sexo",
                          "Laterality_Quotient")],
                 data.frame(ID = df_full$ID), by = "ID")
cat("Age:\n")
print(wilcox.test(Edad ~ HEMISPHERIC, data = demo_df))
print(tapply(demo_df$Edad, demo_df$HEMISPHERIC, median))
cat("\nSex:\n")
print(chisq.test(table(demo_df$Sexo, demo_df$HEMISPHERIC)))
print(prop.table(table(demo_df$Sexo, demo_df$HEMISPHERIC), margin = 2))
cat("\nLaterality Quotient:\n")
print(wilcox.test(Laterality_Quotient ~ HEMISPHERIC, data = demo_df))
print(tapply(demo_df$Laterality_Quotient, demo_df$HEMISPHERIC, median))

cat("\n--- AFC score ---\n")
print(wilcox.test(df_full$afc_score ~ df_full$HEMISPHERIC))
print(kruskal.test(df_full$afc_score ~ df_full$HEMISPHERIC_2))
print(pairwise.wilcox.test(df_full$afc_score, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))
print(tapply(df_full$afc_score, df_full$HEMISPHERIC_2,
             function(x) c(mean = round(mean(x), 3), sd = round(sd(x), 3))))

cat("\n--- IrrLI ---\n")
print(cor.test(df_full$IrrLI, df_full$LI_GEN, method = "spearman"))
print(kruskal.test(df_full$IrrLI ~ df_full$HEMISPHERIC_2))
print(pairwise.wilcox.test(df_full$IrrLI, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))

# ==============================================================================
# 16. SUBGROUP ACCURACY
# ==============================================================================
cat("\n====== CONFIG 4 ACCURACY BY SUBGROUP ======\n")
for (cls in c("Left", "Right", "Bilateral")) {
  idx        <- which(df_full$HEMISPHERIC_2 == cls)
  true_label <- ifelse(cls == "Left", 1, 0)
  correct    <- sum(y_class_s4[idx] == true_label)
  total      <- length(idx)
  cat(sprintf("%-12s: %d / %d correctly predicted (%.1f%%)\n",
              cls, correct, total, (correct / total) * 100))
}

cat("\n====== CONFIG 9 ACCURACY BY SUBGROUP ======\n")
coords_s9  <- coords(roc_s9, "best", best.method = "closest.topleft",
                     ret = "threshold")
y_class_s9 <- ifelse(preds_s9 >= coords_s9$threshold, 1, 0)
for (cls in c("Left", "Right", "Bilateral")) {
  idx        <- which(df_full$HEMISPHERIC_2 == cls)
  true_label <- ifelse(cls == "Left", 1, 0)
  correct    <- sum(y_class_s9[idx] == true_label)
  total      <- length(idx)
  cat(sprintf("%-12s: %d / %d correctly predicted (%.1f%%)\n",
              cls, correct, total, (correct / total) * 100))
}

# ==============================================================================
# 17. INTERPRETABILITY: NODEAI COEFFICIENTS
# ==============================================================================
cat("\n====== NODEAI REGION COEFFICIENTS ======\n")
tr_means <- colMeans(X_nodeai_lq)
tr_sds   <- apply(X_nodeai_lq, 2, sd); tr_sds[tr_sds == 0] <- 1
X_scaled <- scale(X_nodeai_lq, center = tr_means, scale = tr_sds)
fit_final <- cv.glmnet(X_scaled, y, alpha = 0, family = "binomial", nfolds = 10)
coefs     <- coef(fit_final, s = "lambda.min")
region_lookup <- read.csv(paste0(BASE_PATH, "region_lookup_100.csv"))
coef_df <- data.frame(Feature = rownames(coefs)[-1],
                      Coef    = as.numeric(coefs[-1]))
coef_df <- coef_df[!coef_df$Feature %in% c("Laterality_Quotient", "Edad"), ]
feature_nums         <- as.integer(gsub("nodeAI_R", "", coef_df$Feature))
positions            <- match(feature_nums, roi_idx_simple)
coef_df$Abbreviation <- region_lookup$Abbreviation[positions]
coef_df$Full.Name    <- region_lookup$Full.Name[positions]
coef_df$Atlas        <- region_lookup$Atlas[positions]
coef_df$Hemisphere   <- region_lookup$Hemisphere[positions]
coef_df <- coef_df[order(abs(coef_df$Coef), decreasing = TRUE), ]
print(head(coef_df[, c("Abbreviation", "Full.Name", "Atlas",
                       "Hemisphere", "Coef")], 20))

# ==============================================================================
# 18. SUMMARY TABLES
# ==============================================================================
cat("\n====== SUMMARY TABLE: ALL CONFIGS ======\n")
get_metrics <- function(roc_obj, preds, y_vec) {
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = c("threshold", "sensitivity", "specificity"))
  y_class <- ifelse(preds >= coords_obj$threshold, 1, 0)
  cm      <- confusionMatrix(factor(y_class), factor(y_vec), positive = "1")
  list(AUC          = round(as.numeric(auc(roc_obj)), 3),
       Balanced_Acc = round(as.numeric(cm$byClass["Balanced Accuracy"]), 3),
       Sensitivity  = round(as.numeric(cm$byClass["Sensitivity"]), 3),
       Specificity  = round(as.numeric(cm$byClass["Specificity"]), 3),
       Kappa        = round(as.numeric(cm$overall["Kappa"]), 3))
}

roc_list  <- list(roc_s1, roc_s2, roc_s3, roc_s4, roc_s5, roc_s6, roc_s7,
                  roc_s8, roc_s9, roc_s10, roc_s11, roc_s12)
pred_list <- list(preds_s1, preds_s2, preds_s3, preds_s4, preds_s5,
                  preds_s6, preds_s7, preds_s8, preds_s9, preds_s10,
                  preds_s11, preds_s12)
config_labels <- c(
  "1: Full AFC + Ridge",
  "2: AFC PCA + Ridge",
  "3: AFC PCA + LQ + Ridge",
  "4: AFC PCA + LQ + Age + Ridge (PRIMARY)",
  "5: Elastic Net PCA + LQ + Age",
  "6: Random Forest PCA + LQ + Age",
  "7: AFC PCA + LQ + Age + IrrLI + Ridge",
  "8: NodeAI + Ridge",
  "9: NodeAI + LQ + Age + Ridge",
  "10: NodeAI + Homotopic + LQ + Age + Ridge",
  "11: Graph Features + LQ + Age + Ridge",
  "12: All Interpretable + LQ + Age + Ridge")

metrics_list <- mapply(get_metrics, roc_list, pred_list,
                       MoreArgs = list(y_vec = y), SIMPLIFY = FALSE)
results_table <- data.frame(
  Config       = config_labels,
  AUC          = sapply(metrics_list, `[[`, "AUC"),
  Balanced_Acc = sapply(metrics_list, `[[`, "Balanced_Acc"),
  Sensitivity  = sapply(metrics_list, `[[`, "Sensitivity"),
  Specificity  = sapply(metrics_list, `[[`, "Specificity"),
  Kappa        = sapply(metrics_list, `[[`, "Kappa"))
print(results_table)

cat(sprintf("\nNested CV AUC:            %.3f\n", auc(roc_nested)))
cat(sprintf("Plain LR (Config 4):      %.3f\n", auc(roc_glm)))
cat(sprintf("Cross-cohort EZM2->EZM3:  %.3f\n", auc(roc_cross1)))
cat(sprintf("Cross-cohort EZM3->EZM2:  %.3f\n", auc(roc_cross2)))
cat(sprintf("Bilateral removed:        %.3f\n", auc(roc_nobi)))





# ==============================================================================
# 19. PCA BACK-PROJECTION: CONNECTION-LEVEL IMPORTANCE
# ==============================================================================
cat("\n====== PCA BACK-PROJECTION ======\n")

# Recompute PCA on full AFC matrix (same data used for scores)
pca_obj <- prcomp(X_afc, center = TRUE, scale. = TRUE)

# Fit ridge on full dataset using Config 4 features to extract coefficients
tr_means_bp <- colMeans(X_pca_lq_age)
tr_sds_bp   <- apply(X_pca_lq_age, 2, sd); tr_sds_bp[tr_sds_bp == 0] <- 1
X_scaled_bp <- scale(X_pca_lq_age, center = tr_means_bp, scale = tr_sds_bp)
fit_full    <- cv.glmnet(X_scaled_bp, y, alpha = 0,
                         family = "binomial", nfolds = 10)

# Extract only the PC1-PC10 coefficients (exclude LQ, Age intercept)
all_coefs   <- as.numeric(coef(fit_full, s = "lambda.min"))
# coef order: intercept, PC1..PC10, Laterality_Quotient, Edad
ridge_coefs <- all_coefs[2:11]  # PC1 through PC10

# Project back onto connection space
connection_weights <- pca_obj$rotation[, 1:10] %*% ridge_coefs

weight_df <- data.frame(
  connection = afc_cols,
  weight     = as.numeric(connection_weights)
)
weight_df <- weight_df[order(abs(weight_df$weight), decreasing = TRUE), ]

cat("Top 20 connections by absolute weight:\n")
print(head(weight_df, 20))

write.csv(weight_df, paste0(BASE_PATH, "pca_backprojection_weights.csv"),
          row.names = FALSE)
run_enet_cv <- function(X_mat, y_vec, folds,
                        alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  preds <- numeric(length(y_vec))
  best_alphas <- numeric(length(folds))
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_along(y_vec), train_idx)
    tr_means  <- colMeans(X_mat[train_idx, , drop = FALSE])
    tr_sds    <- apply(X_mat[train_idx, , drop = FALSE], 2, sd)
    tr_sds[tr_sds == 0] <- 1
    X_tr <- scale(X_mat[train_idx, , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    X_te <- scale(X_mat[test_idx,  , drop = FALSE],
                  center = tr_means, scale = tr_sds)
    best_alpha <- 0.5; best_cvm <- Inf
    for (a in alpha_grid) {
      set.seed(99)
      fit_inner <- cv.glmnet(X_tr, y_vec[train_idx],
                             alpha = a, family = "binomial", nfolds = 5)
      if (min(fit_inner$cvm) < best_cvm) {
        best_cvm <- min(fit_inner$cvm); best_alpha <- a
      }
    }
    best_alphas[i] <- best_alpha
    fit <- cv.glmnet(X_tr, y_vec[train_idx],
                     alpha = best_alpha, family = "binomial", nfolds = 5)
    preds[test_idx] <- predict(fit, X_te, s = "lambda.min",
                               type = "response")[, 1]
  }
  cat("Best alpha per fold:", round(best_alphas, 2), "\n")
  cat(sprintf("Mean alpha: %.3f (SD: %.3f)\n",
              mean(best_alphas), sd(best_alphas)))
  preds
}

evaluate <- function(preds, y_vec, label) {
  roc_obj    <- roc(y_vec, preds, quiet = TRUE)
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = c("threshold", "sensitivity", "specificity"))
  y_class <- ifelse(preds >= coords_obj$threshold, 1, 0)
  cm      <- confusionMatrix(factor(y_class), factor(y_vec), positive = "1")
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("AUC:               %.3f (95%% CI: %.3f-%.3f)\n",
              auc(roc_obj), ci(roc_obj)[1], ci(roc_obj)[3]))
  cat(sprintf("Balanced Accuracy: %.3f\n", cm$byClass["Balanced Accuracy"]))
  cat(sprintf("Sensitivity:       %.3f\n", cm$byClass["Sensitivity"]))
  cat(sprintf("Specificity:       %.3f\n", cm$byClass["Specificity"]))
  cat(sprintf("Kappa:             %.3f\n", cm$overall["Kappa"]))
  print(cm$table)
  invisible(roc_obj)
}

# ==============================================================================
# 4. CREATE FOLDS ONCE - used by ALL models
# ==============================================================================
set.seed(42)
folds <- createFolds(df_full$HEMISPHERIC, k = 10, returnTrain = TRUE)

# ==============================================================================
# 5. PCA COMPONENT SWEEP
# ==============================================================================
cat("\n====== PCA COMPONENT SWEEP ======\n")
sweep_k   <- c(5, 10, 15, 20, 25, 50, 75, 100)
sweep_k   <- sweep_k[sweep_k <= max_safe_k]
sweep_auc <- numeric(length(sweep_k))
sweep_bil <- numeric(length(sweep_k))

for (j in seq_along(sweep_k)) {
  n       <- sweep_k[j]
  pc_feat <- c(paste0("PC", 1:n), "Laterality_Quotient", "Edad")
  X_n     <- as.matrix(df_full[, pc_feat])
  preds   <- run_ridge_cv(X_n, y, folds)
  roc_obj <- roc(y, preds, quiet = TRUE)
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = "threshold")
  y_class  <- ifelse(preds >= coords_obj$threshold, 1, 0)
  bil_idx  <- which(df_full$HEMISPHERIC_2 == "Bilateral")
  bil_sens <- sum(y_class[bil_idx] == 0) / length(bil_idx)
  sweep_auc[j] <- round(as.numeric(auc(roc_obj)), 3)
  sweep_bil[j] <- round(bil_sens, 3)
  cat(sprintf("k = %3d: AUC = %.3f | Bilateral sens = %.3f\n",
              n, sweep_auc[j], sweep_bil[j]))
}

sweep_results <- data.frame(k = sweep_k, AUC = sweep_auc,
                            Bilateral_sensitivity = sweep_bil)
print(sweep_results)

ggplot(sweep_results, aes(x = k, y = AUC)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_point(data = sweep_results[sweep_results$k == 10, ],
             color = "red", size = 4) +
  scale_x_continuous(breaks = sweep_k) +
  labs(x = "Number of PCA components (k)", y = "Cross-validated AUC") +
  theme_minimal()
ggsave(paste0(BASE_PATH, "figures/pca_sweep.png"),
       width = 7, height = 4, dpi = 300)

# PC correlations with LI_GEN
cat("\nSpearman correlations: PC scores vs LI_GEN\n")
for (k in 1:10) {
  r <- cor.test(df_full[[paste0("PC", k)]], df_full$LI_GEN,
                method = "spearman")
  cat(sprintf("PC%d: rho = %.3f, p = %.4f\n", k, r$estimate, r$p.value))
}

# ==============================================================================
# 6. NESTED CV VALIDATION
# ==============================================================================
cat("\n====== NESTED CV VALIDATION ======\n")
set.seed(42)
outer_folds  <- createFolds(df_full$HEMISPHERIC, k = 10, returnTrain = TRUE)
k_grid       <- sweep_k
nested_preds <- numeric(length(y))
selected_k   <- numeric(length(outer_folds))

for (i in seq_along(outer_folds)) {
  outer_train <- outer_folds[[i]]
  outer_test  <- setdiff(seq_along(y), outer_train)
  inner_folds <- createFolds(y[outer_train], k = 5, returnTrain = TRUE)
  best_k <- 10; best_auc <- 0
  for (k in k_grid) {
    pc_feat <- c(paste0("PC", 1:k), "Laterality_Quotient", "Edad")
    X_k     <- as.matrix(df_full[outer_train, pc_feat])
    inner_preds <- numeric(length(outer_train))
    for (j in seq_along(inner_folds)) {
      inner_train <- inner_folds[[j]]
      inner_test  <- setdiff(seq_along(outer_train), inner_train)
      tr_means <- colMeans(X_k[inner_train, ])
      tr_sds   <- apply(X_k[inner_train, ], 2, sd)
      tr_sds[tr_sds == 0] <- 1
      X_tr <- scale(X_k[inner_train, ], center = tr_means, scale = tr_sds)
      X_te <- scale(X_k[inner_test,  ], center = tr_means, scale = tr_sds)
      fit  <- cv.glmnet(X_tr, y[outer_train][inner_train],
                        alpha = 0, family = "binomial", nfolds = 5)
      inner_preds[inner_test] <- predict(fit, X_te, s = "lambda.min",
                                         type = "response")[, 1]
    }
    k_auc <- as.numeric(auc(roc(y[outer_train], inner_preds, quiet = TRUE)))
    if (k_auc > best_auc) { best_auc <- k_auc; best_k <- k }
  }
  selected_k[i] <- best_k
  pc_feat  <- c(paste0("PC", 1:best_k), "Laterality_Quotient", "Edad")
  X_outer  <- as.matrix(df_full[, pc_feat])
  tr_means <- colMeans(X_outer[outer_train, ])
  tr_sds   <- apply(X_outer[outer_train, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- scale(X_outer[outer_train, ], center = tr_means, scale = tr_sds)
  X_te <- scale(X_outer[outer_test,  ], center = tr_means, scale = tr_sds)
  fit  <- cv.glmnet(X_tr, y[outer_train],
                    alpha = 0, family = "binomial", nfolds = 5)
  nested_preds[outer_test] <- predict(fit, X_te, s = "lambda.min",
                                      type = "response")[, 1]
}
roc_nested <- roc(y, nested_preds, quiet = TRUE)
cat(sprintf("Nested CV AUC: %.3f\n", auc(roc_nested)))
cat("\nSelected k per outer fold:\n")
print(table(selected_k))

# ==============================================================================
# 7. STRATEGY 1 MODELS
# ==============================================================================
cat("\n--- Config 1: Full AFC + Ridge ---\n")
preds_s1 <- run_ridge_cv(X_afc, y, folds)
roc_s1   <- evaluate(preds_s1, y, "Config 1: Full AFC + Ridge")

cat("\n--- Config 2: AFC PCA + Ridge ---\n")
preds_s2 <- run_ridge_cv(X_pca, y, folds)
roc_s2   <- evaluate(preds_s2, y, "Config 2: AFC PCA + Ridge")

cat("\n--- Config 3: AFC PCA + LQ + Ridge ---\n")
preds_s3 <- run_ridge_cv(X_pca_lq, y, folds)
roc_s3   <- evaluate(preds_s3, y, "Config 3: AFC PCA + LQ + Ridge")

cat("\n--- Config 4: AFC PCA + LQ + Age + Ridge (PRIMARY) ---\n")
preds_s4 <- run_ridge_cv(X_pca_lq_age, y, folds)
roc_s4   <- evaluate(preds_s4, y, "Config 4: AFC PCA + LQ + Age + Ridge")

cat("\n--- Config 5: Elastic Net PCA + LQ + Age ---\n")
set.seed(42)
preds_s5 <- run_enet_cv(X_pca_lq_age, y, folds)
roc_s5   <- evaluate(preds_s5, y, "Config 5: Elastic Net PCA + LQ + Age")

cat("\n--- Config 6: Random Forest PCA + LQ + Age ---\n")
preds_s6 <- numeric(length(y))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y), train_idx)
  train_df  <- data.frame(df_full[train_idx, feat_pca_lq_age],
                          y = factor(y[train_idx]))
  test_df   <- data.frame(df_full[test_idx,  feat_pca_lq_age])
  rf_fit    <- randomForest(y ~ ., data = train_df, ntree = 500,
                            importance = FALSE)
  preds_s6[test_idx] <- predict(rf_fit, test_df, type = "prob")[, "1"]
}
roc_s6 <- evaluate(preds_s6, y, "Config 6: Random Forest PCA + LQ + Age")

cat("\n--- Config 7: AFC PCA + LQ + Age + IrrLI + Ridge ---\n")
preds_s7 <- run_ridge_cv(X_pca_lq_age_irr, y, folds)
roc_s7   <- evaluate(preds_s7, y, "Config 7: AFC PCA + LQ + Age + IrrLI")

cat(sprintf("\nNested vs fixed-k difference: %.3f\n",
            auc(roc_s4) - auc(roc_nested)))

# Plain LR reference
cat("\n--- Plain LR (Config 4 features, reference) ---\n")
preds_glm <- numeric(length(y))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- as.data.frame(scale(X_pca_lq_age[train_idx, ],
                              center = tr_means, scale = tr_sds))
  X_te <- as.data.frame(scale(X_pca_lq_age[test_idx, ],
                              center = tr_means, scale = tr_sds))
  fit  <- glm(y_train ~ ., data = cbind(X_tr, y_train = y[train_idx]),
              family = binomial)
  preds_glm[test_idx] <- predict(fit, X_te, type = "response")
}
roc_glm <- roc(y, preds_glm, quiet = TRUE)
cat(sprintf("Plain LR AUC: %.3f\n", auc(roc_glm)))

# ==============================================================================
# 8. STRATEGY 1 DELONG TESTS
# ==============================================================================
cat("\n====== STRATEGY 1 DELONG TESTS ======\n")
cat("Config 4 vs Config 1:\n"); print(roc.test(roc_s4, roc_s1, method = "delong"))
cat("Config 4 vs Config 2:\n"); print(roc.test(roc_s4, roc_s2, method = "delong"))
cat("Config 4 vs Config 3:\n"); print(roc.test(roc_s4, roc_s3, method = "delong"))
cat("Config 4 vs Config 5:\n"); print(roc.test(roc_s4, roc_s5, method = "delong"))
cat("Config 4 vs Config 6:\n"); print(roc.test(roc_s4, roc_s6, method = "delong"))
cat("Config 4 vs Config 7:\n"); print(roc.test(roc_s4, roc_s7, method = "delong"))
cat("Config 4 vs Plain LR:\n"); print(roc.test(roc_s4, roc_glm, method = "delong"))

# ==============================================================================
# 9. PERMUTATION TEST
# ==============================================================================
cat("\n====== PERMUTATION TEST ======\n")
set.seed(42)
n_perm    <- 1000
perm_aucs <- numeric(n_perm)
for (p in 1:n_perm) {
  y_perm       <- sample(y)
  preds_perm   <- run_ridge_cv(X_pca_lq_age, y_perm, folds)
  perm_aucs[p] <- as.numeric(auc(roc(y_perm, preds_perm, quiet = TRUE)))
}
obs_auc <- as.numeric(auc(roc_s4))
cat(sprintf("Observed AUC:        %.3f\n", obs_auc))
cat(sprintf("Permutation p-value: %.4f\n", mean(perm_aucs >= obs_auc)))
cat(sprintf("Null AUC mean:       %.3f (SD: %.3f)\n",
            mean(perm_aucs), sd(perm_aucs)))

ggplot(data.frame(AUC = perm_aucs), aes(x = AUC)) +
  geom_histogram(bins = 40, fill = "#5369A0", color = "white", alpha = 0.8) +
  geom_vline(xintercept = obs_auc, color = "#D64E4E",
             linewidth = 1.2, linetype = "solid") +
  annotate("text", x = obs_auc + 0.005, y = Inf,
           label = sprintf("Observed AUC = %.3f\np < 0.001", obs_auc),
           hjust = 0, vjust = 1.4, size = 3.8, color = "#D64E4E") +
  labs(x = "Permuted AUC", y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
ggsave(paste0(BASE_PATH, "figures/permutation_null.png"),
       width = 6, height = 4, dpi = 300)

# ==============================================================================
# 10. CROSS-COHORT GENERALIZATION
# ==============================================================================
cat("\n====== CROSS-COHORT GENERALIZATION ======\n")
tr_means <- colMeans(X_pca_lq_age[ezm2_idx, ])
tr_sds   <- apply(X_pca_lq_age[ezm2_idx, ], 2, sd)
X_tr <- scale(X_pca_lq_age[ezm2_idx, ], center = tr_means, scale = tr_sds)
X_te <- scale(X_pca_lq_age[ezm3_idx, ], center = tr_means, scale = tr_sds)
fit_c1       <- cv.glmnet(X_tr, y[ezm2_idx], alpha = 0,
                          family = "binomial", nfolds = 5)
preds_cross1 <- predict(fit_c1, X_te, s = "lambda.min",
                        type = "response")[, 1]
roc_cross1   <- roc(y[ezm3_idx], preds_cross1, quiet = TRUE)
cat(sprintf("Train EZM2 -> Test EZM3: AUC = %.3f\n", auc(roc_cross1)))

tr_means2 <- colMeans(X_pca_lq_age[ezm3_idx, ])
tr_sds2   <- apply(X_pca_lq_age[ezm3_idx, ], 2, sd)
X_tr2 <- scale(X_pca_lq_age[ezm3_idx, ], center = tr_means2, scale = tr_sds2)
X_te2 <- scale(X_pca_lq_age[ezm2_idx, ], center = tr_means2, scale = tr_sds2)
fit_c2       <- cv.glmnet(X_tr2, y[ezm3_idx], alpha = 0,
                          family = "binomial", nfolds = 5)
preds_cross2 <- predict(fit_c2, X_te2, s = "lambda.min",
                        type = "response")[, 1]
roc_cross2   <- roc(y[ezm2_idx], preds_cross2, quiet = TRUE)
cat(sprintf("Train EZM3 -> Test EZM2: AUC = %.3f\n", auc(roc_cross2)))

# ==============================================================================
# 11. STRATEGY 2 MODELS - same folds object
# ==============================================================================
cat("\n--- Config 8: NodeAI + Ridge ---\n")
preds_s8 <- run_ridge_cv(X_nodeai, y, folds)
roc_s8   <- evaluate(preds_s8, y, "Config 8: NodeAI + Ridge")

cat("\n--- Config 9: NodeAI + LQ + Age + Ridge ---\n")
preds_s9 <- run_ridge_cv(X_nodeai_lq, y, folds)
roc_s9   <- evaluate(preds_s9, y, "Config 9: NodeAI + LQ + Age + Ridge")

cat("\n--- Config 10: NodeAI + Homotopic + LQ + Age + Ridge ---\n")
preds_s10 <- run_ridge_cv(X_conn_lq, y, folds)
roc_s10   <- evaluate(preds_s10, y,
                      "Config 10: NodeAI + Homotopic + LQ + Age + Ridge")

cat("\n--- Config 11: Graph Features + LQ + Age + Ridge ---\n")
preds_s11 <- run_ridge_cv(X_graph_lq, y, folds)
roc_s11   <- evaluate(preds_s11, y,
                      "Config 11: Graph Features + LQ + Age + Ridge")

cat("\n--- Config 12: All Interpretable + LQ + Age + Ridge ---\n")
preds_s12 <- run_ridge_cv(X_all_interp_lq, y, folds)
roc_s12   <- evaluate(preds_s12, y,
                      "Config 12: All Interpretable + LQ + Age + Ridge")

# ==============================================================================
# 12. STRATEGY 2 DELONG TESTS vs Config 4
# ==============================================================================
cat("\n====== STRATEGY 2 DELONG TESTS vs Config 4 ======\n")
cat("Config 9  vs Config 4:\n"); print(roc.test(roc_s9,  roc_s4, method = "delong"))
cat("Config 10 vs Config 4:\n"); print(roc.test(roc_s10, roc_s4, method = "delong"))
cat("Config 11 vs Config 4:\n"); print(roc.test(roc_s11, roc_s4, method = "delong"))
cat("Config 12 vs Config 4:\n"); print(roc.test(roc_s12, roc_s4, method = "delong"))
cat("Config 10 vs Config 9:\n"); print(roc.test(roc_s10, roc_s9, method = "delong"))

# ==============================================================================
# 13. SENSITIVITY ANALYSES
# ==============================================================================
cat("\n====== BILATERAL REMOVED ======\n")
no_bil_idx <- which(df_full$HEMISPHERIC_2 != "Bilateral")
df_nobi    <- df_full[no_bil_idx, ]
y_nobi     <- ifelse(df_nobi$HEMISPHERIC == "Typical", 1, 0)
X_nobi     <- as.matrix(df_nobi[, feat_pca_lq_age])
set.seed(42)
folds_nobi <- createFolds(df_nobi$HEMISPHERIC, k = 10, returnTrain = TRUE)
preds_nobi <- run_ridge_cv(X_nobi, y_nobi, folds_nobi)
roc_nobi   <- evaluate(preds_nobi, y_nobi,
                       "Bilateral removed: AFC PCA + LQ + Age + Ridge")

# Misclassification breakdown
coords_s4 <- coords(roc_s4, "best", best.method = "closest.topleft",
                    ret = "threshold")
y_class_s4 <- ifelse(preds_s4 >= coords_s4$threshold, 1, 0)
misclassified_atyp <- which(y == 0 & preds_s4 >= coords_s4$threshold)
cat("Misclassified atypical subjects:", length(misclassified_atyp), "\n")
print(table(df_full$HEMISPHERIC_2[misclassified_atyp]))

# Three-class classification
cat("\n====== THREE-CLASS CLASSIFICATION ======\n")
set.seed(42)
folds_3      <- createFolds(df_full$HEMISPHERIC_2, k = 10, returnTrain = TRUE)
preds_3class <- character(length(y_3class))
for (i in seq_along(folds_3)) {
  train_idx <- folds_3[[i]]
  test_idx  <- setdiff(seq_along(y_3class), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- as.data.frame(scale(X_pca_lq_age[train_idx, ],
                              center = tr_means, scale = tr_sds))
  X_te <- as.data.frame(scale(X_pca_lq_age[test_idx, ],
                              center = tr_means, scale = tr_sds))
  train_df <- cbind(X_tr, y = y_3class[train_idx])
  fit_3    <- multinom(y ~ ., data = train_df, trace = FALSE)
  preds_3class[test_idx] <- as.character(predict(fit_3, X_te))
}
cm_3 <- confusionMatrix(factor(preds_3class), factor(y_3class),
                        mode = "everything")
print(cm_3$table)
print(cm_3$byClass[, c("Sensitivity", "Specificity", "Balanced Accuracy")])
cat(sprintf("Overall Accuracy: %.3f\n", cm_3$overall["Accuracy"]))
cat(sprintf("Kappa:            %.3f\n",   cm_3$overall["Kappa"]))

# ==============================================================================
# 14. CONTINUOUS LI PREDICTION
# ==============================================================================
cat("\n====== CONTINUOUS LI PREDICTION ======\n")
y_cont     <- df_full$LI_GEN
preds_cont <- numeric(length(y_cont))
for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  test_idx  <- setdiff(seq_along(y_cont), train_idx)
  tr_means  <- colMeans(X_pca_lq_age[train_idx, ])
  tr_sds    <- apply(X_pca_lq_age[train_idx, ], 2, sd)
  tr_sds[tr_sds == 0] <- 1
  X_tr <- scale(X_pca_lq_age[train_idx, ], center = tr_means, scale = tr_sds)
  X_te <- scale(X_pca_lq_age[test_idx,  ], center = tr_means, scale = tr_sds)
  fit  <- cv.glmnet(X_tr, y_cont[train_idx],
                    alpha = 0, family = "gaussian", nfolds = 5)
  preds_cont[test_idx] <- predict(fit, X_te, s = "lambda.min")[, 1]
}
r_cont <- cor.test(preds_cont, y_cont, method = "spearman")
cat(sprintf("Spearman rho = %.3f, p = %.2e\n",
            r_cont$estimate, r_cont$p.value))
print(tapply(preds_cont, df_full$HEMISPHERIC_2, median))
kw <- kruskal.test(preds_cont ~ df_full$HEMISPHERIC_2)
cat(sprintf("Kruskal-Wallis: chi2 = %.3f, p = %.2e\n",
            kw$statistic, kw$p.value))
print(pairwise.wilcox.test(preds_cont, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))

plot_df <- data.frame(Observed  = y_cont, Predicted = preds_cont,
                      Group = factor(df_full$HEMISPHERIC_2,
                                     levels = c("Left", "Right", "Bilateral")))
group_colors <- c("Left" = "#D64E4E", "Right" = "#2B81AC",
                  "Bilateral" = "#7A51A0")
ggplot(plot_df, aes(x = Observed, y = Predicted, color = Group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40",
              linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = group_colors) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Spearman \u03c1 = %.3f, p < 0.001",
                           r_cont$estimate),
           hjust = 1.05, vjust = -0.8, size = 3.8, color = "grey30") +
  labs(x = "Observed LI (task-based)",
       y = "Predicted LI (AFC PCA + LQ + Age)",
       color = "Lateralization group") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
ggsave(paste0(BASE_PATH, "figures/li_prediction_scatter.png"),
       width = 6, height = 5, dpi = 300)

# ==============================================================================
# 15. GROUP DIFFERENCE ANALYSES
# ==============================================================================
cat("\n====== GROUP DIFFERENCE ANALYSES ======\n")

demo_df <- merge(demo[, c("ID", "HEMISPHERIC", "Edad", "Sexo",
                          "Laterality_Quotient")],
                 data.frame(ID = df_full$ID), by = "ID")
cat("Age:\n")
print(wilcox.test(Edad ~ HEMISPHERIC, data = demo_df))
print(tapply(demo_df$Edad, demo_df$HEMISPHERIC, median))
cat("\nSex:\n")
print(chisq.test(table(demo_df$Sexo, demo_df$HEMISPHERIC)))
print(prop.table(table(demo_df$Sexo, demo_df$HEMISPHERIC), margin = 2))
cat("\nLaterality Quotient:\n")
print(wilcox.test(Laterality_Quotient ~ HEMISPHERIC, data = demo_df))
print(tapply(demo_df$Laterality_Quotient, demo_df$HEMISPHERIC, median))

cat("\n--- AFC score ---\n")
print(wilcox.test(df_full$afc_score ~ df_full$HEMISPHERIC))
print(kruskal.test(df_full$afc_score ~ df_full$HEMISPHERIC_2))
print(pairwise.wilcox.test(df_full$afc_score, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))
print(tapply(df_full$afc_score, df_full$HEMISPHERIC_2,
             function(x) c(mean = round(mean(x), 3), sd = round(sd(x), 3))))

cat("\n--- IrrLI ---\n")
print(cor.test(df_full$IrrLI, df_full$LI_GEN, method = "spearman"))
print(kruskal.test(df_full$IrrLI ~ df_full$HEMISPHERIC_2))
print(pairwise.wilcox.test(df_full$IrrLI, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))

# ==============================================================================
# 16. SUBGROUP ACCURACY
# ==============================================================================
cat("\n====== CONFIG 4 ACCURACY BY SUBGROUP ======\n")
for (cls in c("Left", "Right", "Bilateral")) {
  idx        <- which(df_full$HEMISPHERIC_2 == cls)
  true_label <- ifelse(cls == "Left", 1, 0)
  correct    <- sum(y_class_s4[idx] == true_label)
  total      <- length(idx)
  cat(sprintf("%-12s: %d / %d correctly predicted (%.1f%%)\n",
              cls, correct, total, (correct / total) * 100))
}

cat("\n====== CONFIG 9 ACCURACY BY SUBGROUP ======\n")
coords_s9  <- coords(roc_s9, "best", best.method = "closest.topleft",
                     ret = "threshold")
y_class_s9 <- ifelse(preds_s9 >= coords_s9$threshold, 1, 0)
for (cls in c("Left", "Right", "Bilateral")) {
  idx        <- which(df_full$HEMISPHERIC_2 == cls)
  true_label <- ifelse(cls == "Left", 1, 0)
  correct    <- sum(y_class_s9[idx] == true_label)
  total      <- length(idx)
  cat(sprintf("%-12s: %d / %d correctly predicted (%.1f%%)\n",
              cls, correct, total, (correct / total) * 100))
}

# ==============================================================================
# 17. INTERPRETABILITY: NODEAI COEFFICIENTS
# ==============================================================================
cat("\n====== NODEAI REGION COEFFICIENTS ======\n")
tr_means <- colMeans(X_nodeai_lq)
tr_sds   <- apply(X_nodeai_lq, 2, sd); tr_sds[tr_sds == 0] <- 1
X_scaled <- scale(X_nodeai_lq, center = tr_means, scale = tr_sds)
fit_final <- cv.glmnet(X_scaled, y, alpha = 0, family = "binomial", nfolds = 10)
coefs     <- coef(fit_final, s = "lambda.min")
region_lookup <- read.csv(paste0(BASE_PATH, "region_lookup_100.csv"))
coef_df <- data.frame(Feature = rownames(coefs)[-1],
                      Coef    = as.numeric(coefs[-1]))
coef_df <- coef_df[!coef_df$Feature %in% c("Laterality_Quotient", "Edad"), ]
feature_nums         <- as.integer(gsub("nodeAI_R", "", coef_df$Feature))
positions            <- match(feature_nums, roi_idx_simple)
coef_df$Abbreviation <- region_lookup$Abbreviation[positions]
coef_df$Full.Name    <- region_lookup$Full.Name[positions]
coef_df$Atlas        <- region_lookup$Atlas[positions]
coef_df$Hemisphere   <- region_lookup$Hemisphere[positions]
coef_df <- coef_df[order(abs(coef_df$Coef), decreasing = TRUE), ]
print(head(coef_df[, c("Abbreviation", "Full.Name", "Atlas",
                       "Hemisphere", "Coef")], 20))

# ==============================================================================
# 18. SUMMARY TABLES
# ==============================================================================
cat("\n====== SUMMARY TABLE: ALL CONFIGS ======\n")
get_metrics <- function(roc_obj, preds, y_vec) {
  coords_obj <- coords(roc_obj, "best", best.method = "closest.topleft",
                       ret = c("threshold", "sensitivity", "specificity"))
  y_class <- ifelse(preds >= coords_obj$threshold, 1, 0)
  cm      <- confusionMatrix(factor(y_class), factor(y_vec), positive = "1")
  list(AUC          = round(as.numeric(auc(roc_obj)), 3),
       Balanced_Acc = round(as.numeric(cm$byClass["Balanced Accuracy"]), 3),
       Sensitivity  = round(as.numeric(cm$byClass["Sensitivity"]), 3),
       Specificity  = round(as.numeric(cm$byClass["Specificity"]), 3),
       Kappa        = round(as.numeric(cm$overall["Kappa"]), 3))
}

roc_list  <- list(roc_s1, roc_s2, roc_s3, roc_s4, roc_s5, roc_s6, roc_s7,
                  roc_s8, roc_s9, roc_s10, roc_s11, roc_s12)
pred_list <- list(preds_s1, preds_s2, preds_s3, preds_s4, preds_s5,
                  preds_s6, preds_s7, preds_s8, preds_s9, preds_s10,
                  preds_s11, preds_s12)
config_labels <- c(
  "1: Full AFC + Ridge",
  "2: AFC PCA + Ridge",
  "3: AFC PCA + LQ + Ridge",
  "4: AFC PCA + LQ + Age + Ridge (PRIMARY)",
  "5: Elastic Net PCA + LQ + Age",
  "6: Random Forest PCA + LQ + Age",
  "7: AFC PCA + LQ + Age + IrrLI + Ridge",
  "8: NodeAI + Ridge",
  "9: NodeAI + LQ + Age + Ridge",
  "10: NodeAI + Homotopic + LQ + Age + Ridge",
  "11: Graph Features + LQ + Age + Ridge",
  "12: All Interpretable + LQ + Age + Ridge")

metrics_list <- mapply(get_metrics, roc_list, pred_list,
                       MoreArgs = list(y_vec = y), SIMPLIFY = FALSE)
results_table <- data.frame(
  Config       = config_labels,
  AUC          = sapply(metrics_list, `[[`, "AUC"),
  Balanced_Acc = sapply(metrics_list, `[[`, "Balanced_Acc"),
  Sensitivity  = sapply(metrics_list, `[[`, "Sensitivity"),
  Specificity  = sapply(metrics_list, `[[`, "Specificity"),
  Kappa        = sapply(metrics_list, `[[`, "Kappa"))
print(results_table)

cat(sprintf("\nNested CV AUC:            %.3f\n", auc(roc_nested)))
cat(sprintf("Plain LR (Config 4):      %.3f\n", auc(roc_glm)))
cat(sprintf("Cross-cohort EZM2->EZM3:  %.3f\n", auc(roc_cross1)))
cat(sprintf("Cross-cohort EZM3->EZM2:  %.3f\n", auc(roc_cross2)))
cat(sprintf("Bilateral removed:        %.3f\n", auc(roc_nobi)))





# ==============================================================================
# 19. PCA BACK-PROJECTION: CONNECTION-LEVEL IMPORTANCE
# ==============================================================================
cat("\n====== PCA BACK-PROJECTION ======\n")

# Recompute PCA on full AFC matrix (same data used for scores)
pca_obj <- prcomp(X_afc, center = TRUE, scale. = TRUE)

# Fit ridge on full dataset using Config 4 features to extract coefficients
tr_means_bp <- colMeans(X_pca_lq_age)
tr_sds_bp   <- apply(X_pca_lq_age, 2, sd); tr_sds_bp[tr_sds_bp == 0] <- 1
X_scaled_bp <- scale(X_pca_lq_age, center = tr_means_bp, scale = tr_sds_bp)
fit_full    <- cv.glmnet(X_scaled_bp, y, alpha = 0,
                         family = "binomial", nfolds = 10)

# Extract only the PC1-PC10 coefficients (exclude LQ, Age intercept)
all_coefs   <- as.numeric(coef(fit_full, s = "lambda.min"))
# coef order: intercept, PC1..PC10, Laterality_Quotient, Edad
ridge_coefs <- all_coefs[2:11]  # PC1 through PC10

# Project back onto connection space
connection_weights <- pca_obj$rotation[, 1:10] %*% ridge_coefs

weight_df <- data.frame(
  connection = afc_cols,
  weight     = as.numeric(connection_weights)
)
weight_df <- weight_df[order(abs(weight_df$weight), decreasing = TRUE), ]

cat("Top 20 connections by absolute weight:\n")
print(head(weight_df, 20))

write.csv(weight_df, paste0(BASE_PATH, "pca_backprojection_weights.csv"),
          row.names = FALSE)

# ==============================================================================
# 20. PREDICTED PROBABILITY STRIP PLOT
# ==============================================================================
cat("\n====== PREDICTED PROBABILITY STRIP PLOT ======\n")

plot_strip <- data.frame(
  prob       = preds_s4,
  true_group = factor(df_full$HEMISPHERIC_2,
                      levels = c("Left", "Right", "Bilateral")),
  correct    = factor(
    ifelse(y_class_s4 == y, "Correct", "Misclassified"),
    levels = c("Correct", "Misclassified")
  )
)

group_colors <- c("Left" = "#D64E4E", "Right" = "#2B81AC", "Bilateral" = "#7A51A0")

ggplot(plot_strip, aes(x = prob, y = true_group, color = true_group)) +
  geom_jitter(
    aes(shape = correct, size = correct, alpha = correct),
    height = 0.25, width = 0, seed = 42
  ) +
  geom_vline(
    xintercept = coords_s4$threshold,
    linetype = "dashed", color = "grey30", linewidth = 0.7
  ) +
  annotate(
    "text",
    x     = coords_s4$threshold + 0.02,
    y     = 3.45,
    label = sprintf("Threshold = %.2f", coords_s4$threshold),
    hjust = 0, size = 3.2, color = "grey30"
  ) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = c("Correct" = 16, "Misclassified" = 4)) +
  scale_size_manual(values  = c("Correct" = 1.8, "Misclassified" = 2.5)) +
  scale_alpha_manual(values = c("Correct" = 0.55, "Misclassified" = 0.95)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x     = "Predicted probability (typical lateralization)",
    y     = "True lateralization group",
    shape = NULL,
    size  = NULL,
    alpha = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position   = "bottom",
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )

ggsave(paste0(BASE_PATH, "figures/strip_plot_probabilities.png"),
       width = 7, height = 4, dpi = 300)

# ==============================================================================
# 20b. STRIP PLOT - BILATERAL REMOVED
# ==============================================================================

coords_nobi <- coords(roc_nobi_obj, "best", best.method = "closest.topleft",
                      ret = c("threshold", "sensitivity", "specificity"))

y_class_nobi <- ifelse(preds_nobi >= as.numeric(coords_nobi$threshold[1]), 1, 0)

plot_strip_nobi <- data.frame(
  prob       = preds_nobi,
  true_group = factor(df_nobi$HEMISPHERIC_2,
                      levels = c("Left", "Right")),
  correct    = factor(
    ifelse(y_class_nobi == y_nobi, "Correct", "Misclassified"),
    levels = c("Correct", "Misclassified")
  )
)

ggplot(plot_strip_nobi, aes(x = prob, y = true_group, color = true_group)) +
  geom_jitter(
    aes(shape = correct, size = correct, alpha = correct),
    position = position_jitter(height = 0.25, width = 0, seed = 42)
  ) +
  geom_vline(
    xintercept = as.numeric(coords_nobi$threshold[1]),
    linetype = "dashed", color = "grey30", linewidth = 0.7
  ) +
  annotate(
    "text",
    x     = as.numeric(coords_nobi$threshold[1]) + 0.02,
    y     = 2.4,
    label = sprintf("Threshold = %.2f", as.numeric(coords_nobi$threshold[1])),
    hjust = 0, size = 3.2, color = "grey30"
  ) +
  scale_color_manual(values = c("Left" = "#D64E4E", "Right" = "#2B81AC"),
                     guide = "none") +
  scale_shape_manual(values = c("Correct" = 16, "Misclassified" = 4)) +
  scale_size_manual(values  = c("Correct" = 1.8, "Misclassified" = 2.5)) +
  scale_alpha_manual(values = c("Correct" = 0.55, "Misclassified" = 0.95)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x     = "Predicted probability (typical lateralization)",
    y     = "True lateralization group",
    shape = NULL,
    size  = NULL,
    alpha = NULL,
    title = "Bilateral subjects removed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom",
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA)
  )

ggsave(paste0(BASE_PATH, "figures/strip_plot_nobi.png"),
       width = 7, height = 3.5, dpi = 300)

# ==============================================================================
# 21. ROC CURVE WITH SUBGROUP OVERLAYS
# ==============================================================================
cat("\n====== ROC CURVE WITH SUBGROUP OVERLAYS ======\n")

# Full model ROC (already have roc_s4)
# Build per-subgroup ROC curves
# Left subjects: classified as typical (1) vs atypical (0)
left_idx  <- which(df_full$HEMISPHERIC_2 == "Left")
right_idx <- which(df_full$HEMISPHERIC_2 == "Right")
bil_idx   <- which(df_full$HEMISPHERIC_2 == "Bilateral")

# For subgroup ROCs we compare each group against the rest
roc_left <- roc(
  response  = ifelse(df_full$HEMISPHERIC_2[c(left_idx, right_idx, bil_idx)] == "Left", 1, 0),
  predictor = preds_s4[c(left_idx, right_idx, bil_idx)],
  quiet     = TRUE
)

# More informative: one-vs-rest predicted probability distributions per group
# Build ROC data frames for ggplot
roc_to_df <- function(roc_obj, label) {
  data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    model = label
  )
}

roc_nobi_obj <- roc(y[no_bil_idx], preds_nobi, quiet = TRUE)

roc_df_full  <- roc_to_df(roc_s4,      sprintf("All subjects (AUC = %.3f)", auc(roc_s4)))
roc_df_nobi  <- roc_to_df(roc_nobi_obj, sprintf("Bilateral removed (AUC = %.3f)", auc(roc_nobi_obj)))

roc_combined <- rbind(roc_df_full, roc_df_nobi)
roc_combined$model <- factor(roc_combined$model,
                             levels = c(
                               sprintf("All subjects (AUC = %.3f)", auc(roc_s4)),
                               sprintf("Bilateral removed (AUC = %.3f)", auc(roc_nobi_obj))
                             ))

# Replace the op_point definition with this:
coords_s4_full <- coords(roc_s4, "best", best.method = "closest.topleft",
                         ret = c("threshold", "sensitivity", "specificity"))

op_point <- data.frame(
  FPR   = 1 - as.numeric(coords_s4_full$specificity[1]),
  TPR   = as.numeric(coords_s4_full$sensitivity[1]),
  model = sprintf("All subjects (AUC = %.3f)", auc(roc_s4))
)

label_full <- sprintf("All subjects (AUC = %.3f)", auc(roc_s4))
label_nobi <- sprintf("Bilateral removed (AUC = %.3f)", auc(roc_nobi_obj))

roc_colors <- c("#2B4590", "#D64E4E")
names(roc_colors) <- c(label_full, label_nobi)

ggplot(roc_combined, aes(x = FPR, y = TPR, color = model)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey70", linewidth = 0.6) +
  geom_line(linewidth = 1.1) +
  geom_point(
    data  = op_point,
    aes(x = FPR, y = TPR),
    color = "#2B4590", size = 3.5, shape = 16
  ) +
  annotate(
    "text",
    x     = op_point$FPR + 0.03,
    y     = op_point$TPR - 0.04,
    label = "Operating\nthreshold",
    size  = 3, color = "#2B4590", hjust = 0
  ) +
  scale_color_manual(values = roc_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x     = "False positive rate (1 - Specificity)",
    y     = "True positive rate (Sensitivity)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    legend.position   = c(0.65, 0.15),
    legend.background = element_rect(fill = "white", color = "grey85",
                                     linewidth = 0.4),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )
ggsave(paste0(BASE_PATH, "figures/roc_subgroup_overlay.png"),
       width = 6, height = 6, dpi = 300)



# ==============================================================================
# SENSITIVITY: RETRAIN ON HEMISPHERIC_3 LABELS
# ==============================================================================
cat("\n====== RETRAIN ON HEMISPHERIC_3 LABELS ======\n")

y_h3 <- ifelse(df_full$HEMISPHERIC_3 == "Left", 1, 0)
set.seed(42)
folds_h3  <- createFolds(y_h3, k = 10, returnTrain = TRUE)
preds_h3  <- run_ridge_cv(X_pca_lq_age, y_h3, folds_h3)
roc_h3    <- roc(y_h3, preds_h3, quiet = TRUE)
cat(sprintf("HEMISPHERIC_3 AUC: %.3f (95%% CI: %.3f-%.3f)\n",
            auc(roc_h3), ci(roc_h3)[1], ci(roc_h3)[3]))

# DeLong vs primary
print(roc.test(roc_h3, roc_s4, method = "delong"))

# Subgroup accuracy under HEMISPHERIC_3
coords_h3 <- coords(roc_h3, "best", best.method = "closest.topleft",
                    ret = c("threshold", "sensitivity", "specificity"))
y_class_h3 <- ifelse(preds_h3 >= as.numeric(coords_h3$threshold[1]), 1, 0)

cat("\nSubgroup accuracy (HEMISPHERIC_3 labels, retrained model):\n")
for (cls in c("Left", "Bilateral", "Right")) {
  idx        <- which(df_full$HEMISPHERIC_3 == cls)
  true_label <- ifelse(cls == "Left", 1, 0)
  correct    <- sum(y_class_h3[idx] == true_label)
  total      <- length(idx)
  cat(sprintf("%-12s: %d / %d correctly predicted (%.1f%%)\n",
              cls, correct, total, (correct / total) * 100))
}

cat(sprintf("\nDecision threshold: %.3f\n", as.numeric(coords_h3$threshold[1])))


