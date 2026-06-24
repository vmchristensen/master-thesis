# PCA in R
library(dplyr)

BASE_PATH <- "C:/UPC/Thesis/Rcode/"

# Load data
afc_matrix <- read.csv(paste0(BASE_PATH, "AFC_vectors.csv"))
demo       <- read.csv(paste0(BASE_PATH, "demographics_cleaned_filtered_259.csv"))

# Standardize ID
names(afc_matrix)[names(afc_matrix) == "subject_id"] <- "ID"

# Get 258 IDs (exclude EZM30214 with missing LQ)
ids_258 <- demo$ID[!is.na(demo$Laterality_Quotient)]
cat("Subjects for PCA:", length(ids_258), "\n")

# Filter and align AFC matrix
afc_258 <- afc_matrix[afc_matrix$ID %in% ids_258, ]
afc_258 <- afc_258[match(ids_258, afc_258$ID), ]  # ensure same order as demo
cat("AFC rows:", nrow(afc_258), "\n")

# Recompute PCA
afc_cols  <- paste0("afc_", 1:9900)
X_for_pca <- as.matrix(afc_258[, afc_cols])

cat("Running prcomp on 258 x 9900 matrix...\n")
pca_r <- prcomp(X_for_pca, center = TRUE, scale. = TRUE)
cat("Done.\n")

# Variance explained
var_exp <- summary(pca_r)$importance[2, 1:10]
cat("Variance explained PC1-10:\n")
print(round(var_exp * 100, 2))
cat("Cumulative PC1-10:", round(sum(var_exp) * 100, 2), "%\n")

# Save scores (up to 120 components, same as before)
n_components <- min(120, ncol(pca_r$x))
scores_df    <- as.data.frame(pca_r$x[, 1:n_components])
scores_df$ID <- afc_258$ID

# Save scores with ID as first column
scores_df <- scores_df[, c("ID", paste0("PC", 1:n_components))]
write.csv(scores_df,
          paste0(BASE_PATH, "AFC_PCA_scores_120_R.csv"),
          row.names = FALSE)
cat("Saved: AFC_PCA_scores_120_R.csv\n")
cat("Scores dimensions:", nrow(scores_df), "x", ncol(scores_df), "\n")

# Save rotation matrix for back-projection later
saveRDS(pca_r, paste0(BASE_PATH, "pca_r_object.rds"))
cat("Saved: pca_r_object.rds (full PCA object for back-projection)\n")

# Spearman correlations with LI_GEN to check which PCs carry lateralization signal
li_gen <- demo$LI_GEN[match(afc_258$ID, demo$ID)]
cat("\nSpearman correlations PC1-10 with LI_GEN:\n")
for (k in 1:10) {
  r <- cor.test(pca_r$x[, k], li_gen, method = "spearman")
  cat(sprintf("  PC%d: rho = %+.3f, p = %.4f\n", k, r$estimate, r$p.value))
}