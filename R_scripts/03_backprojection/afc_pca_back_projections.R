# ==============================================================================
# PCA BACK-PROJECTION ANALYSIS
# Goal 1: Connection-level importance map
# Goal 2: PC loadings inspection for PC3, PC4, PC6
#
# PCA computed on n=258 in R (pca_r_object.rds)
# Ridge coefficients extracted on n=258 using AFC_PCA_scores_120_R.csv
# ==============================================================================

library(glmnet)
library(dplyr)
library(ggplot2)
library(gridExtra)

BASE_PATH <- "C:/UPC/Thesis/Rcode/"

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
afc_matrix    <- read.csv(paste0(BASE_PATH, "AFC_vectors.csv"))
afc_pca       <- read.csv(paste0(BASE_PATH, "AFC_PCA_scores_120_R.csv"))
demo          <- read.csv(paste0(BASE_PATH, "demographics_cleaned_filtered_259.csv"))
region_lookup <- read.csv(paste0(BASE_PATH, "region_lookup_100.csv"))

names(afc_matrix)[names(afc_matrix) == "subject_id"] <- "ID"
names(afc_pca)[names(afc_pca)       == "subject_id"] <- "ID"

pc_cols  <- paste0("PC", 1:10)
afc_cols <- paste0("afc_", 1:9900)

# The order used when pca_r_object.rds was computed
ids_258_order <- demo$ID[!is.na(demo$Laterality_Quotient)]
cat("n in PCA order:", length(ids_258_order), "\n")

# ==============================================================================
# 2. BUILD DATASETS
# ==============================================================================
df_259 <- merge(
  demo[, c("ID", "HEMISPHERIC", "LI_GEN", "Laterality_Quotient", "Edad")],
  afc_matrix[, c("ID", afc_cols)],
  by = "ID"
)
df_259 <- df_259[!is.na(df_259$HEMISPHERIC), ]
df_258 <- df_259[!is.na(df_259$Laterality_Quotient), ]
cat("n_259:", nrow(df_259), "| n_258:", nrow(df_258), "\n")

# For ridge coefficients: use saved R PCA scores
df_pca_258 <- merge(
  demo[!is.na(demo$Laterality_Quotient),
       c("ID", "HEMISPHERIC", "Laterality_Quotient", "Edad")],
  afc_pca[, c("ID", pc_cols)],
  by = "ID"
)
df_pca_258 <- df_pca_258[!is.na(df_pca_258$HEMISPHERIC), ]
cat("n for ridge:", nrow(df_pca_258), "\n")

y_258   <- ifelse(df_pca_258$HEMISPHERIC == "Typical", 1, 0)
X_ridge <- as.matrix(df_pca_258[, c(pc_cols, "Laterality_Quotient", "Edad")])

# ==============================================================================
# 3. BUILD AFC PAIR LOOKUP TABLE
# ==============================================================================
n_rois <- 100
pairs  <- data.frame(roi_i = integer(4950), roi_j = integer(4950))
idx <- 1
for (i in 1:(n_rois - 1)) {
  for (j in (i + 1):n_rois) {
    pairs$roi_i[idx] <- i
    pairs$roi_j[idx] <- j
    idx <- idx + 1
  }
}

annotate_pairs <- function(pairs_df, offset, type_label) {
  pairs_df %>%
    mutate(
      afc_index        = row_number() + offset,
      type             = type_label,
      abbrev_i         = region_lookup$Abbreviation[roi_i],
      abbrev_j         = region_lookup$Abbreviation[roi_j],
      fullname_i       = region_lookup$`Full.Name`[roi_i],
      fullname_j       = region_lookup$`Full.Name`[roi_j],
      atlas_i          = region_lookup$Atlas[roi_i],
      atlas_j          = region_lookup$Atlas[roi_j],
      connection_label = paste0(abbrev_i, " - ", abbrev_j,
                                " (", type_label, ")")
    )
}

pair_lookup <- bind_rows(
  annotate_pairs(pairs, offset = 0,    type_label = "within"),
  annotate_pairs(pairs, offset = 4950, type_label = "between")
)
cat("Pair lookup:", nrow(pair_lookup), "connections\n")

# ==============================================================================
# 4. LOAD SAVED PCA OBJECT AND SANITY CHECK
# ==============================================================================
cat("\nLoading saved PCA object...\n")
pca_r <- readRDS(paste0(BASE_PATH, "pca_r_object.rds"))

# Align pca_r scores to df_pca_258 row order using ids_258_order
recomp_pc1 <- pca_r$x[match(df_pca_258$ID, ids_258_order), 1]
saved_pc1  <- afc_pca$PC1[match(df_pca_258$ID, afc_pca$ID)]

rho_check <- cor(saved_pc1, recomp_pc1, use = "complete.obs")
cat(sprintf("Sanity check - cor(saved PC1, aligned PC1): %.4f\n", rho_check))
cat("(Should be |1.0|; sign flip is fine)\n")

pca_obj <- pca_r

# ==============================================================================
# 5. EXTRACT RIDGE COEFFICIENTS ON 258 SUBJECTS
# ==============================================================================
cat("\n====== EXTRACTING RIDGE COEFFICIENTS (n=258) ======\n")

tr_means   <- colMeans(X_ridge)
tr_sds     <- apply(X_ridge, 2, sd)
tr_sds[tr_sds == 0] <- 1
X_ridge_sc <- scale(X_ridge, center = tr_means, scale = tr_sds)

set.seed(42)
fit_full <- cv.glmnet(X_ridge_sc, y_258, alpha = 0,
                      family = "binomial", nfolds = 10)

all_coefs      <- as.numeric(coef(fit_full, s = "lambda.min"))
ridge_pc_coefs <- all_coefs[2:11]  # PC1..PC10

cat("Ridge coefficients for PC1-PC10:\n")
for (i in 1:10) cat(sprintf("  PC%d: %+.5f\n", i, ridge_pc_coefs[i]))

# ==============================================================================
# 6. BACK-PROJECTION: CONNECTION-LEVEL IMPORTANCE
# ==============================================================================
cat("\n====== BACK-PROJECTION: CONNECTION WEIGHTS ======\n")

connection_weights <- as.numeric(pca_obj$rotation[, 1:10] %*% ridge_pc_coefs)

weight_df <- data.frame(afc_index = 1:9900, weight = connection_weights) %>%
  left_join(pair_lookup, by = "afc_index") %>%
  mutate(abs_weight = abs(weight)) %>%
  arrange(desc(abs_weight))

cat("\nTop 20 connections by |weight|:\n")
print(head(weight_df[, c("connection_label", "weight", "type")], 20))

cat("\nMean |weight| by type:\n")
print(weight_df %>%
        group_by(type) %>%
        summarise(mean_abs_weight = round(mean(abs_weight), 6),
                  n_in_top100     = sum(rank(-abs_weight) <= 100)))

write.csv(weight_df,
          paste0(BASE_PATH, "pca_backprojection_labeled.csv"),
          row.names = FALSE)
cat("Saved: pca_backprojection_labeled.csv\n")

top20 <- head(weight_df, 20) %>%
  mutate(connection_label = factor(connection_label,
                                   levels = rev(connection_label)))

ggplot(top20, aes(x = abs_weight, y = connection_label, fill = type)) +
  geom_col() +
  scale_fill_manual(values = c("within" = "#5369A0", "between" = "#D64E4E")) +
  labs(x = "Absolute back-projected weight", y = NULL,
       fill = "Connection type",
       title = "Top 20 AFC connections: classification importance") +
  theme_minimal(base_size = 11) +
  theme(plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

ggsave(paste0(BASE_PATH, "figures/backprojection_top20.png"),
       width = 9, height = 6, dpi = 300)
cat("Saved: figures/backprojection_top20.png\n")

# ==============================================================================
# 7. ROI-LEVEL IMPORTANCE
# ==============================================================================
cat("\n====== ROI-LEVEL IMPORTANCE ======\n")

roi_importance <- bind_rows(
  weight_df %>% select(roi = roi_i, abs_weight, type),
  weight_df %>% select(roi = roi_j, abs_weight, type)
) %>%
  group_by(roi) %>%
  summarise(total_importance = sum(abs_weight),
            mean_importance  = mean(abs_weight),
            n_connections    = n()) %>%
  arrange(desc(total_importance)) %>%
  left_join(region_lookup %>%
              mutate(roi = row_number()) %>%
              select(roi, Abbreviation, `Full.Name`, Atlas),
            by = "roi")

cat("\nTop 20 ROIs by total back-projected importance:\n")
print(head(roi_importance[, c("Abbreviation", "Full.Name",
                              "Atlas", "total_importance")], 20))

write.csv(roi_importance,
          paste0(BASE_PATH, "roi_importance_backprojection.csv"),
          row.names = FALSE)
cat("Saved: roi_importance_backprojection.csv\n")

# ==============================================================================
# 8. PC LOADINGS INSPECTION: PC3, PC4, PC6
# ==============================================================================
cat("\n====== PC LOADINGS INSPECTION ======\n")

target_pcs <- c(3, 4, 6)
pc_labels  <- c(
  "PC3 (rho = +0.256 with LI_GEN, p < 0.001)",
  "PC4 (rho = +0.254 with LI_GEN, p < 0.001)",
  "PC6 (rho = -0.245 with LI_GEN, p = 0.0001)"
)

loadings_results <- list()
plot_list        <- list()

for (p_idx in seq_along(target_pcs)) {
  pc_num   <- target_pcs[p_idx]
  pc_label <- pc_labels[p_idx]
  
  load_df <- data.frame(
    afc_index = 1:9900,
    loading   = pca_obj$rotation[, pc_num]
  ) %>%
    left_join(pair_lookup, by = "afc_index") %>%
    mutate(abs_loading = abs(loading)) %>%
    arrange(desc(abs_loading))
  
  loadings_results[[paste0("PC", pc_num)]] <- load_df
  
  cat(sprintf("\n--- PC%d top 15 connections ---\n", pc_num))
  print(head(load_df[, c("connection_label", "loading", "type")], 15))
  
  cat(sprintf("\nPC%d: within vs between in top 50:\n", pc_num))
  print(head(load_df, 50) %>%
          group_by(type) %>%
          summarise(n = n(),
                    mean_abs_loading = round(mean(abs_loading), 6)))
  
  write.csv(load_df,
            paste0(BASE_PATH, "pc", pc_num, "_loadings_labeled.csv"),
            row.names = FALSE)
  cat(sprintf("Saved: pc%d_loadings_labeled.csv\n", pc_num))
  
  top15 <- head(load_df, 15) %>%
    mutate(connection_label = factor(connection_label,
                                     levels = rev(connection_label)))
  
  plot_list[[p_idx]] <- ggplot(top15,
                               aes(x = loading, y = connection_label, fill = type)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c("within" = "#5369A0",
                                 "between" = "#D64E4E")) +
    labs(x = "Loading", y = NULL, fill = "Type", title = pc_label) +
    theme_minimal(base_size = 10) +
    theme(plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom")
}

g <- arrangeGrob(grobs = plot_list, ncol = 1)
ggsave(paste0(BASE_PATH, "figures/pc_loadings_pc3_pc4_pc6.png"),
       plot = g, width = 9, height = 15, dpi = 300)
cat("Saved: figures/pc_loadings_pc3_pc4_pc6.png\n")

cat("\n====== DONE ======\n")