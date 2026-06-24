library(pROC)

# ==============================================================================
# SETUP
# ==============================================================================
BASE_PATH <- "C:/UPC/Thesis/Rcode/"

demo    <- read.csv(paste0(BASE_PATH, "demographics_cleaned_filtered_259.csv"))
irr_df  <- read.csv(paste0(BASE_PATH, "irrli_259.csv"))
irr_wb  <- read.csv(paste0(BASE_PATH, "INSIDEOUT_wholebrain_tau1.csv"))

# Standardize IDs
names(irr_df)[names(irr_df) == "subject_id"] <- "ID"
names(irr_wb)[names(irr_wb) == "SubjectID"]  <- "ID"
irr_wb$ID <- sub("^sub-", "", irr_wb$ID)

# Build minimal df_full with just what we need
df_full <- merge(
  demo[, c("ID", "HEMISPHERIC", "HEMISPHERIC_2", "LI_GEN",
           "Laterality_Quotient")],
  irr_df[, c("ID", "IrrLI")],
  by = "ID"
)
df_full <- merge(df_full, irr_wb[, c("ID", "IrrWhole")], 
                 by = "ID", all.x = TRUE)
df_full <- df_full[!is.na(df_full$HEMISPHERIC) &
                     !is.na(df_full$Laterality_Quotient), ]

cat("n =", nrow(df_full), "\n")
print(table(df_full$HEMISPHERIC_2))

# Check for missing IrrWhole
cat("Missing IrrWhole:", sum(is.na(df_full$IrrWhole)), "\n")

# ==============================================================================
# WHOLE-BRAIN IRREVERSIBILITY ANALYSES
# ==============================================================================
cat("\nWhole-brain irreversibility by group:\n")
print(tapply(df_full$IrrWhole, df_full$HEMISPHERIC_2,
             function(x) c(mean = round(mean(x, na.rm = TRUE), 4),
                           sd   = round(sd(x,   na.rm = TRUE), 4))))

cat("\nBinary (Typical vs Atypical):\n")
print(wilcox.test(df_full$IrrWhole ~ df_full$HEMISPHERIC))
print(tapply(df_full$IrrWhole, df_full$HEMISPHERIC,
             function(x) c(mean = round(mean(x, na.rm = TRUE), 4),
                           sd   = round(sd(x,   na.rm = TRUE), 4))))

cat("\nKruskal-Wallis (3 groups):\n")
print(kruskal.test(df_full$IrrWhole ~ df_full$HEMISPHERIC_2))

cat("\nPairwise Wilcoxon (Bonferroni):\n")
print(pairwise.wilcox.test(df_full$IrrWhole, df_full$HEMISPHERIC_2,
                           p.adjust.method = "bonferroni"))

cat("\nSpearman correlation with LI_GEN:\n")
print(cor.test(df_full$IrrWhole, df_full$LI_GEN, 
               method = "spearman"))