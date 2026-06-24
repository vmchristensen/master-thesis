# ============================================================
# CREATE TYPICAL AND ATYPICAL NODEAI ASYMMETRY MAPS
# ============================================================

library(RNifti)
library(dplyr)

# ============================================================
# 1. FILE PATHS
# ============================================================

left_atlas_file  <- "atlas/left_composite_atlas.nii"
right_atlas_file <- "atlas/left_flipped_composite_atlas.nii"

region_lookup_file <- "region_lookup_100.csv"

typical_output_file  <- "nodeai_typical_asymmetry_map.nii.gz"
atypical_output_file <- "nodeai_atypical_asymmetry_map.nii.gz"

# ============================================================
# 2. LOAD DATA
# ============================================================

if(!exists("coef_df")){
  stop("coef_df not found in workspace")
}

region_lookup <- read.csv(
  region_lookup_file,
  check.names = FALSE
)

cat("Loaded", nrow(coef_df), "NodeAI coefficients\n")

# ============================================================
# 3. MERGE COEFFICIENTS WITH ATLAS LABELS
# ============================================================

nodeai_weights <- coef_df %>%
  left_join(
    region_lookup %>%
      select(
        Abbreviation,
        `Atlas Label`
      ),
    by = "Abbreviation"
  )

cat("\nMatched labels:\n")
print(sum(!is.na(nodeai_weights$`Atlas Label`)))
cat("of", nrow(nodeai_weights), "\n")

# ============================================================
# 4. LOAD ATLASES
# ============================================================

left_atlas  <- readNifti(left_atlas_file)
right_atlas <- readNifti(right_atlas_file)

cat("\nLeft atlas dimensions:\n")
print(dim(left_atlas))

cat("\nRight atlas dimensions:\n")
print(dim(right_atlas))

# ============================================================
# 5. INITIALIZE MAPS
# ============================================================

typical_map <- array(
  0,
  dim = dim(left_atlas)
)

atypical_map <- array(
  0,
  dim = dim(left_atlas)
)

# ============================================================
# 6. ASSIGN COEFFICIENTS
# ============================================================

for(i in seq_len(nrow(nodeai_weights))) {
  
  label <- nodeai_weights$`Atlas Label`[i]
  beta  <- nodeai_weights$Coef[i]
  
  if(is.na(label)) next
  
  # ----------------------------------------------------------
  # TYPICAL MAP
  #
  # Red = hemisphere favoring typical classification
  # ----------------------------------------------------------
  
  typical_map[left_atlas  == label] <-  beta
  typical_map[right_atlas == label] <- -beta
  
  # ----------------------------------------------------------
  # ATYPICAL MAP
  #
  # Red = hemisphere favoring atypical classification
  # ----------------------------------------------------------
  
  atypical_map[left_atlas  == label] <- -beta
  atypical_map[right_atlas == label] <-  beta
}

# ============================================================
# 7. QC
# ============================================================

cat("\nTypical map range:\n")
print(range(typical_map))

cat("\nAtypical map range:\n")
print(range(atypical_map))

cat("\nCoefficient range:\n")
print(range(nodeai_weights$Coef))

# ============================================================
# 8. SAVE TYPICAL MAP
# ============================================================

typical_img <- left_atlas
typical_img[] <- typical_map

writeNifti(
  typical_img,
  typical_output_file
)

# ============================================================
# 9. SAVE ATYPICAL MAP
# ============================================================

atypical_img <- left_atlas
atypical_img[] <- atypical_map

writeNifti(
  atypical_img,
  atypical_output_file
)

cat("\nSaved:\n")
cat(typical_output_file, "\n")
cat(atypical_output_file, "\n")

# ============================================================
# 10. TOP REGIONS
# ============================================================

cat("\nTop positive coefficients:\n")

nodeai_weights %>%
  arrange(desc(Coef)) %>%
  select(
    Abbreviation,
    Full.Name,
    Coef
  ) %>%
  head(10) %>%
  print()

cat("\nTop negative coefficients:\n")

nodeai_weights %>%
  arrange(Coef) %>%
  select(
    Abbreviation,
    Full.Name,
    Coef
  ) %>%
  head(10) %>%
  print()