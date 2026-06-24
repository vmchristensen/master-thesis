# ============================================================
# CREATE 4 AFC ENDPOINT MAPS (DIVERGING WEIGHTS)
#
# FULL SAMPLE
#   Between  (raw signed weights: +w on left endpoint, -w on right endpoint)
#   Within   (raw signed weights: +w on dominant hemisphere endpoints)
#
# NO BILATERALS
#   Between
#   Within
#
# Positive values = left-dominant / drives typical prediction
# Negative values = right-dominant / drives atypical prediction
# ============================================================

library(dplyr)
library(RNifti)

# ============================================================
# SETTINGS
# ============================================================

BASE_PATH <- "C:/UPC/Thesis/Rcode/"

# ============================================================
# LOAD ATLASES
# ============================================================

left_atlas <- readNifti(
  paste0(
    BASE_PATH,
    "atlas/left_composite_atlas.nii"
  )
)

right_atlas <- readNifti(
  paste0(
    BASE_PATH,
    "atlas/left_flipped_composite_atlas.nii"
  )
)

right_atlas[is.na(right_atlas)] <- 0

# ============================================================
# LOAD LOOKUP
# ============================================================

region_lookup <- read.csv(
  paste0(
    BASE_PATH,
    "region_lookup_100.csv"
  ),
  check.names = FALSE
)

# ============================================================
# LOAD BACKPROJECTIONS
# ============================================================

weight_df <- read.csv(
  paste0(
    BASE_PATH,
    "pca_backprojection_labeled.csv"
  )
)

weight_df_nobi <- read.csv(
  paste0(
    BASE_PATH,
    "backprojection_nobil_labeled.csv"
  )
)

# ============================================================
# CHECK TYPES
# ============================================================

cat("\nFULL TYPES:\n")
print(unique(weight_df$type))

cat("\nNOBI TYPES:\n")
print(unique(weight_df_nobi$type))

# ============================================================
# FUNCTION
# ============================================================

make_endpoint_map <- function(
    df,
    connection_type,
    outfile
){
  
  cat(
    "\n=====================================\n"
  )
  
  cat(
    "Creating:",
    outfile,
    "\n"
  )
  
  tmp <- df %>%
    filter(
      type == connection_type
    )
  
  cat(
    "Connections:",
    nrow(tmp),
    "\n"
  )
  
  # ----------------------------------------------------------
  # EMPTY MAPS
  # ----------------------------------------------------------
  
  left_map <- array(
    0,
    dim = dim(left_atlas)
  )
  
  right_map <- array(
    0,
    dim = dim(right_atlas)
  )
  
  # ----------------------------------------------------------
  # PROJECT CONNECTIONS
  # ----------------------------------------------------------
  
  for(i in seq_len(nrow(tmp))){
    
    # raw signed weight (no abs())
    w <- tmp$weight[i]
    
    left_label <- region_lookup %>%
      filter(
        Abbreviation ==
          tmp$abbrev_i[i]
      ) %>%
      pull(`Atlas Label`)
    
    right_label <- region_lookup %>%
      filter(
        Abbreviation ==
          tmp$abbrev_j[i]
      ) %>%
      pull(`Atlas Label`)
    
    if(length(left_label) == 0 ||
       length(right_label) == 0){
      
      next
      
    }
    
    if(is.na(left_label) ||
       is.na(right_label)){
      
      next
      
    }
    
    left_idx  <- left_atlas  == left_label
    right_idx <- right_atlas == right_label
    
    if(connection_type == "between"){
      
      # abbrev_i is always the left endpoint
      # abbrev_j is always the right endpoint
      # +w on left endpoint, -w on right endpoint
      # so warm = left-dominant, cool = right-dominant
      
      left_map[left_idx] <-
        left_map[left_idx] + w
      
      right_map[right_idx] <-
        right_map[right_idx] - w
      
    } else {
      
      # within: weight = left-left minus right-right asymmetry
      # positive w: left hemisphere dominant -> paint left endpoints
      # negative w: right hemisphere dominant -> paint right endpoints
      
      if(w > 0){
        
        left_map[left_idx] <-
          left_map[left_idx] + w
        
        left_map[left_atlas == right_label] <-
          left_map[left_atlas == right_label] + w
        
      } else {
        
        right_map[right_idx] <-
          right_map[right_idx] + w
        
        right_map[right_atlas == left_label] <-
          right_map[right_atlas == left_label] + w
        
      }
      
    }
    
  }
  
  # ----------------------------------------------------------
  # COMBINE
  # ----------------------------------------------------------
  
  endpoint_map <-
    left_map +
    right_map
  
  # ----------------------------------------------------------
  # SAVE
  # ----------------------------------------------------------
  
  output_img    <- left_atlas
  output_img[]  <- endpoint_map
  
  writeNifti(
    output_img,
    outfile
  )
  
  # ----------------------------------------------------------
  # QC
  # ----------------------------------------------------------
  
  pos_vals <- endpoint_map[endpoint_map > 0]
  neg_vals <- endpoint_map[endpoint_map < 0]
  
  cat(
    "Positive range:",
    if(length(pos_vals) > 0) range(pos_vals) else "none",
    "\n"
  )
  
  cat(
    "Negative range:",
    if(length(neg_vals) > 0) range(neg_vals) else "none",
    "\n"
  )
  
  cat(
    "Nonzero voxels:",
    sum(endpoint_map != 0),
    "\n"
  )
  
  cat(
    "Saved:",
    outfile,
    "\n"
  )
  
}

# ============================================================
# FULL SAMPLE
# ============================================================

make_endpoint_map(
  weight_df,
  "between",
  paste0(
    BASE_PATH,
    "figures/full_between_new.nii.gz"
  )
)

make_endpoint_map(
  weight_df,
  "within",
  paste0(
    BASE_PATH,
    "figures/full_within_new.nii.gz"
  )
)

# ============================================================
# NO BILATERALS
# ============================================================

make_endpoint_map(
  weight_df_nobi,
  "between",
  paste0(
    BASE_PATH,
    "figures/nobi_between_new.nii.gz"
  )
)

make_endpoint_map(
  weight_df_nobi,
  "within",
  paste0(
    BASE_PATH,
    "figures/nobi_within_new.nii.gz"
  )
)

# ============================================================
# FINISHED
# ============================================================

cat(
  "\n=====================================\n"
)

cat(
  "ALL 4 MAPS CREATED\n"
)

cat(
  "=====================================\n"
)