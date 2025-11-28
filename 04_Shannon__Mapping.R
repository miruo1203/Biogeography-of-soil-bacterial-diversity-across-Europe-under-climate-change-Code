# ============================================
# Libraries
# ============================================
library(ranger)
library(terra)
library(dplyr)

# ============================================
# 1. Load pre-trained FRFS Random Forest model (Shannon)
# ============================================
# Change this to the actual path of your FRFS model for shannon_e
shannon_model_path <- "C:/Users/Administrator/Desktop/Shannon_FRFS_model.rds"

rf_frfs_shannon <- readRDS(shannon_model_path)

# Extract selected feature names used in the model
selected_features_shannon <- rf_frfs_shannon$forest$independent.variable.names

cat("Loaded FRFS model for shannon_e\n")
cat("Selected features:\n")
print(selected_features_shannon)

# ============================================
# 2. Load environmental covariate rasters (no tiling)
# ============================================
base_dir <- "D:/Bacterial_mapping/Shannon_initial"

env_folders <- c(
  "C_N",
  "ChannelNetworkBase",
  "Clay",
  "pH",
  "Sand",
  "Silt",
  "TotalNumber_perSquareMeter"
)

cov_list <- list()

for (folder in env_folders) {
  folder_path <- file.path(base_dir, folder)
  tif_files   <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) == 0) {
    message("No raster file found in: ", folder_path)
    next
  }
  
  # Assume a single full-extent raster per folder
  r <- rast(tif_files[1])
  names(r) <- folder
  
  cov_list[[folder]] <- r
  message("Loaded raster: ", folder, " from ", tif_files[1])
}

# Stack all covariates
cov_stack_shannon <- do.call(c, cov_list)

cat("Covariate raster names:\n")
print(names(cov_stack_shannon))

# ============================================
# 3. Match raster layers to model features
# ============================================
feature_names_shannon <- intersect(
  selected_features_shannon,
  names(cov_stack_shannon)
)

if (length(feature_names_shannon) == 0) {
  stop("No matching features between covariates and the shannon_e model.")
}

cat("Features used for prediction:\n")
print(feature_names_shannon)

cov_stack_sub_shannon <- cov_stack_shannon[[feature_names_shannon]]

# Convert to data.frame for prediction
cov_df_shannon <- as.data.frame(cov_stack_sub_shannon, na.rm = FALSE)

# Identify pixels with complete covariate data
valid_idx <- complete.cases(cov_df_shannon)

cat("Total pixels: ", nrow(cov_df_shannon), "\n")
cat("Valid pixels (no NA in features): ", sum(valid_idx), "\n")

# ============================================
# 4. Predict Shannon diversity over full extent
# ============================================
pred_values_shannon <- rep(NA_real_, nrow(cov_df_shannon))

if (sum(valid_idx) > 0) {
  pred_values_shannon[valid_idx] <- predict(
    rf_frfs_shannon,
    data = cov_df_shannon[valid_idx, , drop = FALSE]
  )$predictions
} else {
  warning("No valid pixels for prediction (all rows contain NA).")
}

# ============================================
# 5. Convert predictions to raster and save
# ============================================
template_raster_shannon <- cov_stack_sub_shannon[[1]]

pred_raster_shannon <- rast(template_raster_shannon)
values(pred_raster_shannon) <- pred_values_shannon
names(pred_raster_shannon) <- "Shannon_e_FRFS"

output_dir <- file.path(base_dir, "shannon_FRFS_map")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_path <- file.path(output_dir, "Shannon_e_FRFS_predicted.tif")

writeRaster(pred_raster_shannon, output_path, overwrite = TRUE)
plot(pred_raster_shannon, main = "Shannon_e predicted using FRFS RF")

