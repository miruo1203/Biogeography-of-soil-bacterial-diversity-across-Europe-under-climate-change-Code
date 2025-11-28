# ============================================
# Required Libraries
# ============================================
library(ranger)
library(terra)
library(dplyr)

# ============================================
# Load pre-trained FRFS Random Forest model
# ============================================
model_path <- "C:/Users/Administrator/Desktop/MS_Chao1_FRFS_RF_model.rds"
rf_frfs <- readRDS(model_path)

# Extract selected feature names directly from the model
selected_features <- rf_frfs$forest$independent.variable.names
cat("Loaded FRFS model\n")
cat("Selected features:\n")
print(selected_features)

# ============================================
# Load environmental covariate rasters
# ============================================
# Root directory containing sub-folders of covariates
base_dir <- "D:/Bacterial_mapping/Chao1_initial"

env_folders <- c(
  "C_N",
  "ChannelNetworkBase",
  "Clay",
  "pH",
  "TotalCatchment",
  "CEC",
  "Coarse_fragments",
  "LUCAS_LDCV",
  "MAT",
  "Silt",
  "TotalNumber_perSquareMeter",
  "TexttureUSDA",
  "ValleyDepth"
)

cov_list <- list()

for (folder in env_folders) {
  folder_path <- file.path(base_dir, folder)
  tif_files   <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) == 0) {
    message("No raster file found in: ", folder_path)
    next
  }
  
  r <- rast(tif_files[1])  # assume one full-extent raster per folder
  names(r) <- folder       # layer name must match trained feature name
  cov_list[[folder]] <- r
  message("Loaded raster: ", folder)
}

cov_stack <- do.call(c, cov_list)
cat("Raster layer names:\n")
print(names(cov_stack))

# ============================================
# Match raster layers to FRFS model features
# ============================================
feature_names <- intersect(selected_features, names(cov_stack))

if (length(feature_names) == 0) {
  stop("No matching features between covariates and trained model.")
}

cat("Features used for prediction:\n")
print(feature_names)

cov_stack_sub <- cov_stack[[feature_names]]

# Convert raster pixels to dataframe for prediction
cov_df <- as.data.frame(cov_stack_sub, na.rm = FALSE)

# Pixels with complete data
valid_idx <- complete.cases(cov_df)
cat("Total pixels:", nrow(cov_df), "\n")
cat("Valid pixels:", sum(valid_idx), "\n")

# ============================================
# Predict Chao1 using FRFS model
# ============================================
pred_values <- rep(NA_real_, nrow(cov_df))

if (sum(valid_idx) > 0) {
  pred_values[valid_idx] <- predict(
    rf_frfs,
    data = cov_df[valid_idx, , drop = FALSE]
  )$predictions
}

# ============================================
# Convert predictions back to raster and save
# ============================================
template_raster <- cov_stack_sub[[1]]
pred_raster <- rast(template_raster)
values(pred_raster) <- pred_values
names(pred_raster) <- "Chao1_FRFS"

output_dir <- file.path(base_dir, "Chao1_FRFS_Map")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_path <- file.path(output_dir, "Chao1_FRFS_predicted.tif")
writeRaster(pred_raster, output_path, overwrite = TRUE)



