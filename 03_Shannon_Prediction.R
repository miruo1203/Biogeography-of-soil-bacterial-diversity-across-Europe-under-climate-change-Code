library(ranger)
library(caret)
library(dplyr)

## -----------------------------------------------------------
## 1. Function for model evaluation
## -----------------------------------------------------------
fun.acc <- function(observed, predicted) {
  # Coefficient of determination (explained variance)
  R2 <- 1 - sum((observed - predicted)^2) /
    sum((observed - mean(observed))^2)
  
  # Standard error of prediction (SEP)
  SEP2 <- mean((observed - predicted)^2)
  SEP  <- sqrt(SEP2)
  
  # Bias
  bias <- mean(predicted) - mean(observed)
  
  # Residual variance (bias-corrected SEP)
  SEP2c <- sum((predicted - bias - observed)^2) / length(observed)
  SEPc  <- sqrt(SEP2c)
  
  # Ratio of performance to deviation (RPD)
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance (RPIQ)
  q    <- quantile(observed)
  IQ   <- q[4] - q[2]
  RPIQ <- IQ / SEP
  
  # Concordance correlation coefficient (CCC)
  mx  <- mean(observed)
  my  <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed - mx) * (predicted - my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  data.frame(
    R2          = R2,
    concordance = ccc,
    RMSE        = SEP,
    bias        = bias,
    RPD         = RPD,
    RPIQ        = RPIQ,
    row.names   = NULL
  )
}

## -----------------------------------------------------------
## 2. Forward Recursive Feature Selection (FRFS) using ranger
## -----------------------------------------------------------
FRFS <- function(data, var.res, method = "rf", early.stop = FALSE) {
  # data: data.frame containing response (var.res) and predictors
  # var.res: name of the response variable (character)
  # method: currently only "rf" (random forest) is supported
  # early.stop: if TRUE, stop when a local RMSE minimum is reached
  
  require(caret)
  require(ranger)
  
  if (method != "rf") {
    stop("Error: The current version only supports random forest (method = 'rf').")
  }
  
  # Identify response and predictors
  col.var.res <- which(colnames(data) == var.res)
  var.exp     <- colnames(data)[-col.var.res]
  
  # Initial RF to get variable importance
  set.seed(666)
  rf.temp <- ranger(
    x          = data[, -col.var.res],
    y          = data[,  col.var.res],
    num.trees  = 500,
    importance = "permutation"
  )
  
  vip        <- as.data.frame(rf.temp$variable.importance)
  var.choose <- rownames(vip)[which(vip[, 1] == max(vip[, 1]))]
  var.left   <- var.exp[!(var.exp %in% var.choose)]
  
  # Determine optimal variables
  set.seed(666)
  fitControl <- trainControl(method = "cv", number = 5)  # 5-fold CV
  Grid <- expand.grid(
    mtry          = seq(1, length(var.choose), 1),
    min.node.size = 5,
    splitrule     = "variance"
  )
  formula <- as.formula(
    paste0(var.res, " ~ ", paste(var.choose, collapse = "+"))
  )
  rf.temp <- train(
    formula,
    data      = data,
    method    = "ranger",
    trControl = fitControl,
    tuneGrid  = Grid
  )
  
  RMSE.all <- data.frame(No.var = 1:length(var.exp), RMSE = NA_real_)
  RMSE.all[length(var.choose), "RMSE"] <- min(rf.temp$results$RMSE)
  
  # Forward recursive feature selection loop
  for (iter in 2:length(var.exp)) {
    RMSE.temp <- data.frame(Var = var.left, RMSE = NA_real_)
    
    for (i in seq_along(var.left)) {
      var.test <- c(var.choose, var.left[i])
      
      set.seed(666)
      fitControl <- trainControl(method = "cv", number = 5)
      Grid <- expand.grid(
        mtry          = seq(1, length(var.test), 1),
        min.node.size = 5,
        splitrule     = "variance"
      )
      formula <- as.formula(
        paste0(var.res, " ~ ", paste(var.test, collapse = "+"))
      )
      rf.temp <- train(
        formula,
        data      = data,
        method    = "ranger",
        trControl = fitControl,
        tuneGrid  = Grid
      )
      RMSE.temp[i, "RMSE"] <- min(rf.temp$results$RMSE)
    }
    
    if (min(RMSE.temp$RMSE) < min(RMSE.all$RMSE, na.rm = TRUE)) {
      # Add the best variable
      best.var   <- RMSE.temp$Var[which.min(RMSE.temp$RMSE)]
      var.choose <- c(var.choose, best.var)
      var.left   <- var.exp[!(var.exp %in% var.choose)]
      RMSE.all[length(var.choose), "RMSE"] <- min(RMSE.temp$RMSE)
    } else {
      if (isTRUE(early.stop)) {
        RMSE.all[length(var.choose) + 1, "RMSE"] <- min(RMSE.temp$RMSE)
        break
      } else {
        best.var   <- RMSE.temp$Var[which.min(RMSE.temp$RMSE)]
        var.choose <- c(var.choose, best.var)
        var.left   <- var.exp[!(var.exp %in% var.choose)]
        RMSE.all[length(var.choose), "RMSE"] <- min(RMSE.temp$RMSE)
      }
    }
  }
  
  # Final selected variables
  if (!early.stop) {
    var.final <- var.choose[1:which.min(RMSE.all$RMSE)]
  } else {
    var.final <- var.choose
  }
  
  list(var.order = var.choose, var.final = var.final, RMSE = RMSE.all)
}

## -----------------------------------------------------------
## 3. Variable importance (percentage, sorted)
## -----------------------------------------------------------
importance_percentage_sorted <- function(rf_model) {
  importance_values      <- rf_model$variable.importance
  importance_percentages <- importance_values / sum(importance_values) * 100
  
  importance_df <- data.frame(
    Variable  = names(importance_percentages),
    Importance = as.numeric(importance_percentages),
    row.names = NULL
  )
  
  importance_df %>% arrange(desc(Importance))
}

## -----------------------------------------------------------
## 4. Data import and preprocessing
## -----------------------------------------------------------
setwd("C:/Users/Administrator/Desktop")

data <- read.csv("MS_LUCAS.csv")

data$Landcover_ <- as.factor(data$Landcover_)
data$texttureUS <- as.factor(data$texttureUS)

set.seed(666)

row.cal <- sample(1:nrow(data), round(nrow(data) * 0.75))
row.cal <- sort(row.cal)

train_data      <- data[row.cal, ]
validation_data <- data[-row.cal, ]

## -----------------------------------------------------------
## 5. Full RF model for Shannon diversity (shannon_e)
## -----------------------------------------------------------
rf_full <- ranger(
  x         = data[row.cal, 21:87],
  y         = data$shannon_e[row.cal],
  num.trees = 500
)

# Save full model to disk
saveRDS(rf_full, "MS_shannon_e_full_RF_model.rds")

# Validation predictions
validation_results_full <- data.frame(
  obs  = data$shannon_e[-row.cal],
  pred = predict(rf_full, data[-row.cal, 21:87])$predictions
)

fun.acc(validation_results_full$obs, validation_results_full$pred)

write.csv(
  validation_results_full,
  file      = "MS_shannon_e_validation_results.csv",
  row.names = FALSE
)

# Training predictions
train_results_full <- data.frame(
  obs  = data$shannon_e[row.cal],
  pred = predict(rf_full, data[row.cal, 21:87])$predictions
)

write.csv(
  train_results_full,
  file      = "MS_shannon_e_train_results.csv",
  row.names = FALSE
)

importance_full <- importance_percentage_sorted(rf_full)
print(importance_full)

## -----------------------------------------------------------
## 6. FRFS-based RF model for Shannon diversity (shannon_e)
## -----------------------------------------------------------
# Here column 16 is assumed to be "shannon_e"
frfs_data <- data[, c(16, 21:87)]

frfs <- FRFS(frfs_data, var.res = "shannon_e", early.stop = TRUE)

# Selected variables
frfs$var.final
frfs$var.order

# Indices of selected variables in the original data
var.sel <- which(colnames(data) %in% frfs$var.order)

# RF with FRFS-selected predictors
rf_frfs <- ranger(
  x          = data[row.cal, var.sel],
  y          = data$shannon_e[row.cal],
  num.trees  = 500,
  importance = "permutation"
)

# Save FRFS RF model and selected features
saveRDS(rf_frfs,       "MS_shannon_e_FRFS_RF_model.rds")
saveRDS(frfs$var.final, "MS_shannon_e_FRFS_selected_features.rds")

importance_frfs <- importance_percentage_sorted(rf_frfs)
print(importance_frfs)

# Validation predictions (FRFS model)
validation_results_frfs <- data.frame(
  obs  = data$shannon_e[-row.cal],
  pred = predict(rf_frfs, data[-row.cal, var.sel])$predictions
)

fun.acc(validation_results_frfs$obs, validation_results_frfs$pred)

write.csv(
  validation_results_frfs,
  file      = "MS_shannon_e_FRFS_validation_results.csv",
  row.names = FALSE
)

# Training predictions (FRFS model)
train_results_frfs <- data.frame(
  obs  = data$shannon_e[row.cal],
  pred = predict(rf_frfs, data[row.cal, var.sel])$predictions
)

write.csv(
  train_results_frfs,
  file      = "MS_shannon_e_FRFS_train_results.csv",
  row.names = FALSE
)
