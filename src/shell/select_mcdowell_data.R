library(dplyr)
library(maptools)
library(maps)
library(rworldmap)
library(countrycode)
library(ggplot2)

#### Loading the data
# Mcdowell 2023
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/OlsenP_2023/"
mcdowell  <- read.csv(paste0(dir,"Final_filtered_data.csv"),  header = TRUE) %>% dplyr::filter(!is.na(OBJECTID))
colnames(mcdowell)[8] <- "lat"
colnames(mcdowell)[9] <- "lon"
# mcdowell_latlon  <- read.csv(paste0(dir,"Final_filtered_data.csv"),  header = TRUE) %>%
#   dplyr::select(Number:OlsenP, Pred_Corrected, lat = "Lat", lon = "Long","Continent" ) %>%
#   dplyr::filter(!is.na(OBJECTID))

# Continents
table(mcdowell$Continent)

# P content
mcdowell$OlsenP


# Subsetting Africa
mcdowell_africa <- mcdowell[which(mcdowell$Continent == "Africa"), ]
#mcdowell_africa <- mcdowell_africa[,-c(1,2,3,"BIOME_NAME","UID","Original.Database":"Depth","X","Countryf":"Continent","Predicted":"Residualsa")]
#mcdowell_africa <- subset(mcdowell_africa, select = c(OlsenP, EVI_01:SubSoilGro, Biom, Depth.end,SoilSGF:Croplandf,Scode:allpca17))

mcdowell_africa <- mcdowell_africa[,-c(1:3.108,110:116,118,123,125,143:147)]
mcdowell_africa <- subset(mcdowell_africa, select = -c(Country,Scode))
mcdowell_africa$OlsenP <- log(mcdowell_africa$OlsenP)

str(mcdowell_africa)


# Split the data into training and testing sets
set.seed(123)
train_idx <- sample(nrow(mcdowell_africa), size = floor(0.8 * nrow(mcdowell_africa)), replace = FALSE)
train_data <- mcdowell_africa[train_idx, ]
test_data <- mcdowell_africa[-train_idx, ]

# Select the variables to use as predictors
predictors <- c("EVI_01", "EVI_02", "EVI_03", "lat", "lon", "runoff1", "prec_01", "prec_02", "prec_03", "prec_04",
                "prec_05", "prec_06", "prec_07", "prec_08", "prec_09", "prec_10", "prec_11", "prec_12", "tavg_01",
                "tavg_02", "tavg_03", "tavg_04", "tavg_05", "tavg_06", "tavg_07", "tavg_08", "tavg_09", "tavg_10",
                "tavg_11", "tavg_12", "tmax_01", "tmax_02", "tmax_03", "tmax_04", "tmax_05", "tmax_06", "tmax_07",
                "tmax_08", "tmax_09", "tmax_10", "tmax_11", "tmax_12", "tmin_01", "tmin_02", "tmin_03", "tmin_04",
                "tmin_05", "tmin_06", "tmin_07", "tmin_08", "tmin_09", "tmin_10", "tmin_11", "tmin_12", "PC_CrownCl",
                "bio01", "bio02", "bio03", "bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11",
                "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")

# Create the training and testing marker and phenotype matrices
train_pheno <- train_data$OlsenP
train_marker <- as.matrix(train_data[, predictors])
test_marker <- as.matrix(test_data[, predictors])

# Use the RNN_func to predict OlsenP
prediction <- RNN_func(train_pheno, train_marker, test_marker)

# Calculate Permutation Importance
importance_scores <- calculate_permutation_importance(prediction$model, train_marker_reshaped, train_pheno, mae)

# Print the importance scores for each feature
importance_df <- data.frame(Feature = colnames(train_marker), Importance = importance_scores)
importance_df <- importance_df[order(-importance_df$Importance), ]
print(importance_df)


# Select the top features based on the importance scores
all_features <- c("prec_10", "bio19", "runoff1", "bio18", "prec_06", "prec_12", "prec_01", "bio11", "tavg_02", "bio08",
                  "prec_03", "tavg_06", "tmin_08", "bio07", "tmax_06", "tavg_04", "tavg_08", "tavg_09", "tavg_11",
                  "tmax_08", "tmax_11", "tmin_03", "tmin_06", "tmin_07", "tmin_12", "bio01", "bio02", "bio10",
                  "tavg_05", "tavg_07", "tavg_10", "tavg_12", "tmax_01", "tmax_07", "tmin_01", "tmin_10", "bio03",
                  "bio05", "bio06", "tmax_12", "tmin_09", "bio09", "tmax_04", "tmax_09", "tmin_02", "tavg_01",
                  "tmin_05", "tmin_11", "tmax_05", "bio14", "tmax_03", "tmax_10", "bio15", "tmax_02", "tmin_04",
                  "tavg_03", "prec_02", "PC_CrownCl", "prec_08", "bio17", "prec_09", "prec_05", "bio13", "lat",
                  "prec_11", "bio04", "prec_04", "bio12", "prec_07", "bio16", "EVI_03", "EVI_02", "EVI_01", "lon")

# Add the response variable to the list of selected features
selected_features <- c("OlsenP", all_features)

# Extract the selected features from the mcdowell_africa dataframe
mcdowell_africa_selected <- subset(mcdowell_africa, select = selected_features)


# Random forrest
# Split the data into training and testing sets
set.seed(42)  # for reproducibility
train_idx <- sample(1:nrow(mcdowell_africa), 0.8 * nrow(mcdowell_africa))  # 80% for training
train_data <- mcdowell_africa[train_idx, ]
test_data <- mcdowell_africa[-train_idx, ]

# Train the Random Forest model
rf_model <- randomForest(OlsenP ~ ., data = train_data, ntree = 500, mtry = floor(sqrt(ncol(train_data) - 1)), importance = TRUE)


# Predict on the test data
predictions <- predict(rf_model, newdata = test_data)

# Predict on the training data
train_predictions <- predict(rf_model, newdata = train_data)
# Calculate the correlation for training set
train_correlation <- cor(train_data$OlsenP, train_predictions)
cat("Training set correlation:", train_correlation, "\n")

# Calculate the correlation for testing set
test_correlation <- cor(test_data$OlsenP, predictions)
cat("Testing set correlation:", test_correlation)



# Evaluate the performance using mean squared error (MSE)
mse <- mean((test_data$OlsenP - predictions)^2)
cat("Mean Squared Error (MSE):", mse)


# View the predicted values
head(prediction$val_predicted)

test_predicted <- prediction$val_predicted
train_predicted <- prediction$train_predicted

cor(train_pheno,train_predicted)
cor(test_data$OlsenP, test_predicted)








library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(caret)
library(randomForest)
library(rrBLUP)
library(tensorflow)
library(keras)
library(kerasR)
library(doParallel)
library(foreach)
library(vroom)
library(dplyr)
library(reticulate)
library(e1071)
library(glmnet)
library(gbm)
library(xgboost)

RNN_func <- function(train_pheno, train_marker, test_marker) {
  # Define the RNN model
  model <- keras_model_sequential()
  model %>% 
    layer_lstm(units = 64, input_shape = c(1, ncol(train_marker))) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_dense(units = 1, activation = "linear")
  
  # Compile the model
  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adam(learning_rate = 0.001)
  )
  
  # Reshape the input data
  train_marker <- array_reshape(train_marker, c(nrow(train_marker), 1, ncol(train_marker)))
  test_marker <- array_reshape(test_marker, c(nrow(test_marker), 1, ncol(test_marker)))
  
  # Fit the model to the training data
  history <- model %>% fit(
    x = train_marker,
    y = train_pheno,
    epochs = 50,
    batch_size = 32,
    validation_split = 0.2,
    verbose = 2
  )
  
  # Make predictions on the test and train data
  val_predicted <- model %>% predict(test_marker)
  train_predicted <- model %>% predict(train_marker)
  return_value = list("val_predicted"=val_predicted, "train_predicted"=train_predicted, "model"=model)
  # Return the test predictions
  return(return_value)
}

calculate_permutation_importance <- function(model, X, y, metric) {
  original_score <- metric(y, model %>% predict(X))
  scores <- numeric(dim(X)[3])
  
  for (i in 1:dim(X)[3]) {
    X_permuted <- X
    X_permuted[, 1, i] <- sample(X_permuted[, 1, i])
    permuted_score <- metric(y, model %>% predict(X_permuted))
    scores[i] <- original_score - permuted_score
  }
  
  return(scores)
}


# Calculate Permutation Importance
importance_scores <- calculate_permutation_importance(prediction$model, train_marker_reshaped, train_pheno, mae)

# Print the importance scores for each feature
importance_df <- data.frame(Feature = colnames(train_marker), Importance = importance_scores)
importance_df <- importance_df[order(-importance_df$Importance), ]
print(importance_df)




DNN_func = function (train_pheno, train_geno, test_geno) {
  train_geno_mnum = train_geno + 1
  train_geno_mnum = tf$cast(train_geno_mnum, tf$float32) / 3
  test_geno_mnum = test_geno + 1
  test_geno_mnum <- tf$cast(test_geno_mnum, tf$float32) /3
  batchs = round(dim(train_geno_mnum)[1]/20)
  if ((batchs %% 2) == 1){batchs = batchs + 1}
  model <- keras_model_sequential()
  model <- model %>%
    layer_dense(units = 256, activation = "relu", input_shape = dim(train_geno_mnum)[2], kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.3, l2 = 0.3)) %>%
    layer_dropout(rate = 0.03) %>%
    layer_dense(units = 128, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.2, l2 = 0.2)) %>%
    layer_dropout(rate = 0.02) %>%
    layer_dense(units = 64, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dropout(rate = 0.01) %>%
    layer_dense(units = 32, activation = "relu", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_batch_normalization(batch_size = batchs) %>%
    layer_dense(units = 16, activation = "linear", kernel_initializer = 'orthogonal', kernel_regularizer=regularizer_l1_l2(l1 = 0.0, l2 = 0)) %>%
    layer_dense(units = 1, activation = "linear")
  model <- model %>% compile(loss = "mse", optimizer = tf$keras$optimizers$legacy$Adamax(learning_rate=0.001, decay = 0.0003), metrics = c("mean_absolute_error", "mean_squared_error"))
  #model %>% compile(loss = "mse", optimizer = optimizer_adamax(lr=0.001, decay = 0.0003), metrics = c(metric_r2_score, metric_cor))
  tensorflow::tf$random$set_seed(200)
  tf$compat$v1$set_random_seed(200)
  set.seed(200)
  tensorflow::tf$random$set_seed(200)
  history <- model %>%  
    fit(train_geno_mnum, train_pheno, epochs = 250, batch_size = 10, validation_split = 0.2, verbose = 0,
        callbacks = list(callback_early_stopping(patience = 30), callback_reduce_lr_on_plateau(factor = 0.1)))
  train_pred <- model %>% predict(train_geno_mnum, batch_size = batchs)
  val_pred <- model %>% predict(test_geno_mnum, batch_size = batchs)
  return_value = list("val_predicted"=val_pred, "train_predicted"=train_pred, "model"=model, "batch"=batchs)
  return(return_value)
}



DNN_func <- function(train_data, test_data) {
  # Preprocess data
  train_x <- train_data[, 1:(ncol(train_data) - 1)]
  train_y <- train_data[, ncol(train_data)]
  test_x <- test_data[, 1:(ncol(test_data) - 1)]
  test_y <- test_data[, ncol(test_data)]
  
  # Scaling the data
  train_x <- scale(train_x)
  test_x <- scale(test_x)
  
  # Build the model
  model <- keras_model_sequential() %>%
    layer_dense(units = 256, activation = 'relu', input_shape = ncol(train_x)) %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = 64, activation = 'relu') %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = 32, activation = 'relu') %>%
    layer_batch_normalization() %>%
    layer_dense(units = 16, activation = 'relu') %>%
    layer_dense(units = 1, activation = 'sigmoid')
  
  # Compile the model
  model %>% compile(
    loss = 'mean_squared_error',
    optimizer = optimizer_adam(learning_rate = 0.001 + 1e-7),
    metrics = list('mean_absolute_error')
  )
  
  # Train the model
  history <- model %>% fit(
    train_x, train_y,
    epochs = 100,
    batch_size = 32,
    validation_split = 0.2,
    callbacks = list(callback_early_stopping(patience = 10))
  )
  
  # Evaluate the model
  train_predicted <- model %>% predict(train_x)
  test_predicted <- model %>% predict(test_x)
  
  cat("Training Mean Absolute Error: ", mean(abs(train_predicted - train_y)), "\n")
  cat("Testing Mean Absolute Error: ", mean(abs(test_predicted - test_y)), "\n")
  
  return(list(train_predicted = train_predicted, test_predicted = test_predicted))
}


