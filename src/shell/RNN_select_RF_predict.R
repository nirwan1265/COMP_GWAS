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


