### Supervised Classificition Machine Learning Prototyping Auto-Pipeline
###
### Usage: For Prototyping of ML modelsand to serve as guide code for increased functionality
### Trains glmnet, xgboost, and random forest models

#Author Han 
#updated
#V0.0 - 2/15/2021 - organized to be modular
#V0.1 - 2/16/2021 - organized to be a function to be run. 
#V0.2 - 2/17/2021 - MARCC tested? not yet

#### Output
#Creates directory with saved models, performance metrics, and plots of:
# 1. Receiver operating characteristic curve (ROC)
# 2. Precision Recall Curve (PRC)
# 3. Isotonic regression calibrated vs uncalibrated plots
# 4. Individual of 1-3 and plots including blox plots of metrics themselves. 
# 5. CSV summary of all model metrics for all three algorithms with 95% Confidence Intervals. 

#### Input parameters & required files: ####

## Directory inputs (code will double check whether these directories exist and files within it are available)
#
#code_dir:                  Main github pull directory - ensures all required scripts are in 
#                           one place - "./eICU_featureization/" 
#                           
#save_dir:                  Location where script will create a new folder titled experiement_name
#                           and save all model related files to. - ideally, if you are a contributer to the github,
#                           do not make pull directory your save_dir (unless you want to make exception to git push)
#
## File Related Inputs: 
#
#merge_identifier:          The identifier/key column feature space and label space can merge on. Please refer to 
#                           the example file with "patientunitstayid" as the merge identifier. 
#
#feature_dir:               Directory and name of the CSV feature space. Must include merge identifier as a column (ideally
#                           the first column)
#
#label_dir:                 Directory and name of the CSV label space. Must include merge identifier as the first column. 
#
#outcomes:                  A c() character list of the class_1 (what the model should predict for) and class_0. 
#                           (ie) outcomes = c("dead", "alive") [c(class_1, class_0)] (ordering matters)
#                           ^ all metrics will be in relation to predicting for the "dead" (class_1) prediction task. [REQUIRED]
#
## Model Parameter Inputs (most parameter inputs can be left blank)
#
#experiment_name:           name the experiment (creates a new directory in save_dir with the experiment_name) [REQUIRED]
#                           Please select a unique intuitive name or files and folders will be overwritten 
#
#num_outer_loop:            specifies how many resamplings of the training and testing set to perform [DEFAULT = 1]
#
#cross_validation_k:        specifies how many k-fold cross validation to perform [DEFAULT = 10]. 
#
#cross_validation_repeats:  specifies how many fold resamplings to repeat [DEFAULT = 1]
#
#random_seed:               pick a random seed for repoducibility [DEFAULT = 1206]
#
#
## Extra_parameters
#
#Already_trained:           [TRUE/FALSE] Sometimes, you just want to plot again or re-run existing files. Set as true
#                           and ensure the save_dir and experiment_name point to the correct folder with .RData trained
#                           models to replot and resave CSVs. 

#Example Code (copy and paste this into a new script and source in the files in github directory)

# build_prototype(code_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/", 
#                 save_dir = "/storage/eICU/eICU_feature_extract/",
#                 
#                 merge_identifier = "patientunitstayid",
#                 feature_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_feature_space.csv",
#                 label_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_label.csv",
#                 outcomes = c("Expired", "Alive"),
#                 experiment_name = "test"
#                 )


#function below - do not alter - if so, submit an issue and make a fork/branch to work on the corrections

build_prototype <- function(code_dir, 
                            save_dir, 
                            
                            merge_identifier,
                            feature_dir,
                            label_dir,
                            outcomes,
                            
                            experiment_name, 
                            num_outer_loop,
                            cross_validation_k,
                            cross_validation_repeats,
                            random_seed,
                            Already_trained
                            ) {
  
  #directory checks
  if (!dir.exists(code_dir)) {
    stop("the CODE directory is not valid and/or does not exist. please double check")
  }
  
  if (!dir.exists(save_dir)) {
    stop("the SAVE directory is not valid and/or does not exist. please double check")
  }
  
  #check of required source scripts in directory
  if (dir.exists(code_dir)) {
    setwd(code_dir)
    if (!file.exists("required_custom_functions.R")) {
      stop("required_custom_functions.R is not within the code dir. Make sure code_dir is correct.\n
            If correct, do not rename or move code out of the code_dir and ensure required_custom_functions.R\n
            is in the directory")
    } else {
      source("required_custom_functions.R")
    }
    
    if (!file.exists("Performance_metric_plotting.R")) {
      stop("Performance_metric_plotting.R is not within the code dir. Make sure code_dir is correct.\n
            If correct, do not rename or move code out of the code_dir and ensure required_custom_functions.R\n
            is in the directory") 
    } else {
      source("Performance_metric_plotting.R")
    }
  }
  
  #experiment name checks
  if (missing(experiment_name)) {
    stop("Experiment name required - pick a unique name. files will be overwitten if otherwise.")
  }
  
  if (!is.character(experiment_name)) {
    experiment_name <- as.character(experiment_name)
    message("Experiment name was not entered as a character so was converted into a character.\nIf this was unintentional, double check")
  }
  
  #file related checks
  if (missing(merge_identifier)) {
    stop("Merge Identifier required - a key or common column between feature space and label space for merge purposes.")
  }
  
  if (!is.character(merge_identifier)) {
    message("Make sure the merge identifier is assinged as a character")
  }
  
  if (!file.exists(feature_dir)) {
    stop("feature space file is not within the feature_dir. Make sure feature_dir. is correct.") 
  }
  
  if (!file.exists("Performance_metric_plotting.R")) {
    stop("label space file is not within the label_dir. Make sure label_dir. is correct.") 
  }
  
  
  #parameter checks
  if (missing(num_outer_loop)) {
    num_outer_loop <- 1
  }
  
  if (missing(cross_validation_k)) {
    cross_validation_k <- 2
  }
  
  if (missing(cross_validation_repeats)) {
    cross_validation_repeats <- 1
  }
  
  if (missing(random_seed)) {
    random_seed <- 1206
  }
  
  if (missing(Already_trained)) {
    Already_trained <- FALSE
  }
  

package_list <- c("tidyverse", "data.table", "ranger", "randomForest", "caret", "MLmetrics", "Metrics", 
                  "EvaluationMeasures", "Rmisc", "xgboost", "glmnet")
load_packages(package_list)



#save_dir
dir <- paste0(save_dir, "/", experiment_name)
dir.create(dir)
setwd(dir)

table(read.csv(label_dir)$label) #used to identify what to label each outcome. Model always predicts for class_1. 

#outcome_names
bad_outcome <- outcomes[[1]] #class_1 - normally want to test for the unfavorable/unwanted outcome. 
#Class 1 is what you are predicting for. 
good_outcome <- outcomes[[2]]  #class_0

print(bad_outcome)
print(good_outcome)



############################# DO NOT CHANGE BELOW UNLESS IT BREAKS ####################################################


#Initializing some data frames to append to per loop. 
validation_results <- rbind()
validation_ROC <- rbind()
validation_PR <- rbind()

uncalibrated_test_ROC <- rbind()
calibrated_test_ROC <- rbind()

uncalibrated_test_probs_c1 <- rbind() 
uncalibrated_test_results_c1 <- rbind()
calibration_test_plot_c1 <- rbind() 
calibrated_test_results_c1 <- rbind() 
calibrated_test_probs_c1 <- rbind() 

test_class_balance <- rbind()



set.seed(random_seed)
for(j in 1:num_outer_loop) { 
  print(paste0("Outer--",j))
  
  #load in feature as full_feature_space and label
  full_feature_space <- read.csv(feature_dir)
  colnames(full_feature_space)[1] <- merge_identifier
  
  label <- read.csv(label_dir)
  colnames(label)[1] <- merge_identifier
  
  #converts existing labels according to what was specified in good_outcome and bad_outcome as class_1 and class_0
  label$label <- as.character(label$label)
  label$label[which(label$label == good_outcome)] <- "class_0"
  label$label[which(label$label == bad_outcome)] <-  "class_1"
  label$label <- as.factor(label$label)
  
  label <- label %>% dplyr::select(merge_identifier, label)
  
  #combine with label. ** label is assumed to be a factor with column header "label"
  full_feature_space <- merge(label, full_feature_space, by = merge_identifier)
  
  #only use complete cases. If it hasn't been imputed at this point, it is ok to remove. 
  full_feature_space <- full_feature_space[complete.cases(full_feature_space), ]
  
  df <- full_feature_space
  df$label <- as.factor(df$label)
  
  #remove_patientunitstayid
  rownames(df) <- df[, 1]
  df <- df[,-1]
  
  #Inner Loop
  intrain <- createDataPartition(y = df$label, p = 0.80, list = FALSE)
  
  #writes training space and testing space as csvs for future reproducible work if needed. 
  
  if (!Already_trained) {
  
  training <- df[intrain,]
  write.csv(training, paste0(dir, "/training_", j, ".csv"), row.names = F)
  testing <- df[-intrain,]
  write.csv(testing, paste0(dir, "/testing_", j, ".csv"), row.names = F)
  
  } else {
    training <- fread(paste0(dir, "/training_", j, ".csv"))
    
    training$label <- as.factor(training$label)
    
    testing <- fread(paste0(dir, "/testing_", j, ".csv"))
    
    testing$label <- as.factor(testing$label)
  }
  
  identifier <- experiment_name
  
  temp_balance <- c(unlist(as.list(table(testing$label))), identifier, j, "eICU")
  names(temp_balance)[c(3,4,5)] <- c("file", "outer_fold", "type")
  test_class_balance <- rbind(test_class_balance, temp_balance)
  
  label_balance <-prop.table(table(training$label))
  lower_class <- names(label_balance)[which.min(label_balance)]
  
  model_weights <- ifelse(training$label == lower_class,
                          label_balance[[1]],
                          label_balance[[2]])
  

  control = trainControl(method="repeatedcv",
                         number=cross_validation_k,
                         repeats=cross_validation_repeats,
                         summaryFunction = twoClassSummary,
                         classProbs = T,
                         savePredictions = T,
                         allowParallel = T)

  if (!Already_trained) {
  
    modelxgboost <- caret::train(label~., data = training, method="xgbTree", metric = 'ROC',trControl=control,
                                 weights = model_weights)
    print(paste0(identifier, "-", j, "-xg done"))
    save(modelxgboost,  file = paste(identifier,'_modelXG_iter_',j, '.Rdata', sep = ''))
    
    
    modelrf <- caret::train(label~., data = training, method="ranger", metric = 'ROC',trControl=control,
                            weights = model_weights)
    print(paste0(identifier, "-", j, "-rf done"))
    save(modelrf,  file = paste(identifier, '_modelRF_iter_',j, '.Rdata', sep = ''))
    
    modelglm <- caret::train(label~., data = training, method="glmnet", metric = 'ROC',trControl=control,
                             weights = model_weights)
    print(paste0(identifier, "-", j, "-glmnet done"))
    save(modelglm,  file = paste(identifier, '_modelGLM_iter_',j, '.Rdata', sep = ''))
  } else {
    print(paste0("loading trained models for outer loop: ", j) )
    load(file = paste(identifier, '_modelXG_iter_',j, '.Rdata', sep = ''))
    load(file = paste(identifier, '_modelRF_iter_',j, '.Rdata', sep = ''))
    load(file = paste(identifier, '_modelGLM_iter_',j, '.Rdata', sep = ''))
  
  }
  
  
  #####compile validation results. #####
  val_results <- resamples(list(xg = modelxgboost,
                                rf = modelrf,
                                glm = modelglm
  ))
  
  val_results <- val_results$values
  val_results$file <- identifier
  val_results$outer_fold <- j
  val_results <- val_results %>% dplyr::select(file, outer_fold, everything())
  validation_results <- rbind(validation_results, val_results)
  
  ##### Validation ROC probabilities #####
  xg_val_prob <- modelxgboost$pred
  xg_val_prob <- xg_val_prob[which(xg_val_prob$nrounds == modelxgboost$bestTune$nrounds & 
                                     xg_val_prob$eta == modelxgboost$bestTune$eta &
                                     xg_val_prob$max_depth == modelxgboost$bestTune$max_depth &
                                     xg_val_prob$gamma == modelxgboost$bestTune$gamma &
                                     xg_val_prob$colsample_bytree == modelxgboost$bestTune$colsample_bytree &
                                     xg_val_prob$min_child_weight == modelxgboost$bestTune$min_child_weight &
                                     xg_val_prob$subsample == modelxgboost$bestTune$subsample), ]
  rf_val_prob <- modelrf$pred
  rf_val_prob <- rf_val_prob[which(rf_val_prob$mtry == modelrf$bestTune$mtry &
                                     rf_val_prob$splitrule == modelrf$bestTune$splitrule &
                                     rf_val_prob$min.node.size == modelrf$bestTune$min.node.size), ]
  glm_val_prob <- modelglm$pred
  glm_val_prob <- glm_val_prob[which(glm_val_prob$alpha == modelglm$bestTune$alpha &
                                       glm_val_prob$lambda == modelglm$bestTune$lambda), ]
  folds <- unique(xg_val_prob$Resample)
  allprobs <- list(xg_val_prob, rf_val_prob, glm_val_prob)
  allprob_names <- c("xg", "rf", "glm")
  
  temp_count = 1
  for (val_prob in allprobs) {
    temp_name <- allprob_names[temp_count]
    
    model_ROC <- rbind()
    fold_count <- 1
    for (fold in folds) {
      temp <- val_prob[which(val_prob$Resample == fold), ] %>% dplyr::select(obs, class_1)
      temppred <- temp$class_1
      tempobs <- 1 * (as.character(temp$obs) == "class_1")
      eval_range <- seq(0, 1, 0.01) 
      
      tempTPR <- c()
      tempFPR <- c()
      for (eval in eval_range) {
        evaluate <- 1 * (temppred > eval)
        tempTPR <- c(tempTPR, EvaluationMeasures.TPR(Real = tempobs, Predicted = evaluate))
        tempFPR <- c(tempFPR, EvaluationMeasures.FPR(Real = tempobs, Predicted = evaluate))
      }
      tempROC <- cbind("file" = rep(identifier, length(eval_range)), 
                       "outer_fold" = rep(j, length(eval_range)),
                       "inner_fold" = rep(fold_count, length(eval_range)), 
                       "model" = temp_name, 
                       "TPR" = tempTPR, 
                       "FPR" = tempFPR)
      
      validation_ROC <- rbind(validation_ROC, tempROC)
      fold_count <- fold_count + 1
    }
    
    temp_count = temp_count + 1
  }
  
  
  ###### Validation Precision recall ######
  temp_count = 1
  for (val_prob in allprobs) {
    temp_name <- allprob_names[temp_count]
    
    model_PR <- rbind()
    fold_count <- 1
    for (fold in folds) {
      temp <- val_prob[which(val_prob$Resample == fold), ] %>% dplyr::select(obs, class_1)
      temppred <- temp$class_1
      tempobs <- 1 * (as.character(temp$obs) == "class_1")
      eval_range <- seq(0, 1, 0.01) 
      
      tempP <- c()
      tempR <- c()
      for (eval in eval_range) {
        evaluate <- 1 * (temppred > eval)
        tempP <- c(tempP, EvaluationMeasures.Precision(Real = tempobs, Predicted = evaluate))
        tempR <- c(tempR, EvaluationMeasures.Recall(Real = tempobs, Predicted = evaluate))
      }
      tempPR <- cbind("file" = rep(identifier, length(eval_range)), 
                      "outer_fold" = rep(j, length(eval_range)),
                      "inner_fold" = rep(fold_count, length(eval_range)), 
                      "model" = temp_name, 
                      "Precision" = tempP, 
                      "Recall" = tempR)
      
      validation_PR <- rbind(validation_PR, tempPR)
      fold_count <- fold_count + 1
    }
    
    temp_count = temp_count + 1
  }
  
  
  
  ##### Uncalibrated test set results class 1 #####
  allmodels <- list(modelxgboost, modelrf, modelglm)
  allmodel_names <- c("xg", "rf", "glm")
  temp_count <- 1
  
  for (model in allmodels) {
    model_name <- allmodel_names[temp_count]
    model_prediction <- predict(model, newdata = testing, type = "raw")
    model_prediction <- relevel(model_prediction, "class_1")
    testing$label <- relevel(testing$label, "class_1")
    test_confusion <- confusionMatrix(model_prediction, as.factor(testing$label))
    test_accuracy <- as.data.frame(test_confusion$overall)[1,]
    
    model_prediction_prob <- predict(model, newdata = testing, type = "prob")
    model_prediction_prob <- as.data.frame(cbind(testing$label,  model_prediction_prob, model_prediction))
    colnames(model_prediction_prob) <- c("obs", "class_0","class_1", "pred")
    
    twoclass <- twoClassSummary(model_prediction_prob, lev = levels(model_prediction_prob$obs))
    pr <- prSummary(model_prediction_prob, lev = levels(model_prediction_prob$obs))
    
    temp_uncalib_test_results <- c(identifier, j, model_name, test_accuracy, twoclass, pr)
    names(temp_uncalib_test_results) <- c("file", "outer_fold", "model", "Acc", "AUC", "Sens", "Spec", "PR_AUC", "Precision", "Recall", "F1" )
    uncalibrated_test_results_c1 <- rbind(uncalibrated_test_results_c1, temp_uncalib_test_results)
    
    #uncalibrated_test_probs
    model_prediction_prob$file <- identifier
    model_prediction_prob$outer_fold <- j
    model_prediction_prob$model <- model_name
    uncalibrated_test_probs_c1 <- rbind(uncalibrated_test_probs_c1, model_prediction_prob)
    
    temp_count <- temp_count + 1
  }
  
  
  #####Calibrated test set results class 1#####
  allmodels <- list(modelxgboost, modelrf, modelglm)
  allmodel_names <- c("xg", "rf", "glm")
  temp_count <- 1
  
  for (model in allmodels) {
    model_name <- allmodel_names[temp_count]
    model_prediction <- predict(model, newdata = testing, type = "raw")
    model_prediction <- relevel(model_prediction, "class_1")
    testing$label <- relevel(testing$label, "class_1")
    model_prediction_prob <- predict(model, newdata = testing, type = "prob")
    model_prediction_prob <- as.data.frame(cbind(testing$label,  model_prediction_prob, model_prediction))
    colnames(model_prediction_prob) <- c("obs", "class_0","class_1", "pred")
    
    cal_plot_data <- calibration(obs ~ class_1, data = model_prediction_prob, cuts = 10)$data
    cal_plot_data$type <- "uncalibrated"
    
    # plot(cal_plot_data$midpoint, cal_plot_data$Percent)
    
    
    temp_probs <- allprobs[[temp_count]] %>% dplyr::select(obs, class_0, class_1)
    temp_weights <- allprobs[[temp_count]] %>% dplyr::select(weights)
    temp_probs$obs <- ordered(temp_probs$obs, levels = c("class_0", "class_1"))
    lr_model = glm(obs ~ ., data = temp_probs,
                   family = quasibinomial)
    
    pred_probs <- temp_probs$class_1
    pred_obs <- 1*(temp_probs$obs == "class_1")
    
    idx <- duplicated(pred_probs)
    pred_probs <- pred_probs[!idx]
    pred_obs <- pred_obs[!idx]
    
    iso.model <- isoreg(pred_probs, pred_obs)
    
    lr_probs <- fit.isoreg(iso.model, model_prediction_prob$class_1)
    
    class_lr_probs = as.data.frame(cbind(as.character(testing$label), lr_probs))
    class_lr_probs$lr_probs <- as.numeric(as.character(class_lr_probs$lr_probs))
    colnames(class_lr_probs) <- c("obs", "class_1")
    class_lr_probs$obs <- ordered(class_lr_probs$obs, levels = c("class_1", "class_0"))
    
    cal_lr_plot_data = calibration(obs ~ class_1,
                                   data = class_lr_probs, class = "class_1")$data
    cal_lr_plot_data$type <- "calibrated"
    
    # plot(cal_lr_plot_data$midpoint, cal_lr_plot_data$Percent)
    
    temp_cal_plot_data <- rbind(cal_plot_data, cal_lr_plot_data)
    temp_cal_plot_data$model <- model_name
    temp_cal_plot_data$outer_fold <- j
    temp_cal_plot_data$file <- identifier
    
    calibration_test_plot_c1 <- rbind(calibration_test_plot_c1, temp_cal_plot_data)
    
    
    #find optimal cut off points that maximizes youden's index (youden's j statistics) by evaluating different thresholds
    #on the training set. 
    
    calib_train_probs <- fit.isoreg(iso.model, temp_probs$class_1)
    
    eval_range <- unique(calib_train_probs)
    eval_range <- eval_range[order(eval_range)]
    
    temp_sens <- c()
    temp_spec <- c()
    youdens <- c()
    for (eval in eval_range) {
      evaluate <- 1 * (calib_train_probs >= eval)
      train_obs <- 1 * (as.character(temp_probs$obs) == "class_1")
      
      tsens <- EvaluationMeasures.Sensitivity(Real = train_obs, Predicted = evaluate)
      tspec <- EvaluationMeasures.Specificity(Real = train_obs, Predicted = evaluate)
      
      temp_sens <- c(temp_sens, tsens)
      temp_spec <- c(temp_spec, tspec)
      
      youden <- tsens + tspec - 1
      youdens <- c(youdens, youden)
    }
    
    selected_sens <- temp_sens[which.max(youdens)]
    selected_spec <- temp_spec[which.max(youdens)]
    selected_eval <- eval_range[which.max(youdens)]
    # plot((1-temp_spec), temp_sens, 'l')
    # points((1-selected_spec), selected_sens, type = "p")
    
    #COMPUTE NEW STATISTICS FOR CALIBRATED MODEL
    
    model_prediction <- 1 * (class_lr_probs$class_1 >= selected_eval)
    model_prediction[model_prediction == 1] <- "class_1"
    model_prediction[model_prediction == "0"] <- "class_0"
    model_prediction <- as.factor(model_prediction)
    test_confusion <- confusionMatrix(model_prediction, as.factor(testing$label))
    test_accuracy <- as.data.frame(test_confusion$overall)[1,]
    
    class_lr_probs$class_0 <- 1 - class_lr_probs$class_1
    model_prediction_prob <- class_lr_probs %>% dplyr::select(class_0, class_1)
    model_prediction_prob <- as.data.frame(cbind(testing$label,  model_prediction_prob, as.character(model_prediction)))
    colnames(model_prediction_prob) <- c("obs", "class_0","class_1", "pred")
    model_prediction_prob$pred <- as.factor(model_prediction_prob$pred)
    model_prediction_prob$pred <- relevel(model_prediction_prob$pred, "class_1")
    
    twoclass <- twoClassSummary(model_prediction_prob, lev = levels(model_prediction_prob$obs))
    pr <- prSummary(model_prediction_prob, lev = levels(model_prediction_prob$obs))
    
    temp_calib_test_results <- c(identifier, j, model_name, test_accuracy, twoclass, pr)
    names(temp_calib_test_results) <- c("file", "outer_fold", "model", "Acc", "AUC", "Sens", "Spec", "PR_AUC", "Precision", "Recall", "F1" )
    calibrated_test_results_c1 <- rbind(calibrated_test_results_c1, temp_calib_test_results)
    
    #calibrated_test_probs
    model_prediction_prob$file <- identifier
    model_prediction_prob$outer_fold <- j
    model_prediction_prob$model <- model_name
    calibrated_test_probs_c1 <- rbind(calibrated_test_probs_c1, model_prediction_prob)
    
    
    temp_count <- temp_count + 1
  }
  
}

#saving all files for later use and use for plotting function. 
write.csv(validation_results, "validation_results.csv", row.names = F)
write.csv(validation_ROC, "validation_ROC.csv", row.names = F)
write.csv(validation_PR, "validation_PR.csv", row.names = F)

write.csv(uncalibrated_test_probs_c1, "uncalibrated_test_probs_c1.csv", row.names = F)
write.csv(uncalibrated_test_results_c1, "uncalibrated_test_results_c1.csv", row.names = F)

write.csv(calibration_test_plot_c1, "calibration_test_plot_c1.csv", row.names = F)
write.csv(calibrated_test_results_c1, "calibrated_test_results_c1.csv", row.names = F)
write.csv(calibrated_test_probs_c1, "calibrated_test_probs_c1.csv", row.names = F)

write.csv(test_class_balance, "test_class_balance.csv", row.names = F)

performance_metric_plotting(experiment_name = experiment_name, saved_file_location = dir)

}
