### Forward Selection feature selection code - 
###
### Usage: Selects features optimal for prediction task based on a greedy yet reliable forward selection method
###

#Author Han 
#updated
#V0.0 - 5/23/2021 - organized to be modular - released on Github

#### Output
#Creates directory with experiment name + FS_results with saved plots and CSV files:
# 1. AUROC per selected feature plot
# 2. AIC per selected feature plot
# 3. Relative likelihood per selected feature plot
# 4. CSV summary of AIC, AUROC, and relative likelihood per feature added.  
# 5. CSV of feature names selected
# 6. Final Selected features only CSV files to be used instead of the original feature for model training. (main result)

#### Input parameters & required files: ####

## Directory inputs (code will double check whether these directories exist and files within it are available)
#
#code_dir:                  Main github pull directory - ensures all required scripts are in 
#                           one place - "./eICU_featureization/" [REQUIRED]
#
#feature_dir:               Directory and name of the CSV feature space. Must include merge identifier as a column (ideally
#                           the first column) [REQUIRED]
#
#label_dir:                 Directory and name of the CSV label space. Must include merge identifier as the first column. 
#                           [REQUIRED]
#                           
#save_dir:                  Location where script will create a new folder titled experiement_name + FS_results
#                           and save all model related files to. [REQUIRED]
#
#outcomes:                  A c() character list of the class_1 (what the model should predict for) and class_0. 
#                           (ie) outcomes = c("dead", "alive") [c(class_1, class_0)] (ordering matters)
#                           ^ all metrics will be in relation to predicting for the "dead" (class_1) prediction task. [REQUIRED]
#
#
#merge_identifier:          The identifier/key column feature space and label space can merge on. Please refer to 
#                           the example file with "patientunitstayid" as the merge identifier. [OPTIONAL - default = "patientunitstayid"]
#
#
#experiment_name:           name the experiment (creates a new directory in save_dir with the experiment_name) [REQUIRED]
#                           Please select a unique intuitive name or files and folders will be overwritten 
#
#relative_like_threshold:   relative likelihood threshold to determine when to stop adding features. [OPTIONAL - default = 1]
#
#


################# example: ###################
# 
# code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")
# 
# forward_selection(code_dir = code_dir,
#                   feature_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_feature_space_4VIF.csv",
#                   label_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_label.csv",
#                   outcomes =  c("Expired", "Alive"),
#                   relative_like_threshold = 1,
#                   save_dir = "~/Downloads/",
#                   experiment_name = "same_name"
# )


forward_selection <- function(code_dir,
                              feature_dir, 
                              label_dir,
                              outcomes,
                              merge_identifier = "patientunitstayid",
                              experiment_name,
                              relative_like_threshold = 1, 
                              save_dir) {
  
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
  }
  
  #directory checks
  if (!dir.exists(save_dir)) {
    stop("the SAVE directory is not valid and/or does not exist. please double check")
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
    print("Merge Identifier missing- assumed to be patientunitstayid")
    merge_identifier = "patientunitstayid"
  }
  
  if (!is.character(merge_identifier)) {
    message("Make sure the merge identifier is assinged as a character")
  }
  
  if (!file.exists(feature_dir)) {
    stop("feature space file is not within the feature_dir. Make sure feature_dir. is correct.") 
  }
  
  if (!file.exists(label_dir)) {
    stop("label space file is not within the label_dir. Make sure label_dir is correct.") 
  }
  
  
  #parameter checks
  if (missing(relative_like_threshold)) {
    print("missing relative_like_threshold - set as default of 1. At higher thresholds, more feature will be selected. ")
    relative_like_threshold <- 1
  }
  

  package_list <- c("tidyverse", "data.table", "caret", "MLmetrics", "Metrics", 
                    "EvaluationMeasures", "Rmisc", "glmnet")
  load_packages(package_list)
  
  
  
  save_dir = paste0(save_dir, "/", experiment_name, "_FS_results/")
  print("files to be save into: ")
  print(save_dir)
  dir.create(save_dir)


  df = read.csv(feature_dir)
  df_original = as.data.frame(df)
  label = read.csv(label_dir)
  
  #remove unwanted features and standardize data
  df <- df[complete.cases(df), ]
  pids <- df %>% dplyr::select(all_of(merge_identifier))
  
  df <- BBmisc::normalize(df[,-1], method = "standardize", range = c(0, 1))
  df <- cbind(pids, df)
  
  #combine labels if separate
  df <- merge(label, df, by = merge_identifier)
  
  #turn label into binary 0, 1
  df$label <- as.character(df$label)
  df$label <- (df$label == outcomes[1]) * 1
  
  
  #split into training and testing for LM 
  set.seed(1206)
  intrain <- createDataPartition(y = df$label, p = 0.80, list = FALSE)
  training <- df[intrain,]
  train_pids <- training %>%  dplyr::select(patientunitstayid)
  training <- training[,-1]
  
  testing <- df[-intrain,]
  test_pids <- testing %>%  dplyr::select(patientunitstayid)
  testing <- testing[,-1]
  
  
  #first step to reduce feature number. If too large of a dimension, can't compute VIF. 
  #lets first do a correlation analysis to remove highly correlated variables first. 
  #
  training <- Filter(function(x) sd(x) != 0, training)
  
  correlationMatrix <- cor(training)
  # print(correlationMatrix)
  # heatmap(correlationMatrix, keep.dendro = F)
  highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.90, exact = T, names = T)
  # print(highlyCorrelated)
  
  training <- training %>% select(-all_of(highlyCorrelated))
  
  
  ######  FIRST PASS TO REMOVE ALIAS (DEPENDENT VARIABLES) ######
  # Build the linear model
  model1 <- lm(label ~., data = training)
  # Make predictions
  predictions <- model1 %>% predict(testing)
  # Model performance
  RMSE = RMSE(predictions, testing$label)
  R2 = R2(predictions, testing$label)
  
  ####VIF analysis dataframe to store iteration. 
  VIF_analysis <- as.data.frame(matrix(ncol = 5))
  colnames(VIF_analysis) <- c("RMSE", "R2", "num_features", "removed_feature", "VIF_removed")
  VIF_analysis$RMSE[1] <- RMSE
  VIF_analysis$R2[1] <- R2
  VIF_analysis$num_features[1] <- ncol(training) - 1
  VIF_analysis$removed_feature[1] <- NA
  
  # #compute alias of features in model.
  alias_ <- alias(model1)
  model1_alias <- alias_$Complete
  
  if (!is_empty(model1_alias)) {
    alias_names <- names(rowSums(model1_alias))
    
    #remove alias features from training data for the next iterations.
    training2 <- training[, -which(colnames(training) %in% alias_names)]
    # training2 <- training
    
    vif_features <- training2
  
  } else {
    vif_features <- training
  }
  
  vif_features$label <- as.factor(vif_features$label)
  vif_features$label <- NULL
  
  vif_features <- cbind(train_pids, vif_features)
  vif_features <- merge(label, vif_features, by = merge_identifier)
  vif_features$label = as.character(vif_features$label)
  vif_features$label[which(vif_features$label == outcomes[1])] = "class_1"
  vif_features$label[which(vif_features$label == outcomes[0])] = "class_0"
  
  VIF_names <- as.data.frame(colnames(vif_features)[-1])
  VIF_names[,1] <- as.character(VIF_names[,1])
  VIF_names$order <- 1
  VIF_names$category <- "none"
  colnames(VIF_names)[1] <- "x"
  
  ordering <- unique(VIF_names$order)
  
  vif_features <- as.data.frame(vif_features)
  threshold <- relative_like_threshold
  current_selected <- c()
  main_AIC <- rbind()
  for (i in 1:length(ordering)) {
    print(i) 
    
    temp_names <- VIF_names[which(VIF_names$order == i), ]
    temp_names <- temp_names[-which(temp_names$x == "label"), ]
    
    # temp_features <- vif_features %>% select("label", temp_names$x)
    temp_names_upper <- temp_names
    while (nrow(temp_names_upper) != 0) {
      
      temp_names <- temp_names_upper
      
      AIC1 <- rbind()
      while (nrow(temp_names) != 0) {
        # print(nrow(temp_names))
        # print(temp_names)
        
        if (temp_names$category[1] != "none") {
          variables <- c(temp_names$x[which(temp_names$category == temp_names$category[1])], current_selected)
          temp_features <- vif_features[, which(colnames(vif_features) %in% c("label", variables)), drop = F]
          temp_names <- temp_names[-which(temp_names$x %in% variables), , drop = F]
          
        } else {
          variables <- c(temp_names$x[1], current_selected)
          subset_variables <- c("label", variables)
          temp_features <- vif_features[, which(colnames(vif_features) %in% subset_variables), drop = F]
          temp_names <- temp_names[-which(temp_names$x %in% variables), , drop = F]
        }
        
        # print(variables)
        temp_features$label = as.factor(temp_features$label)
        control = trainControl(method="repeatedcv",
                               number=10,
                               repeats=1,
                               summaryFunction = twoClassSummary,
                               classProbs = T,
                               savePredictions = T)
        
        modelGlm <- caret::train(label~., data = temp_features, method="glm", metric = "ROC", trControl = control)
        
        tempAIC <- data.frame("variables" = paste(variables, collapse = ","), "aic" = modelGlm$finalModel$aic, "auc" = modelGlm$results$ROC)
        AIC1 <- rbind(AIC1, tempAIC)
        AIC1$variables <- as.character(AIC1$variables)
        
      }
      
      temp_selected <- c(as.character(AIC1$variables[which.min(AIC1$aic)]))
      raw_selected <- temp_selected
      
      if (temp_selected %like% ",") {
        # print("has multiple features")
        temp_selected <- strsplit(temp_selected, ",")[[1]]
      }
      
      selected_AIC <- AIC1[which(AIC1$variables %in% raw_selected), ]
      
      ## evaluate AIC ## 
      if (!is.null(main_AIC)) {
        prior_selected_AIC.value <- main_AIC$aic[nrow(main_AIC)]
        selected_AIC.value <- selected_AIC$aic
        
        rel_likelihood <- exp((prior_selected_AIC.value - selected_AIC.value)/2)
      } else {
        rel_likelihood <- threshold + 1
      }
      
      if (rel_likelihood > threshold) {
        selected_AIC$rel_like <- rel_likelihood
        current_selected <- unique(c(current_selected, temp_selected))
        temp_names_upper <- temp_names_upper[-which(temp_names_upper$x %in% temp_selected), ]
        print(selected_AIC)
        main_AIC <- rbind(main_AIC, selected_AIC)
      } else {
        temp_names_upper <- temp_names_upper[-which(temp_names_upper$x %in% temp_names_upper$x), ]
    
        print("below threshold - exited")
      }
      
      
    }
    
    print("observe current_selected and main_AIC")
    
  }
  
  current_selected_df = as.data.frame(current_selected, stringsAsFactors = F)
  colnames(current_selected_df) = "features"
  
  write.csv(current_selected_df, paste0(save_dir, "FS_selected_features_names.csv"), row.names = F)
  write.csv(main_AIC, paste0(save_dir, "FS_info.csv"), row.names = F)
  
  
  jpeg(file= paste0(save_dir, "/AIC_plot.jpeg"))
  plot(seq(1, nrow(main_AIC)), main_AIC$aic, 'l', xlab = "feature #", ylab = "AIC", main = paste0("Forward Selection AIC: threshold > ", relative_like_threshold))
  points(seq(1, nrow(main_AIC)), main_AIC$aic)
  dev.off()
  
  jpeg(file= paste0(save_dir, "/AUROC_plot.jpeg"))
  plot(seq(1, nrow(main_AIC)), main_AIC$auc, 'l', xlab = "feature #", ylab = "AUROC", main = "Forward Selection AUROC: threshold > 0.5")
  points(seq(1, nrow(main_AIC)), main_AIC$auc)
  dev.off()
  
  jpeg(file= paste0(save_dir, "/relative_likelihood_plot.jpeg"))
  plot(seq(1, nrow(main_AIC[-2,])), main_AIC$rel_like[-2], 'l', xlab = "feature #", ylab = "AUROC", main = "Forward Selection rel_like: threshold > 0.5")
       # ylim = c(0.1, 10))
  points(seq(1, nrow(main_AIC[-2,])), main_AIC$rel_like[-2])
  dev.off()
  
  var_list <- current_selected_df$features
  df_original =  df_original[, which(colnames(df_original) %in% c(merge_identifier, var_list))]
  df_original = as.data.frame(df_original)
  df_original =  df_original[, which(colnames(df_original) %in% c(merge_identifier, var_list))]
  write.csv(df_original, paste0(save_dir, "FS_feature_space.csv"), row.names = F)
}




