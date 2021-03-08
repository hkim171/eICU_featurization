### Feature Ranking Functions Auto-Pipeline
###
### Usage: feature ranking functions for predictor interpretability. 
### must be run after build_prototype is run. 

####Creates feature rankings for: 
### Completed
### 1. Random forest impurity ranking from build_prototype trained RData
###     This is only compatible with v0.3 build_prototype or higher. 
### 2. Random forest tree depth ranking
###
### 3. XGboost trained model default feature importance (uses gradient boosting feature importance ranking)
### 4. XGBoost Shapley feature ranking and summary. 
###
### In progress
### 5. Glmnet beta coefficients from build_prototype trained RData
### 6. Univariate GLM beta coefficients computed from top random forest ranked features
###
### All of the above rankings will be visualized and plots saved into the same plot directory as build_protoype. 

#Author Han 
#updated
#V0.0 - 3/2/2021 - random_forest_rank completed. 
#V0.1 - 3/3/2021 - GLM ranking in progress


########################################## RANDOM FOREST FEATURE RANKINGS #######################################################
#### Output
#Creates directory with saved ranking CSV and plots of each algorithm's feature ranking. 

#### RANDOM FOREST Input parameters & required files: ####

## Directory inputs (code will double check whether these directories exist and files within it are available)
#
#experiment_folder_dir      location of your experiment_name folder (parent folder of the folder built_prototype created)
#
#code_dir:                  Main github pull directory - ensures all required scripts are in 
#                           one place - "./eICU_featureization/" 
#                           

## Model Parameter Inputs (most parameter inputs can be left blank)
#
#experiment_name:           name the experiment (must be the same as the build_prototype as it will read from the
#                           experiment folder) [REQUIRED]
#
#num_outer_loop:            Same outerloop number as in build_protoype
#
#how_many_top_features:     Number of top ranked features to visualize. This is important to easily view the plot texts. 
#
## Extra_parameters
#
#use_trained_rf:            [TRUE/FALSE] 
#                           TRUE plots impurity based Random forest ranking from the trained Rdata
#                           saved from build_prototype
#                           FALSE trains and plots new tree depth based ranking using the package
#                           randomforestSRC. 



#example function calls
# experiment_folder_dir <- "/storage/eICU/eICU_feature_extract/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")
# source(paste0(code_dir, "/feature_ranking.R"))
# 
# random_forest_rank(experiment_folder_dir <- experiment_folder_dir, 
#                    code_dir <- code_dir,
#                    experiment_name = "test",
#                    num_outer_loop = 5,
#                    how_many_top_features = 50,
#                    use_trained_rf = TRUE)
# 
# random_forest_rank(experiment_folder_dir <- experiment_folder_dir, 
#                    code_dir <- code_dir,
#                    experiment_name = "test",
#                    num_outer_loop = 5,
#                    how_many_top_features = 50,
#                    use_trained_rf = FALSE)

########################################## RANDOM FOREST FEATURE RANKINGS #######################################################

random_forest_rank <- function(experiment_folder_dir, code_dir, experiment_name, num_outer_loop, how_many_top_features, use_trained_rf) {
  
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
  
  #set proper directory then make sure all random forest trained files exist. 
  
  if (missing(use_trained_rf)) {
    use_trained_rf <- FALSE
  }
  
  if (use_trained_rf) {
    rf_object_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/")
    if (dir.exists(rf_object_dir)) {
      setwd(rf_object_dir)
      for (i in 1:num_outer_loop) {
        rf_obj_name <- paste0(experiment_name, "_modelRF_iter_", i, ".Rdata")
        
        if (!file.exists(rf_obj_name)) {
          stop(paste0("make sure that experiment_folder_dir is the folder that the experiment_name folder sits.\nAlso 
               double check that ", rf_obj_name, " exists in ", rf_object_dir))
        } else {
          print(paste0(rf_obj_name, " found"))
        }
      }
    }
  }
  
  if (missing(num_outer_loop)) {
    stop("num_outer_loop must be specified - this is the parameter num_outer_loop for the build_prototype function")
  }
  
  
  package_list <- c("randomForestSRC", "tidyverse", "data.table", "ranger", "randomForest", "caret", "MLmetrics", "Metrics", 
                    "EvaluationMeasures", "Rmisc", "xgboost", "glmnet")
  load_packages(package_list)
  
  save_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/plots/")
 
  #if use_trained_rf is true, it will utilize the trained RF with ranking based on impurities of each node for each decision tree.
  #This compiles the results of all the outer folds and averages them with CI to present the comprehensive multi-fold feature
  #ranking result. 
  if (use_trained_rf) {
    
    for (i in 1:num_outer_loop) {
      setwd(rf_object_dir)
      rf_obj_name <- paste0(experiment_name, "_modelRF_iter_", i, ".Rdata")
      load(rf_obj_name)
      
      tempimp <- varImp(modelrf)$importance #this is based on Random forest impurity measurements. 
      tempimp <- rownames_to_column(tempimp)
      
      if (i == 1) {
        temp_rf_ranks <- tempimp
        temp_rf_ranks_long <- tempimp
      } else {
        temp_rf_ranks <- merge(temp_rf_ranks, tempimp, by = "rowname", all = T)
        temp_rf_ranks_long <- rbind(temp_rf_ranks_long, tempimp)
      }
 
    }
    
    if(how_many_top_features > nrow(temp_rf_ranks)) {
      how_many_top_features <- nrow(temp_rf_ranks)
    }
    
    temp_rf_ranks$avg <- rowSums(temp_rf_ranks[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    temp_rf_ranks <- temp_rf_ranks[order(-temp_rf_ranks$avg), ]
    
    temp_rf_ranks_subset <- temp_rf_ranks_long[which(temp_rf_ranks_long$rowname %in% temp_rf_ranks$rowname[c(1,1:how_many_top_features)]), ]
    
    # png(paste0(save_dir, "/top_", how_many_top_features,"_features_trained.png"), width = 7, height = 10,units = "in", res = 600)
    p <- ggplot(temp_rf_ranks_subset, aes(x = fct_reorder(rowname, Overall , .desc=FALSE), y = Overall, fill = factor(rowname))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," trained RF impurity ranked features\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("RF ranking based on node Impurity")
    # dev.off()
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_trained_impurity_ranked.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(temp_rf_ranks) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(temp_rf_ranks, paste0(save_dir,"/trained_RF_impurity_rankings.csv"), row.names = F)
    
  } else {
    
    if (!use_trained_rf) {
      rf_object_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/")
      if (dir.exists(rf_object_dir)) {
        setwd(rf_object_dir)
        for (i in 1:num_outer_loop) {
          rf_obj_name <- paste0("training_", i, ".csv")
          
          if (!file.exists(rf_obj_name)) {
            stop(paste0("make sure that experiment_folder_dir is the folder that the experiment_name folder sits.\nAlso 
               double check that ", rf_obj_name, " exists in ", rf_object_dir))
          } else {
            print(paste0(rf_obj_name, " found"))
          }
        }
      }
    }
    
    setwd(rf_object_dir)

    for (i in 1:num_outer_loop) {
      
      print(paste0("outer loop ", i, ": training and ranking using Random Forest SRC"))
      
      dat <- read.csv(paste0("training_", i, ".csv"))
      
      rf= rfsrc(label~.,   data = dat, 
                ntree = 5000, splitrule = 'gini',
                na.action = "na.omit")
      print(rf)
      
      max.Subtree = max.subtree(rf, conservative = F)  
      
      if (i == 1) {
        allvardepth = sort(max.Subtree$order[, 1]);
        allvardepth.df = data.frame(Variable=names(allvardepth),MinDepthMaxSubtree=allvardepth,row.names = NULL)
        
        allvardepth.all_iterR <- allvardepth.df
        allvardepth.all_iterM <- allvardepth.df

      } else {
        allvardepth = sort(max.Subtree$order[, 1]);
        allvardepth.df = data.frame(Variable=names(allvardepth),MinDepthMaxSubtree=allvardepth,row.names = NULL)
        
        allvardepth.all_iterM = merge(allvardepth.all_iterM, allvardepth.df, by = "Variable")
        allvardepth.all_iterR = rbind(allvardepth.all_iterR, allvardepth.df)
      }
      
    }
    
    allvardepth.all_iterM$avg <- rowSums(allvardepth.all_iterM[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    allvardepth.all_iterM <- allvardepth.all_iterM[order(allvardepth.all_iterM$avg), ]
    
    allvardepth.all_iterR_subset <- allvardepth.all_iterR[which(allvardepth.all_iterR$Variable %in% allvardepth.all_iterM$Variable[1:60]), ]
    
    p <- ggplot(allvardepth.all_iterR_subset, aes(x = fct_reorder(Variable, MinDepthMaxSubtree , .desc=TRUE), y = MinDepthMaxSubtree, fill = factor(Variable))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," RF SRC tree depth ranked features\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("RF ranking based on Feature Tree Depth")
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_trained_depth_ranked.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(allvardepth.all_iterM) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(allvardepth.all_iterM, paste0(save_dir,"/trained_RF_depth_rankings.csv"), row.names = F)
    
  }

}

########################################## RANDOM FOREST DONE #######################################################

########################################## XGBOOST FEATURE RANKINGS #######################################################
#### Output
#Creates directory with saved ranking CSV and plots of each algorithm's feature ranking. 

#### XGBOOST Input parameters & required files: ####

## Directory inputs (code will double check whether these directories exist and files within it are available)
#
#experiment_folder_dir      location of your experiment_name folder (parent folder of the folder built_prototype created)
#
#code_dir:                  Main github pull directory - ensures all required scripts are in 
#                           one place - "./eICU_featureization/" 
#                           

## Model Parameter Inputs (most parameter inputs can be left blank)
#
#experiment_name:           name the experiment (must be the same as the build_prototype as it will read from the
#                           experiment folder) [REQUIRED]
#
#num_outer_loop:            Same outerloop number as in build_protoype
#
#how_many_top_features:     Number of top ranked features to visualize. This is important to easily view the plot texts. 
#
## Extra_parameters
#
#shap:                      [TRUE/FALSE] - defaults to FALSE
#                           TRUE plots trained XGBOOST feature importance
#                           FALSE plots SHAPforxgboost feature importance and summary


#example function calls
# experiment_folder_dir <- "/storage/eICU/eICU_feature_extract/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")
# # source(paste0(code_dir, "/feature_ranking.R"))
# 
# XG_rank(experiment_folder_dir <- experiment_folder_dir,
#                    code_dir <- code_dir,
#                    experiment_name = "test",
#                    num_outer_loop = 5,
#                    how_many_top_features = 50,
#                    shap = TRUE)
# 
# XG_rank(experiment_folder_dir <- experiment_folder_dir,
#                    code_dir <- code_dir,
#                    experiment_name = "test",
#                    num_outer_loop = 5,
#                    how_many_top_features = 50,
#                    shap = FALSE)

#XGboost rankings from gradient boosted feature importance from trained model using build_prototype

#will always use trained RF - shap = T, shap = F will be used to integrate shapley feature explanations. 

XG_rank <- function(experiment_folder_dir, code_dir, experiment_name, num_outer_loop, how_many_top_features, shap) {
  
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
  
  #set proper directory then make sure all XG trained files exist. 
  xg_object_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/")
  if (dir.exists(xg_object_dir)) {
    setwd(xg_object_dir)
    for (i in 1:num_outer_loop) {
      xg_obj_name <- paste0(experiment_name, "_modelXG_iter_", i, ".Rdata")
      
      if (!file.exists(xg_obj_name)) {
        stop(paste0("make sure that experiment_folder_dir is the folder that the experiment_name folder sits.\nAlso 
             double check that ", xg_obj_name, " exists in ", xg_object_dir))
      } else {
        print(paste0(xg_obj_name, " found"))
      }
    }
  }

  
  if (missing(num_outer_loop)) {
    stop("num_outer_loop must be specified - this is the parameter num_outer_loop for the build_prototype function")
  }
  
  
  package_list <- c("randomForestSRC", "tidyverse", "data.table", "ranger", "randomForest", "caret", "MLmetrics", "Metrics", 
                    "EvaluationMeasures", "Rmisc", "xgboost", "glmnet")
  load_packages(package_list)
  
  save_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/plots/")
  
  #first is using built in Variable importance from the xgboost training schema using Caret - 
  if (!shap) {
    
    for (i in 1:num_outer_loop) {
      setwd(xg_object_dir)
      xg_obj_name <- paste0(experiment_name, "_modelXG_iter_", i, ".Rdata")
      load(xg_obj_name)
      
      tempimp <- varImp(modelxgboost)$importance #this is based on Random forest impurity measurements. 
      tempimp <- rownames_to_column(tempimp)
      
      if (i == 1) {
        temp_xg_ranks <- tempimp
        temp_xg_ranks_long <- tempimp
      } else {
        temp_xg_ranks <- merge(temp_xg_ranks, tempimp, by = "rowname", all = T)
        temp_xg_ranks_long <- rbind(temp_xg_ranks_long, tempimp)
      }
      
    }
    
    if(how_many_top_features > nrow(temp_xg_ranks)) {
      how_many_top_features <- nrow(temp_xg_ranks)
    }
    
    temp_xg_ranks$avg <- rowSums(temp_xg_ranks[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    temp_xg_ranks <- temp_xg_ranks[order(-temp_xg_ranks$avg), ]
    
    temp_xg_ranks_subset <- temp_xg_ranks_long[which(temp_xg_ranks_long$rowname %in% temp_xg_ranks$rowname[c(1,1:how_many_top_features)]), ]
    
    # png(paste0(save_dir, "/top_", how_many_top_features,"_features_trained.png"), width = 7, height = 10,units = "in", res = 600)
    p <- ggplot(temp_xg_ranks_subset, aes(x = fct_reorder(rowname, Overall , .desc=FALSE), y = Overall, fill = factor(rowname))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," trained XGBoost features\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("XGBoost ranking")
    # dev.off()
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_XGBoost_rankings.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(temp_xg_ranks) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(temp_xg_ranks, paste0(save_dir,"/trained_XGBoost_rankings.csv"), row.names = F)
    
  } else {
    
    package_list <- c("SHAPforxgboost")
    load_packages(package_list)
    
    for (i in 1:num_outer_loop) {
      setwd(xg_object_dir)
      xg_obj_name <- paste0(experiment_name, "_modelXG_iter_", i, ".Rdata")
      load(xg_obj_name)
      
      final_model <- modelxgboost$finalModel
      train <- modelxgboost$trainingData

      shap_values <- shap.values(xgb_model = final_model, X_train = as.matrix(train[,-1]))
      tempimp <- as.data.frame(shap_values$mean_shap_score)
      tempimp <- rownames_to_column(tempimp)
      
      
      if (i == 1) {
        temp_xg_ranks <- tempimp
        temp_xg_ranks_long <- tempimp
      } else {
        temp_xg_ranks <- merge(temp_xg_ranks, tempimp, by = "rowname", all = T)
        temp_xg_ranks_long <- rbind(temp_xg_ranks_long, tempimp)
      }
      
    }
    
    if(how_many_top_features > nrow(temp_xg_ranks)) {
      how_many_top_features <- nrow(temp_xg_ranks)
    }
    
    temp_xg_ranks$avg <- rowSums(temp_xg_ranks[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    temp_xg_ranks <- temp_xg_ranks[order(-temp_xg_ranks$avg), ]
    
    temp_xg_ranks_subset <- temp_xg_ranks_long[which(temp_xg_ranks_long$rowname %in% temp_xg_ranks$rowname[c(1,1:how_many_top_features)]), ]
    colnames(temp_xg_ranks_subset)[2] <- "Overall"
    
    # png(paste0(save_dir, "/top_", how_many_top_features,"_features_trained.png"), width = 7, height = 10,units = "in", res = 600)
    p <- ggplot(temp_xg_ranks_subset, aes(x = fct_reorder(rowname, Overall , .desc=FALSE), y = Overall, fill = factor(rowname))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," trained XGBoost SHAP VALUE Ranking\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("XGBoost SHAP ranking")
    # dev.off()
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_XGBoost_SHAP_rankings.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(temp_xg_ranks) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(temp_xg_ranks, paste0(save_dir,"/trained_XGBoost_SHAP_rankings.csv"), row.names = F)
    
    top_features <- unique(temp_xg_ranks_subset$rowname)

    for (i in 1:num_outer_loop) {
      setwd(xg_object_dir)
      xg_obj_name <- paste0(experiment_name, "_modelXG_iter_", i, ".Rdata")
      load(xg_obj_name)
      
      final_model <- modelxgboost$finalModel
      train <- modelxgboost$trainingData
      
      shap_long <- shap.prep(xgb_model = final_model, X_train = as.matrix(train[,-1]))
      shap_long <- shap_long[which(shap_long$variable %in% top_features), ] %>% droplevels(.)

      
      if (i == 1) {
        temp_xg_ranks_long <- shap_long
      } else {
        temp_xg_ranks_long <- rbind(temp_xg_ranks_long, shap_long)
      }
      
    }
     
    temp_xg_ranks_long_test <- temp_xg_ranks_long
    
    new_averages <- aggregate(temp_xg_ranks_long$mean_value, list(temp_xg_ranks_long$variable), mean)
    colnames(new_averages)[1] <- "variable"
    temp_xg_ranks_long <- merge(temp_xg_ranks_long, new_averages, by = "variable")
    
    temp_xg_ranks_long <- temp_xg_ranks_long[,-6]
    colnames(temp_xg_ranks_long)[6] <- "mean_value"
    
    temp_xg_ranks_long$variable <- fct_reorder(temp_xg_ranks_long$variable, -temp_xg_ranks_long$mean_value)
    
    p <- shap.plot.summary(temp_xg_ranks_long)
    p <- p + ggtitle(paste0(experiment_name, ": top ",how_many_top_features," trained XGBoost SHAP Summary \nSina plots over ", num_outer_loop, " outer fold iterations"))

    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_SHAP_Summary_plot.png"), width = 7, height = 10, dpi = 600, units = "in")

    
  }
  
}



########################################## GLM FEATURE RANKINGS #######################################################
#glmnet coefficients from trained model

A <- coef(modelglm$finalModel, modelglm$finalModel$lambdaOpt)
names <- as.data.frame(as.character(A@Dimnames[[1]]))
colnames(names) <- "features"
names$features <- as.character(names$features)

names$coefficient <- NA
names$coefficient[A@i + 1] <- A@x

#glm coefficients from curated list of top random forest features


#glm coefficients from curated list of top xgboost features

GLM_rank <- function(experiment_folder_dir, code_dir, experiment_name, num_outer_loop, how_many_top_features, use_trained_rf) {
  
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
  
  #set proper directory then make sure all random forest trained files exist. 
  
  if (missing(use_trained_rf)) {
    use_trained_rf <- FALSE
  }
  
  if (use_trained_rf) {
    rf_object_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/")
    if (dir.exists(rf_object_dir)) {
      setwd(rf_object_dir)
      for (i in 1:num_outer_loop) {
        rf_obj_name <- paste0(experiment_name, "_modelRF_iter_", i, ".Rdata")
        
        if (!file.exists(rf_obj_name)) {
          stop(paste0("make sure that experiment_folder_dir is the folder that the experiment_name folder sits.\nAlso 
               double check that ", rf_obj_name, " exists in ", rf_object_dir))
        } else {
          print(paste0(rf_obj_name, " found"))
        }
      }
    }
  }
  
  if (missing(num_outer_loop)) {
    stop("num_outer_loop must be specified - this is the parameter num_outer_loop for the build_prototype function")
  }
  
  
  package_list <- c("randomForestSRC", "tidyverse", "data.table", "ranger", "randomForest", "caret", "MLmetrics", "Metrics", 
                    "EvaluationMeasures", "Rmisc", "xgboost", "glmnet")
  load_packages(package_list)
  
  save_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/plots/")
  
  #if use_trained_rf is true, it will utilize the trained RF with ranking based on impurities of each node for each decision tree.
  #This compiles the results of all the outer folds and averages them with CI to present the comprehensive multi-fold feature
  #ranking result. 
  if (use_trained_rf) {
    
    for (i in 1:num_outer_loop) {
      setwd(rf_object_dir)
      rf_obj_name <- paste0(experiment_name, "_modelRF_iter_", i, ".Rdata")
      load(rf_obj_name)
      
      tempimp <- varImp(modelrf)$importance #this is based on Random forest impurity measurements. 
      tempimp <- rownames_to_column(tempimp)
      
      if (i == 1) {
        temp_rf_ranks <- tempimp
        temp_rf_ranks_long <- tempimp
      } else {
        temp_rf_ranks <- merge(temp_rf_ranks, tempimp, by = "rowname", all = T)
        temp_rf_ranks_long <- rbind(temp_rf_ranks_long, tempimp)
      }
      
    }
    
    if(how_many_top_features > nrow(temp_rf_ranks)) {
      how_many_top_features <- nrow(temp_rf_ranks)
    }
    
    temp_rf_ranks$avg <- rowSums(temp_rf_ranks[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    temp_rf_ranks <- temp_rf_ranks[order(-temp_rf_ranks$avg), ]
    
    temp_rf_ranks_subset <- temp_rf_ranks_long[which(temp_rf_ranks_long$rowname %in% temp_rf_ranks$rowname[c(1,1:how_many_top_features)]), ]
    
    # png(paste0(save_dir, "/top_", how_many_top_features,"_features_trained.png"), width = 7, height = 10,units = "in", res = 600)
    p <- ggplot(temp_rf_ranks_subset, aes(x = fct_reorder(rowname, Overall , .desc=FALSE), y = Overall, fill = factor(rowname))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," trained RF impurity ranked features\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("RF ranking based on node Impurity")
    # dev.off()
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_trained_impurity_ranked.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(temp_rf_ranks) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(temp_rf_ranks, paste0(save_dir,"/trained_RF_impurity_rankings.csv"), row.names = F)
    
  } else {
    
    if (!use_trained_rf) {
      rf_object_dir <- paste0(experiment_folder_dir, "/", experiment_name, "/")
      if (dir.exists(rf_object_dir)) {
        setwd(rf_object_dir)
        for (i in 1:num_outer_loop) {
          rf_obj_name <- paste0("training_", i, ".csv")
          
          if (!file.exists(rf_obj_name)) {
            stop(paste0("make sure that experiment_folder_dir is the folder that the experiment_name folder sits.\nAlso 
               double check that ", rf_obj_name, " exists in ", rf_object_dir))
          } else {
            print(paste0(rf_obj_name, " found"))
          }
        }
      }
    }
    
    setwd(rf_object_dir)
    
    for (i in 1:num_outer_loop) {
      
      print(paste0("outer loop ", i, ": training and ranking using Random Forest SRC"))
      
      dat <- read.csv(paste0("training_", i, ".csv"))
      
      rf= rfsrc(label~.,   data = dat, 
                ntree = 5000, splitrule = 'gini',
                na.action = "na.omit")
      print(rf)
      
      max.Subtree = max.subtree(rf, conservative = F)  
      
      if (i == 1) {
        allvardepth = sort(max.Subtree$order[, 1]);
        allvardepth.df = data.frame(Variable=names(allvardepth),MinDepthMaxSubtree=allvardepth,row.names = NULL)
        
        allvardepth.all_iterR <- allvardepth.df
        allvardepth.all_iterM <- allvardepth.df
        
      } else {
        allvardepth = sort(max.Subtree$order[, 1]);
        allvardepth.df = data.frame(Variable=names(allvardepth),MinDepthMaxSubtree=allvardepth,row.names = NULL)
        
        allvardepth.all_iterM = merge(allvardepth.all_iterM, allvardepth.df, by = "Variable")
        allvardepth.all_iterR = rbind(allvardepth.all_iterR, allvardepth.df)
      }
      
    }
    
    allvardepth.all_iterM$avg <- rowSums(allvardepth.all_iterM[, c(2:(num_outer_loop + 1))]) / num_outer_loop
    
    allvardepth.all_iterM <- allvardepth.all_iterM[order(allvardepth.all_iterM$avg), ]
    
    allvardepth.all_iterR_subset <- allvardepth.all_iterR[which(allvardepth.all_iterR$Variable %in% allvardepth.all_iterM$Variable[1:60]), ]
    
    p <- ggplot(allvardepth.all_iterR_subset, aes(x = fct_reorder(Variable, MinDepthMaxSubtree , .desc=TRUE), y = MinDepthMaxSubtree, fill = factor(Variable))) + geom_boxplot() + coord_flip() + theme(legend.position = "none") +
      ggtitle(paste0(experiment_name, ": top ",how_many_top_features," RF SRC tree depth ranked features\nBoxplot over ", num_outer_loop,
                     " outer fold iterations")) +
      xlab("variable names") + ylab("RF ranking based on Feature Tree Depth")
    
    setwd(save_dir)
    ggsave(plot = p, filename = paste0("top_", how_many_top_features,"_features_trained_depth_ranked.png"), width = 7, height = 10, dpi = 600, units = "in")
    
    colnames(allvardepth.all_iterM) <- c("feature_name", paste0("outer_loop_", seq(from = 1, to = num_outer_loop, by = 1)), "mean")
    write.csv(allvardepth.all_iterM, paste0(save_dir,"/trained_RF_depth_rankings.csv"), row.names = F)
    
  }
  
}



