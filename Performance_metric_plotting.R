#All in one Plotting 

#replaces infinites in a dataframe with replace_with
replace_infinite <- function(df, replace_with) {
  df <- do.call(data.frame,lapply(df, function(x) replace(x, is.infinite(x),replace_with)))
  return(df)
}

#replaces NAN in a dataframe with replace_with
replace_nan <- function(df, replace_with) {
  df <- do.call(data.frame,lapply(df, function(x) replace(x, is.nan(x),replace_with)))
  return(df)
}

#replaces NAs in a dataframe with replace_with
replace_NA <- function(df, replace_with) { 
  df[is.na(df)] = replace_with
  return(df)
}

#Counts number of NAs in a dataframe
check_na <- function(df) {
  table(is.na(df))
}

#counts number of Inf and -Inf in a dataframe
check_infinite <- function(df){
  yes <- 0
  no <- 0
  for (i in 1:ncol(df)) {
    tempyes <- sum(is.infinite(df[,i]))
    yes <- yes + tempyes
    no <- no + (length(df[,1]) - yes)
  }
  df <- data.frame(no, yes, row.names = "")
  colnames("FALSE", "TRUE")
  print(as.data.frame(df))
}


performance_metric_plotting <- function(experiment_name, saved_file_location) {

library(tidyverse)
library(randomForestSRC)
library(ranger)
library(randomForest)
library(caret)
library(MLmetrics)
library(Metrics)
library(EvaluationMeasures)
library(Rmisc)
library(data.table)
library(caretEnsemble)

# source("~/useful_functions.R")

validation_results <- rbind()#done
validation_ROC <- rbind()#done
validation_PR <- rbind()#done

uncalibrated_test_ROC <- rbind()
uncalibrated_test_probs_c1 <- rbind() #done
uncalibrated_test_results_c1 <- rbind()#done

calibration_test_plot_c1 <- rbind() #done
calibrated_test_results_c1 <- rbind() #done
calibrated_test_probs_c1 <- rbind() #done

test_class_balance <- rbind()

main_id <- experiment_name
    
identifier <- experiment_name
subset <- experiment_name
dir <- saved_file_location
setwd(dir)  

identifiers <- c(experiment_name)
subsets <- c(experiment_name)

##### read in saved files #####

#validation_results table
temp_validation_results <- read.csv("validation_results.csv", stringsAsFactors = F)
temp_validation_results$identifier <- identifier
temp_validation_results$subset <- subset
validation_results <- rbind.fill(validation_results, temp_validation_results)
rm(temp_validation_results)

#validation ROC table
temp_validation_ROC <- read.csv("validation_ROC.csv", stringsAsFactors = F)
temp_validation_ROC$identifier <- identifier
temp_validation_ROC$subset <- subset
validation_ROC <- rbind.fill(validation_ROC, temp_validation_ROC)
rm(temp_validation_ROC)

#validation PR table
temp_validation_PR <- read.csv("validation_PR.csv", stringsAsFactors = F)
temp_validation_PR$identifier <- identifier
temp_validation_PR$subset <- subset
validation_PR <- rbind.fill(validation_PR, temp_validation_PR)
rm(temp_validation_PR)

#uncalibrated_test_probs_c1
temp_uncalibrated_test_probs_c1 <- read.csv("uncalibrated_test_probs_c1.csv", stringsAsFactors = F)
temp_uncalibrated_test_probs_c1$identifier <- identifier
temp_uncalibrated_test_probs_c1$subset <- subset
uncalibrated_test_probs_c1 <- rbind.fill(uncalibrated_test_probs_c1, temp_uncalibrated_test_probs_c1)
rm(temp_uncalibrated_test_probs_c1)

#uncalibrated_test_results_c1
temp_uncalibrated_test_results_c1 <- read.csv("uncalibrated_test_results_c1.csv", stringsAsFactors = F)
temp_uncalibrated_test_results_c1$identifier <- identifier
temp_uncalibrated_test_results_c1$subset <- subset
uncalibrated_test_results_c1 <- rbind.fill(uncalibrated_test_results_c1, temp_uncalibrated_test_results_c1)
rm(temp_uncalibrated_test_results_c1)

#calibration_test_plot_c1
temp_calibration_test_plot_c1 <- read.csv("calibration_test_plot_c1.csv", stringsAsFactors = F)
temp_calibration_test_plot_c1$identifier <- identifier
temp_calibration_test_plot_c1$subset <- subset
calibration_test_plot_c1 <- rbind.fill(calibration_test_plot_c1, temp_calibration_test_plot_c1)
rm(temp_calibration_test_plot_c1)

#calibrated_test_results_c1
temp_calibrated_test_results_c1 <- read.csv("calibrated_test_results_c1.csv", stringsAsFactors = F)
temp_calibrated_test_results_c1$identifier <- identifier
temp_calibrated_test_results_c1$subset <- subset
calibrated_test_results_c1 <- rbind.fill(calibrated_test_results_c1, temp_calibrated_test_results_c1)
rm(temp_calibrated_test_results_c1)

#calibrated_test_probs_c1
temp_calibrated_test_probs_c1 <- read.csv("calibrated_test_probs_c1.csv", stringsAsFactors = F)
temp_calibrated_test_probs_c1$identifier <- identifier
temp_calibrated_test_probs_c1$subset <- subset
calibrated_test_probs_c1 <- rbind.fill(calibrated_test_probs_c1, temp_calibrated_test_probs_c1)
rm(temp_calibrated_test_probs_c1)

#test_class_balance
temp_test_class_balance <- read.csv("test_class_balance.csv", stringsAsFactors = F)
temp_test_class_balance$identifier <- identifier
temp_test_class_balance$subset <- subset
test_class_balance <- rbind.fill(test_class_balance, temp_test_class_balance)
rm(temp_test_class_balance)

folds_outer <- length(unique(calibrated_test_results_c1$outer_fold))

#############################################################################################

dir <- paste0(saved_file_location, "/plots/")
dir.create(dir)
setwd(dir)  


xval_results_backup <- validation_results 

##### Plots ##### 
##### /---Calibration plot: test result for class 1 (bad outcome)##### 
calibration_test_plot <- as.data.frame(calibration_test_plot_c1)
calibrated_test_results <- as.data.frame(calibrated_test_results_c1)
uncalibrated_test_results <- as.data.frame(uncalibrated_test_results_c1)
uncalibrated_test_probs <- as.data.frame(uncalibrated_test_probs_c1)
uncalibrated_test_probs$type <- "uncalibrated"
calibrated_test_probs <- as.data.frame(calibrated_test_probs_c1)
calibrated_test_probs$type <- "calibrated"

test_probs <- rbind.fill(uncalibrated_test_probs, calibrated_test_probs)
balance <- as.data.frame(test_class_balance)
balance <- balance[which(balance$type == "eICU"), ] %>% group_by(file) %>% dplyr::slice(which.min(outer_fold))
balance$class_0 <- as.numeric(as.character(balance$class_0))
balance$class_1 <- as.numeric(as.character(balance$class_1))
balance$class_0_prop <- trunc(((balance$class_0) / (balance$class_0 + balance$class_1))*10^2)/10^2
balance$class_1_prop <- trunc(((balance$class_1) / (balance$class_0 + balance$class_1))*10^2)/10^2

savedir <- paste0(dir)
dir.create(savedir)
setwd(savedir)  

for (identifier in identifiers) {
  
  per_identifier_plot <- list()
  
  for (subset in subsets) {
    allmodel_names <- unique(calibration_test_plot[which(calibration_test_plot$subset == subset), ]$model)
    
    for (model_name in allmodel_names) {
      calib <- calibration_test_plot[which(calibration_test_plot$identifier == identifier & 
                                             calibration_test_plot$model == model_name &
                                             calibration_test_plot$subset == subset), ]
      
      
      
      temp_probs <- test_probs[which(test_probs$identifier == identifier &
                                       test_probs$model == model_name&
                                       test_probs$subset == subset), ]
      temp_balance <- balance[which(balance$identifier == identifier &
                                      balance$subset == subset), ]
      
      temp_probs <- temp_probs %>% dplyr::select(class_1, type, outer_fold)
      temp_probs$type_fold <- paste0(temp_probs$type, "_", temp_probs$outer_fold)
      
      
      br <- seq(0,100,by=10)
      
      hist_data <- rbind()
      
      for (typefold in unique(temp_probs$type_fold)) {
        tprobs <- temp_probs[which(temp_probs$type_fold == typefold), ]
        freq <- hist(tprobs$class_1 * 100, breaks = br, include.lowest = T, plot = F)
        counts <- (freq$counts/sum(freq$counts)) * 100
        mids <- freq$mids
        
        tempdata <- cbind(mids, counts, rep(tprobs$type[1], length(mids)))
        hist_data <- rbind(hist_data, tempdata)
        
      }
      hist_data <- as.data.frame(hist_data, stringAsFactors = F)
      hist_data$mids <- as.numeric(as.character(hist_data$mids))
      hist_data$counts <- as.numeric(as.character(hist_data$counts))
      colnames(hist_data)[3] <- "type"
      
      hist_data <- summarySE(hist_data, measurevar="counts", groupvars=c("type", "mids"))
      
      calib_plot <- ggplot(calib, aes(x = midpoint, y = Percent, group = as.factor(type), color = as.factor(type))) + geom_smooth(method = "loess") +
        geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = 2) +
        labs(subtitle = paste0(identifier, "--",subset,"--", model_name,"\nclass_1: " , temp_balance$class_1, " (",temp_balance$class_1_prop,")", "; class_0: ", temp_balance$class_0, " (",temp_balance$class_0_prop,")")) + 
        scale_color_discrete(name = "") + theme(
          plot.title = element_text(size = 20, face = "bold"),
          plot.subtitle = element_text(size = 12),
          plot.caption = element_text(size = 11)
        ) +
        ylim(0,100) +
        xlim(0,100)
      
      ggsave(paste0("eICU_calib_plot_", identifier, "_", subset, "_", model_name, ".png"), width = 6, height = 5, dpi = 300, units = "in")
      
      
      calib_results <- calibrated_test_results[which(calibrated_test_results$identifier == identifier & 
                                                       calibrated_test_results$model == model_name &
                                                       calibrated_test_results$subset == subset), ]
      calib_results$type <- "calibrated"
      
      uncalib_results <- uncalibrated_test_results[which(uncalibrated_test_results$identifier == identifier & 
                                                           uncalibrated_test_results$model == model_name &
                                                           uncalibrated_test_results$subset == subset), ]
      uncalib_results$type <- "uncalibrated"
      
      test_results <- rbind.fill(calib_results, uncalib_results)
      
      metrics <- c("Acc", "AUC", "Sens", "Spec")
      
      temp_plots <- list()
      for (metric in metrics) {
        
        temp_plots[[metric]] <- ggplot(test_results, aes(x = as.numeric(as.character(.data[[metric]])), y = type, group = factor(type), color = factor(type))) + geom_boxplot() + 
          ylab(NULL) + xlab(NULL) + ggtitle(metric) + theme(legend.position = "none") + xlim(0.5, 1)
        
      }
      
      x <- gridExtra::grid.arrange(temp_plots$Acc, temp_plots$AUC, temp_plots$Sens, temp_plots$Spec, ncol = 2)
      # png(paste0(identifier, "_", model_name, "_test_calibration.png"), width = 16, height = 6, units = "in", res = 300)
      per_identifier_plot[[model_name]] <- gridExtra::grid.arrange(calib_plot, x, ncol = 2)
      temp_plot_perid <- gridExtra::grid.arrange(calib_plot, x, ncol = 2)
      # dev.off()
      
      png(paste0(identifier, "_", subset, "_", model_name, "_test_calibration_c1.png"), width = 12, height = 5, units = "in", res = 300)
      gridExtra::grid.arrange(temp_plot_perid, ncol = 1)
      dev.off()
      
    }
    
  }
  
  png(paste0(identifier, "_test_calibration_c1.png"), width = 11, height = 18, units = "in", res = 300)
  gridExtra::grid.arrange(per_identifier_plot$xg, per_identifier_plot$rf, per_identifier_plot$glm, ncol = 1)
  dev.off()
  
}


##### ROCPLOTS #####

##### test ROC plots: calibrated ####

savedir <- paste0(dir)
dir.create(savedir)
setwd(savedir)  

test_probs <- calibrated_test_probs_c1
test_ROC <- c()
for (identifier in identifiers) {
  for (subset in subsets) {
    allmodel_names <- unique(test_probs[which(test_probs$subset == subset), ]$model)
    for (model_name in allmodel_names) {
      for (outer in unique(test_probs$outer_fold)) {
        
        temp_probs <- test_probs[which(test_probs$identifier == identifier &
                                         test_probs$model == model_name &
                                         test_probs$outer_fold == outer &
                                         test_probs$subset == subset), ]
        temppred <- temp_probs$class_1
        tempobs <- 1 * (as.character(temp_probs$obs) == "class_1")
        eval_range <- seq(0, 1, 0.01) 
        tempTPR <- c()
        tempFPR <- c()
        
        for (eval in eval_range) {
          evaluate <- 1 * (temppred > eval)
          tempTPR <- c(tempTPR, EvaluationMeasures.TPR(Real = tempobs, Predicted = evaluate))
          tempFPR <- c(tempFPR, EvaluationMeasures.FPR(Real = tempobs, Predicted = evaluate))
        }
        tempROC <- cbind("file" = rep(identifier, length(eval_range)), 
                         "outer_fold" = rep(outer, length(eval_range)),
                         "model" = model_name, 
                         "TPR" = tempTPR, 
                         "FPR" = tempFPR,
                         "subset" = rep(subset, length(eval_range)),
                         "inner_fold" = rep(1, length(eval_range)))
        
        test_ROC <- rbind(test_ROC, tempROC)
        
      }
    }
  }
}

validation_ROC <- as.data.frame(test_ROC)
validation_ROC$TPR <- as.numeric(as.character(validation_ROC$TPR))
validation_ROC$FPR <- as.numeric(as.character(validation_ROC$FPR))

innerfolds <- as.numeric(unique(validation_ROC$inner_fold))
outerfolds <- as.numeric(unique(validation_ROC$outer_fold))

validation_results <- as.data.frame(calibrated_test_results_c1)
reformatted_validation_results <- rbind()

per_identifier_plot <- list()
all_FPR_TPR_for_models <- rbind()

for (identifier in identifiers) {
  
  for (subset in subsets) {
    all_FPR_TPR <- rbind()
    all_long_TPR_FPR <- rbind()
    temp_results <- validation_results[which(validation_results$identifier == identifier &
                                               validation_results$subset == subset), ]
    allmodel_names <- unique(validation_ROC[which(validation_ROC$subset == subset), ]$model)  
    for (model_name in allmodel_names) {
      
      TPR <- cbind()
      FPR <- cbind()
      long_TPR_FPR <- rbind()
      for (outer in outerfolds) {
        for (inner in innerfolds) {
          tempROC <- validation_ROC[which(validation_ROC$file == identifier &
                                            validation_ROC$model == model_name &
                                            validation_ROC$outer_fold == outer &
                                            validation_ROC$inner_fold == inner &
                                            validation_ROC$subset == subset),]
          
          TPR <- cbind(TPR, tempROC$TPR)
          FPR <- cbind(FPR, tempROC$FPR)
          long_TPR_FPR <- rbind(long_TPR_FPR, tempROC)
          
        }
      }
      
      FPR <- as.data.frame(FPR)
      FPR$mean <- rowMeans(FPR)
      # FPR$mean <- apply(FPR, 1, mean)
      FPR$sd <- apply(FPR, 1, sd)
      FPR$sem <- FPR$sd/sqrt(folds_outer)
      
      FPR$CI_lower <- FPR$mean + qt((1-0.95)/2, df=folds_outer)*FPR$sem
      FPR$CI_upper <- FPR$mean - qt((1-0.95)/2, df=folds_outer)*FPR$sem
      
      TPR <- as.data.frame(TPR)
      TPR$mean <- rowMeans(TPR)
      # TPR$mean <- apply(TPR, 1, mean)
      TPR$sd <- apply(TPR, 1, sd)
      TPR$sem <- TPR$sd/sqrt(folds_outer)
      
      TPR$CI_lower <- TPR$mean + qt((1-0.95)/2, df=folds_outer)*TPR$sem
      TPR$CI_upper <- TPR$mean - qt((1-0.95)/2, df=folds_outer)*TPR$sem
      
      FPR_TPR <- data.frame(FPR$mean, FPR$CI_lower, FPR$CI_upper, TPR$mean, TPR$CI_lower, TPR$CI_upper)
      FPR_TPR$model <- model_name
      
      long_TPR_FPR$model <- model_name
      long_TPR_FPR$outer_inner <- paste0(long_TPR_FPR$outer_fold, "_", long_TPR_FPR$inner_fold)
      
      
      all_FPR_TPR <- rbind(all_FPR_TPR, FPR_TPR)
      all_long_TPR_FPR <- rbind(all_long_TPR_FPR, long_TPR_FPR)
    }
    
    all_FPR_TPR$subset <- subset
    all_FPR_TPR$identifier <- identifier
    all_FPR_TPR_for_models <- rbind(all_FPR_TPR_for_models, all_FPR_TPR)
    
    ROC_plot <- ggplot(data=all_FPR_TPR, aes(x=FPR.mean, y=TPR.mean, group = model, color = model)) + geom_line(size=0.6, alpha=0.8) +
      geom_ribbon(aes(x=FPR.mean, y=TPR.mean,ymin=TPR.CI_lower, ymax=TPR.CI_upper, xmin = FPR.CI_lower, xmax=FPR.CI_upper, group = model, fill= model, color = model), alpha = 0.2, show.legend = F, linetype=0) +
      geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = 2) + xlab("False Positive Rate") + ylab("True Positive Rate") +
      labs(subtitle = paste0(subset,": 25 outer fold--", identifier, "--", model_name)) + 
      scale_color_discrete(name = "") + theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 11),
        axis.title = element_text(size = 12)
      ) 
    
    ggsave(paste0("eICU_ROC_plot_", identifier, "_", subset, "_", model_name, ".png"), width = 6, height = 5, dpi = 300, units = "in")
    
    
    metrics <- c("Acc", "AUC", "Sens", "Spec")
    
    temp_plots <- list()
    for (metric in metrics) {
      
      temp_plots[[metric]] <- ggplot(temp_results, aes(x = as.numeric(as.character(.data[[metric]])), y = model, group = factor(model), color = factor(model))) + geom_boxplot() + 
        ylab(NULL) + xlab(NULL) + ggtitle(metric) + theme(legend.position = "none") + xlim(0.5, 1)
      
    } 
    
    x <- gridExtra::grid.arrange(temp_plots$Sens, temp_plots$AUC, temp_plots$Spec, temp_plots$Acc, ncol = 2)
    
    per_identifier_plot[[identifier]] <- gridExtra::grid.arrange(ROC_plot, x, ncol = 2)
    
  }
  
}

png(paste0("eICU_Full_Calibrated_Test_ROC_results.png"), width = 8, height = 3, units = "in", res = 300)
gridExtra::grid.arrange(per_identifier_plot[[1]], ncol = 1)
dev.off()


##### test Precision Recall:calibrated#####

savedir <- paste0(dir)
dir.create(savedir)
setwd(savedir)  

test_probs <- calibrated_test_probs_c1
test_PR <- c()
for (identifier in identifiers) {
  for (subset in subsets) {
    allmodel_names <- unique(test_probs[which(test_probs$subset == subset), ]$model)
    for (model_name in allmodel_names) {
      for (outer in unique(test_probs$outer_fold)) {
        
        temp_probs <- test_probs[which(test_probs$identifier == identifier &
                                         test_probs$model == model_name &
                                         test_probs$outer_fold == outer &
                                         test_probs$subset == subset), ]
        temppred <- temp_probs$class_1
        tempobs <- 1 * (as.character(temp_probs$obs) == "class_1")
        eval_range <- seq(0, 1, 0.01) 
        tempP <- c()
        tempR <- c()
        
        for (eval in eval_range) {
          evaluate <- 1 * (temppred > eval)
          tempP <- c(tempP, EvaluationMeasures.Precision(Real = tempobs, Predicted = evaluate))
          tempR <- c(tempR, EvaluationMeasures.Recall(Real = tempobs, Predicted = evaluate))
        }
        tempPR <- cbind("file" = rep(identifier, length(eval_range)), 
                        "outer_fold" = rep(outer, length(eval_range)),
                        "model" = model_name, 
                        "Precision" = tempP, 
                        "Recall" = tempR,
                        "subset" = rep(subset, length(eval_range)),
                        "inner_fold" = rep(1, length(eval_range)))
        
        test_PR <- rbind(test_PR, tempPR)
        
      }
    }
  }
}

validation_PR <- as.data.frame(test_PR)
validation_PR$Precision <- as.numeric(as.character(validation_PR$Precision))
validation_PR$Recall <- as.numeric(as.character(validation_PR$Recall))

innerfolds <- as.numeric(unique(validation_PR$inner_fold))
outerfolds <- as.numeric(unique(validation_PR$outer_fold))

validation_results <- as.data.frame(calibrated_test_results_c1)
reformatted_validation_results <- rbind()

per_identifier_plot <- list()
for (identifier in identifiers) {
  for (subset in subsets) {
    
    temp_results <- validation_results[which(validation_results$identifier == identifier &
                                               validation_results$subset == subset), ]
    allmodel_names <- unique(validation_ROC[which(validation_ROC$subset == subset), ]$model)  
    
    temp_results_reformat <- rbind()
    all_P_R <- rbind()
    all_long_P_R <- rbind()
    for (model_name in allmodel_names) {
      
      Prec <- cbind()
      Rec <- cbind()
      long_P_R <- rbind()
      for (outer in outerfolds) {
        for (inner in innerfolds) {
          tempROC <- validation_PR[which(validation_PR$file == identifier &
                                           validation_PR$model == model_name &
                                           validation_PR$outer_fold == outer &
                                           validation_PR$inner_fold == inner &
                                           validation_PR$subset == subset),]
          
          Prec <- cbind(Prec, tempROC$Precision)
          Rec <- cbind(Rec, tempROC$Recall)
          long_P_R <- rbind(long_P_R, tempROC)
          
        }
      }
      
      Rec <- as.data.frame(Rec)
      Rec <- replace_nan(Rec, 1)
      # Rec$mean <- apply(Rec, 1, mean)
      Rec$mean <- rowMeans(Rec)
      Rec$sd <- apply(Rec, 1, sd)
      Rec$sem <- Rec$sd/sqrt(folds_outer)
      
      Rec$CI_lower <- Rec$mean + qt((1-0.95)/2, df=folds_outer)*Rec$sem
      Rec$CI_upper <- Rec$mean - qt((1-0.95)/2, df=folds_outer)*Rec$sem
      
      Prec <- as.data.frame(Prec)
      Prec <- replace_nan(Prec, 1)
      Prec$mean <- rowMeans(Prec)
      # Prec$mean <- apply(Prec, 1, mean)
      Prec$sd <- apply(Prec, 1, sd)
      Prec$sem <- Prec$sd/sqrt(folds_outer)
      
      Prec$CI_lower <- Prec$mean + qt((1-0.95)/2, df=folds_outer)*Prec$sem
      Prec$CI_upper <- Prec$mean - qt((1-0.95)/2, df=folds_outer)*Prec$sem
      
      P_R <- data.frame(Rec$mean, Rec$CI_lower, Rec$CI_upper, Prec$mean, Prec$CI_lower, Prec$CI_upper)
      P_R$model <- model_name
      
      long_P_R$model <- model_name
      long_P_R$outer_inner <- paste0(long_P_R$outer_fold, "_", long_P_R$inner_fold)
      
      
      all_P_R <- rbind(all_P_R, P_R)
      all_long_P_R <- rbind(all_long_P_R, long_P_R)
    }
    
    all_long_P_R <- replace_nan(all_long_P_R, 1)
    
    balance <- test_class_balance[which(test_class_balance$identifier == identifier &
                                          test_class_balance$subset == subset &
                                          test_class_balance$type == "eICU"),]
    balance$prop1 <- balance$class_0 / (balance$class_0 + balance$class_1) 
    balance$prop2 <- balance$class_1 / (balance$class_0 + balance$class_1) 
    
    lower_class <- mean(balance$prop2)
    
    
    
    ROC_plot <- ggplot(data=all_P_R, aes(x=Rec.mean, y=Prec.mean, group = model, color = model)) + geom_path(size=0.6, alpha=0.8) +
      geom_ribbon(aes(x=Rec.mean, y=Prec.mean,ymin=Prec.CI_lower, ymax=Prec.CI_upper, xmin = Rec.CI_lower, xmax=Rec.CI_upper, group = model, fill= model, color = model), alpha = 0.2, show.legend = F, linetype=0) +
    
      geom_hline(yintercept = lower_class, color = "grey50", linetype = 2) + xlab("Recall") + ylab("Precision") +
      labs(subtitle = paste0(subset, "-25 outer fold--", identifier, "--", model_name))+ 
      scale_color_discrete(name = "") + theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 11),
        axis.title = element_text(size = 12)
      ) +
      ylim(0,1)
    
    ggsave(paste0("eICU_PR_plot_", identifier, "_", subset, ".png"), width = 6, height = 5, dpi = 300, units = "in")
    
    metrics <- c("PR_AUC", "Precision", "Recall", "F1")
    
    temp_plots <- list()
    for (metric in metrics) {
      
      temp_plots[[metric]] <- ggplot(temp_results, aes(x = as.numeric(as.character(.data[[metric]])), y = model, group = factor(model), color = factor(model))) + geom_boxplot() + 
        ylab(NULL) + xlab(NULL) + ggtitle(metric) + theme(legend.position = "none") + xlim(0.5, 1)
      
    }
    
    
    x <- gridExtra::grid.arrange(temp_plots$Precision, temp_plots$PR_AUC, temp_plots$Recall, temp_plots$F1, ncol = 2)
    
    per_identifier_plot[[identifier]] <- gridExtra::grid.arrange(ROC_plot, x, ncol = 2)
    
  }
}

png(paste0("Full_Calibrated_Test_PR_results.png"), width = 8, height = 3, units = "in", res = 300)
# gridExtra::grid.arrange(per_identifier_plot[[1]], per_identifier_plot[[2]], per_identifier_plot[[3]], ncol = 1)
gridExtra::grid.arrange(per_identifier_plot[[1]], ncol = 1)
dev.off()



###### results table.  #####
# summary(calibrated_test_results_c1)

setwd(dir)  

compiled_data_table <- rbind()

for (identifier in identifiers){
  for (subset in subsets) {
    df <- calibrated_test_results_c1[which(calibrated_test_results_c1$identifier == identifier &
                                             calibrated_test_results_c1$subset == subset), ]
    
    model_names <- unique(df$model)
    
    for (model in model_names) {
      df_model <- df[which(df$model == model), ]
      
      #check for columns with NAs
      df_model <- df_model[ , colSums(is.na(df_model)) == 0]
      
      data_cols <- colnames(df_model)
      data_cols <- data_cols[-which(data_cols %in% c("file", "outer_fold", "model", "identifier", "subset"))]
      
      temp_row <- c()
      
      for (colname in data_cols) {
        one_col <- df_model %>% dplyr::select(colname)
        
        decimals <- 2
        
        upper_int <- format(round(CI(one_col[,1], ci = 0.95)[1], digits = decimals), nsmall = decimals)
        lower_int <- format(round(CI(one_col[,1], ci = 0.95)[3], digits = decimals), nsmall = decimals)
        mean <- format(round(CI(one_col[,1], ci = 0.95)[2], digits = decimals), nsmall = decimals)
        
        result_string <- paste0(mean, " (", upper_int, ", ", lower_int, ")")
        
        temp_row <- c(temp_row, result_string)
      }
      
      temp_row <- as.data.frame(t(temp_row), stringsAsFactors = F)
      temp_row$model <- model
      temp_row$identifier <- identifier
      temp_row$subset <- subset
      colnames(temp_row) <- c(data_cols, "model", "identifier", "subset")
      
      compiled_data_table <- rbind.fill(compiled_data_table, temp_row)
      compiled_data_table <- compiled_data_table %>% dplyr::select("model", "identifier", "subset", everything())
      
    }
  }
}

write.csv(compiled_data_table, "eicu_model_results_formatted.csv", row.names = F)

}

#isoreg fit
fit.isoreg <- function(iso, x0) 
{
  o = iso$o
  if (is.null(o)) 
    o = 1:length(x)
  x = iso$x[o]
  y = iso$yf
  ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
  min.x <- min(x)
  max.x <- max(x)
  adjusted.knots <- iso$iKnots[c(1, which(iso$yf[iso$iKnots] > 0))]
  fits = sapply(seq(along = x0), function(i) {
    j = ind[i]
    
    # Handles the case where unseen data is outside range of the training data
    if (is.na(j)) {
      if (x0[i] > max.x) j <- length(x)
      else if (x0[i] < min.x) j <- 1
    }
    
    # Find the upper and lower parts of the step
    upper.step.n <- min(which(adjusted.knots > j))
    upper.step <- adjusted.knots[upper.step.n]
    lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )
    
    # Pefrom a liner interpolation between the start and end of the step
    denom <- x[upper.step] - x[lower.step] 
    denom <- ifelse(denom == 0, 1, denom)
    val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
    
    # Ensure we bound the probabilities to [0, 1]
    val <- ifelse(val > 1, max.x, val)
    val <- ifelse(val < 0, min.x, val)
    val <- ifelse(is.na(val), max.x, val) # Bit of a hack, NA when at right extreme of distribution
    val
  })
  fits
}
