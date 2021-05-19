library(tidyverse)
library(randomForestSRC)
library(ranger)
library(randomForest)
library(caret)
library(MLmetrics)
require(Metrics)
source("~/useful_functions.R")



dir <- "~/RO1_2020/Data/VIF_1_15_21/"
dir.create(dir)
setwd(dir)

#load dataset. 
combined <- read.csv("../lab_mGCS_sofa_pts_imputed_features_excluded_1_11_21.csv")
combined <- combined[, c(1:10)]
colnames(combined)[1] <- "patientunitstayid"
identifier <- "MVtask1"
date <- Sys.Date()

label <- read.csv("../labels_excluded_1_11_21.csv")
colnames(label)[1] <- "patientunitstayid"

label <- label %>% dplyr::select(patientunitstayid, label)

#remove unwanted features and standardize data
combined <- combined[complete.cases(combined), ]
pids <- combined %>% dplyr::select(patientunitstayid)

combined <- BBmisc::normalize(combined[,-1], method = "standardize", range = c(0, 1))
combined <- cbind(pids, combined)


#if labels are separate, load labels

#combine labels if separate
combined <- merge(label, combined, by = "patientunitstayid")

#turn label into binary 0, 1
df <- combined
table(df$label)

df$label <- (df$label == "class1") * 1

#remove patient identifier if patient identifiers exist
# df <- df[,-1]

#split into training and testing for LM 
set.seed(1206)
intrain <- createDataPartition(y = df$label, p = 0.80, list = FALSE)
training <- df[intrain,]
train_pids <- training %>%  dplyr::select(patientunitstayid)
training <- training[,-1]

# table(training$label)
testing <- df[-intrain,]
test_pids <- testing %>%  dplyr::select(patientunitstayid)
testing <- testing[,-1]


#first step to reduce feature number. If too large of a dimension, can't compute VIF. 
#lets first do a correlation analysis to remove highly correlated variables first. 
#
training <- Filter(function(x) sd(x) != 0, training)


correlationMatrix <- cor(training)
print(correlationMatrix)
heatmap(correlationMatrix, keep.dendro = F)
highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.90, exact = T, names = T)
print(highlyCorrelated)

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
vif_features <- merge(label, vif_features, by = "patientunitstayid")

VIF_names <- as.data.frame(colnames(vif_features)[-1])
VIF_names[,1] <- as.character(VIF_names[,1])
VIF_names$order <- 1
VIF_names$category <- "none"
colnames(VIF_names)[1] <- "x"

ordering <- unique(VIF_names$order)


threshold <- 0.5
current_selected <- c()
main_AIC <- rbind()
for (i in 1:length(ordering)) {
  print(i) 
  
  temp_names <- VIF_names[which(VIF_names$order == i), ]
  
  # temp_features <- vif_features %>% select("label", temp_names$x)
  temp_names_upper <- temp_names
  while (nrow(temp_names_upper) != 0) {
    
    temp_names <- temp_names_upper
    
    AIC1 <- rbind()
    while (nrow(temp_names) != 0) {
      print(nrow(temp_names))
      print(temp_names)
      
      if (temp_names$category[1] != "none") {
        variables <- c(temp_names$x[which(temp_names$category == temp_names$category[1])], current_selected)
        temp_features <- vif_features[, which(colnames(vif_features) %in% c("label", variables)), drop = F]
        temp_names <- temp_names[-which(temp_names$x %in% variables), , drop = F]
        
      } else {
        variables <- c(temp_names$x[1], current_selected)
        temp_features <- vif_features[, which(colnames(vif_features) %in% c("label", variables)), drop = F]
        temp_names <- temp_names[-which(temp_names$x %in% variables), , drop = F]
      }
      
      print(variables)
      
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


# 
# write.csv(current_selected, "~/RO1_2020/Data/MV_excluded_2hr_1_12_21//ALL_AIC_reduced_11621.csv", row.names = F)
# write.csv(main_AIC, "~/RO1_2020/Data/MV_excluded_2hr_1_12_21//ALL_AIC_reduced_AIC_values_11621.csv", row.names = F)

plot(seq(1, nrow(main_AIC)), main_AIC$aic, 'l', xlab = "feature #", ylab = "AIC", main = "Forward Selection AIC: threshold > 0.5")
points(seq(1, nrow(main_AIC)), main_AIC$aic)

plot(seq(1, nrow(main_AIC)), main_AIC$auc, 'l', xlab = "feature #", ylab = "AUROC", main = "Forward Selection AUROC: threshold > 0.5")
points(seq(1, nrow(main_AIC)), main_AIC$auc)

plot(seq(1, nrow(main_AIC[-2,])), main_AIC$rel_like[-2], 'l', xlab = "feature #", ylab = "AUROC", main = "Forward Selection rel_like: threshold > 0.5", ylim = c(0.1, 10))
points(seq(1, nrow(main_AIC[-2,])), main_AIC$rel_like[-2])


######## remake feature space
selected_threshold <- 1

selected_main_AIC <- rbind()
var_list <- c()

count = 1
thresh_eval <- TRUE
while (thresh_eval) {
  
  temp <- main_AIC[count, ]
  temp_var <- temp$variables
  temp_var <- strsplit(temp_var, ",")[[1]]
  
  if (temp$rel_like < selected_threshold) {
    thresh_eval <- FALSE
    
  } else {
    selected_main_AIC <- rbind(selected_main_AIC, temp)
    var_list <- unique(c(var_list, temp_var))
    
  }
  
  count = count + 1
}


combined <- read.csv("../lab_mGCS_sofa_pts_imputed_features_excluded_1_11_21.csv")
colnames(combined)[1] <- "patientunitstayid"

combined <- combined[, which(colnames(combined) %in% c("patientunitstayid", var_list))]
# write.csv(combined, "~/RO1_2020/lab_mGCS_sofa_pts_imputed_features_excluded_1_18_21_FS.csv", row.names = F)





