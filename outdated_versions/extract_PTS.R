##@#### OUTDATED


#EXTRACT features from all PTS data. 
#1. Simple statistics 
#2. Aggregate statistics
#3. PCA components

#author:Han
#date: 3 10 2020

# pts_name = "hr"
# lower_hr <- 24
# upper_hr <- 48
# per_pid_boolean = F
# table_df <- hr
# patientunitsatyids = patientunitstayids

extract_PTS <- function(patientunitstayids, table_df, pts_name, lower_hr, upper_hr, per_pid_boolean) {
  source("~/useful_functions.R")
  require(tidyverse)
  require(plotly)
  library(caret)
  library(doParallel)
  
  if (missing(lower_hr)) {
    lower_hr <- -Inf
  }
  
  if (missing(upper_hr)) {
    upper_hr <- Inf
  }
  
  if (missing(per_pid_boolean)) {
    per_pid_boolean == FALSE
  }
  
  if (per_pid_boolean == TRUE) {
    patientunitstayids <- c(patientunitstayids)
  }
  
  #filter for PTS that meets full length criteria
  #
  #New criteria: ONLY PTS with enough data will be eligible for feature extraction since many features now depend on 
  #the fact that here needs to be enough data. 
  maxoffsets <- table_df %>% group_by(patientunitstayid) %>% distinct(.) %>% 
    dplyr::summarise(max_ = max(offset, na.rm = T))
  max_criteria_pids <- maxoffsets$patientunitstayid[which(maxoffsets$max_ >= (24 * 60))]
  table_df <- table_df[which(table_df$patientunitstayid %in% max_criteria_pids), ] %>% droplevels(.) %>% distinct(.)
  
  # filter for pids and time frame
  table_df <- table_df[which(table_df$patientunitstayid %in% patientunitstayids), ] %>% droplevels(.) %>% distinct(.)
  table_df <- table_df[which(table_df$offset > (lower_hr * 60) & table_df$offset < (upper_hr * 60)), ] %>% droplevels(.) %>% distinct(.)
  
  features <- as.data.frame(patientunitstayids)
  colnames(features) <- c("patientunitstayid")
  id <- "patientunitstayid"
  
  
  ### extract features---------------------------------------------------------------------------------------
  
  
  #simple features--------------------------------------------------------
  meanperpatient = table_df %>% group_by(patientunitstayid) %>% distinct(.) %>% 
    dplyr::summarise(mean_ = mean(value, na.rm = T))
  
  sdperpatient = table_df %>% group_by(patientunitstayid) %>% distinct(.) %>% 
    dplyr::summarise(sd_ = sd(value, na.rm = T))
  
  minperpatient = table_df %>% group_by(patientunitstayid) %>% distinct(.) %>% 
    dplyr::summarise(min_ = min(value, na.rm = T))
  
  maxperpatient = table_df %>% group_by(patientunitstayid) %>% distinct(.) %>% 
    dplyr::summarise(max_ = max(value, na.rm = T))
  
  #autocorrelation, decomposed variances-----------------------------------
  Upid <- unique(table_df$patientunitstayid)
  autocorpid <- c()
  autocor <- c()
  var_seasonal <- c()
  var_trend <- c()
  var_noise <- c()
  
  for (i in 1:length(Upid)){
    tempPid <- Upid[i]
    tempdf <- table_df[table_df$patientunitstayid == tempPid,] %>% distinct(.)
    tempacf <- acf(tempdf$value, lag.max = 24, type = "correlation", na.action = na.pass, plot = F)$acf
    tempacf <- tempacf[!is.na(tempacf)]
    autocorpid <- c(autocorpid, tempPid)
    autocor <- c(autocor, sum(tempacf))
    print(i)
    
    ts_tempdf = ts(tempdf$value, frequency = 24)
    
    if (length(ts_tempdf) > 47) {
    decompose_df = decompose(ts_tempdf, "additive")
    var_seasonal <- c(var_seasonal, var(decompose_df$seasonal, na.rm = T))
    var_trend <- c(var_trend, var(decompose_df$trend, na.rm = T))
    var_noise <- c(var_noise, var(decompose_df$random, na.rm = T))
    } else {
      var_seasonal <- c(var_seasonal,NA)
      var_trend <- c(var_trend, NA)
      var_noise <- c(var_noise, NA)
    }

  }
  
  autocor <- data.frame(autocorpid, autocor, var_seasonal, var_trend, var_noise)
  colnames(autocor) <- c("patientunitstayid", "autocor", "var_seasonal", "var_trend", "var_noise")
  
  stats <- cbind(meanperpatient, minperpatient[,2], maxperpatient[,2], sdperpatient[,2], autocor[,c(2:5)])
  
  varnames <- c("patientunitstayid", "mean", "min", "max", "sd","autocor", "var_seasonal", "var_trend", "var_noise")
  varnames <- c("patientunitstayid", paste0(pts_name, "_", varnames[-1], sep = ""))
  
  colnames(stats) <- varnames
  
  features <- merge(features, stats, by = "patientunitstayid", all = T)
  

  #Windowed features--------------------------------------------------------
  timeframe <- 2 * 60
  
  timeframe_df <- data.frame(matrix(ncol = 9, nrow = length(Upid)))
  varnames <- c("patientunitstayid", "mean_mean", "sd_mean", "mean_min", "mean_max", 
                "mean_delta", "min_delta", "max_delta", "sd_delta")
  colnames(timeframe_df) <- varnames

  
  for (i in 1:length(Upid)) {
    tempPid <- Upid[i]
    tempdf <- table_df[table_df$patientunitstayid == tempPid,] 
    
    maxtime <- max(tempdf$offset)
    
    tempmean <- c()
    tempmin <- c()
    tempmax <- c()
    tempdelta <- c()
    for (j in 1: ceiling(maxtime/timeframe)){
      tempsubset <- tempdf[tempdf$offset >= (j * timeframe - timeframe) & tempdf$offset < (j*timeframe), ] 
      tempmean <- c(tempmean, mean(tempsubset$value, na.rm = T))
      tempmin <- c(tempmin, min(tempsubset$value, na.rm = T))
      tempmax <- c(tempmax, max(tempsubset$value, na.rm = T))
      
      temp <- tail(tempsubset$value, -1) - head(tempsubset$value, -1)
      tempdelta <- c(tempdelta, sum(abs(temp[is.finite(temp)])))
      print(paste0("--",j))
    }
    
    timeframe_df$patientunitstayid[i] <- tempPid
    timeframe_df$mean_mean[i] <- mean(tempmean[is.finite(tempmean)], na.rm = T, )
    timeframe_df$sd_mean[i] <- sd(tempmean[is.finite(tempmean)], na.rm = T)
    timeframe_df$mean_min[i] <- mean(tempmin[is.finite(tempmin)], na.rm = T)
    timeframe_df$mean_max[i] <- mean(tempmax[is.finite(tempmax)], na.rm = T)
    timeframe_df$mean_delta[i] <- mean(tempdelta[is.finite(tempdelta)], na.rm = T)
    timeframe_df$sd_delta[i] <- sd(tempdelta[is.finite(tempdelta)], na.rm = T)
    timeframe_df$min_delta[i] <- min(tempdelta[is.finite(tempdelta)], na.rm = T)
    timeframe_df$max_delta[i] <- max(tempdelta[is.finite(tempdelta)], na.rm = T)
    print(i)
  }

  varnames <- c("patientunitstayid", paste0(pts_name, "_", varnames[-1], sep = ""))
  colnames(timeframe_df) <- varnames
  
  features <- merge(features, timeframe_df,  by = "patientunitstayid", all = T)
  # saveRDS(hr_features, "/home/han/TBI/dec_19/PTS/hr_features.Rds")
  
  
  # #PCA_features ----------------------------------------------------------------------------
  # library(factoextra)
  # library(FactoMineR)
  # 
  # pca_df <- data.frame(matrix(ncol = 288, nrow = length(Upid)))
  # for (i in 1:length(Upid)) {
  #   tempPid <- Upid[i]
  #   tempdf <- table_df[table_df$patientunitstayid == tempPid,] 
  #   pca_df[i, 1] <- tempPid
  #   
  #   for (j in 2:288) {
  #     pca_df[i,j] <- tempdf$value[j - 1]
  #     print(paste0("--",j))
  #   }
  #   print(i)
  # }
  # 
  # rownames(pca_df) <- pca_df[,1]
  # 
  # df.pca <- prcomp(pca_df[, -1], scale = FALSE)
  # 
  # # #visualize the pca data. 
  # # print(df.pca)
  # # summary(df.pca)
  # # pcaCharts(df.pca)
  # ind <- get_pca_ind(df.pca)
  # pca_10 <- as.data.frame(ind$coord[, c(1:10)])
  # 
  # pca_features <- cbind(pca_df[,1], pca_10)
  # 
  # varnames <- paste0(pts_name, "_pca_", as.character(seq(from = 1, to = 10, by = 1)))
  # 
  # colnames(pca_features) <- c("patientunitstayid", varnames)
  # 
  # features <- merge(features, pca_features,  by = "patientunitstayid", all = T)
  
  
  #POLVAR, PermEn-------------------------------------------------------------------------
  #uses the polvar function in useful_function.R
  #uses the PermEn function in useful_function.R
  
  
  P <- data.frame(matrix(ncol = 6, nrow = length(Upid)))
  D2 <- data.frame(matrix(ncol = 9, nrow = length(Upid)))
  D4 <- data.frame(matrix(ncol = 7, nrow = length(Upid)))
  for (i in 1:length(Upid)) {
    tempPid <- Upid[i]
    tempdf <- table_df[table_df$patientunitstayid == tempPid,] 
    
    #extracting those that were predictive in cardiac arrest study. 
    p55 <- polvar(tempdf$value, 5, 5)
    p54 <- polvar(tempdf$value, 5, 4)
    p53 <- polvar(tempdf$value, 5, 3)
    pe42 <- permEN(tempdf$value, 4, 2)
    pe52 <- permEN(tempdf$value, 5, 2)
    
    p <- c(tempPid, p55, p54, p53, pe42, pe52)
    P[i, ] <- p
    
    db2 <- as.data.frame(t(c(tempPid, db(tempdf$value, "haar"))))
    D2[i, ] <- db2
    
    db4 <- as.data.frame(t(c(tempPid, db(tempdf$value, "d4"))))
    D4[i, ] <- db4
    
    print(i)
  }
  
  Pvarnames <- c("patientunitstayid", "polvar5,5", "polvar5,4", "polvar5,3", "permEn4,2", "permEn5,2")
  Pvarnames <- c("patientunitstayid", paste0(pts_name, "_", Pvarnames[-1], sep = ""))
  colnames(P) <- Pvarnames
  
  D2varnames <-  c("patientunitstayid", "db2_W1", "db2_W2", "db2_W3", "db2_W4", "db2_W5", "db2_W6", "db2_W7", "db2_W8" )
  D2varnames <- c("patientunitstayid", paste0(pts_name, "_", D2varnames[-1], sep = ""))
  colnames(D2) <- D2varnames
  
  D4varnames <- c("patientunitstayid", "db4_W1", "db4_W2", "db4_W3", "db4_W4", "db4_W5", "db4_W6")
  D4varnames <- c("patientunitstayid", paste0(pts_name, "_", D4varnames[-1], sep = ""))
  colnames(D4) <- D4varnames
  
  features <- merge(features, P,  by = "patientunitstayid", all = T)
  features <- merge(features, D2,  by = "patientunitstayid", all = T)
  features <- merge(features, D4,  by = "patientunitstayid", all = T)
  
  
  
  return(features)
}


#how to run and things i ran using function above. 
# hr <- readRDS("~/eICU/imputed_3_9_20/hr_filled.Rds")
# hr <- hr %>% distinct(.)
# 
# table_df <- hr
# lower_hr <- 0
# upper_hr <- 24
# per_pid_boolean <- FALSE
# pre_stats <- FALSE
# patientunitstayids <- readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1]
# pts_name <- 'hr'
# 
# hr_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/hr_filled.Rds"), pts_name = "hr", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(hr_features, "~/eICU/imputed_3_9_20/hr_features.Rds")
# 
# sao2_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/sao2_filled.Rds"), pts_name = "sao2", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(sao2_features, "~/eICU/imputed_3_9_20/sao2_features.Rds")
# 
# resp_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/resp_filled.Rds"), pts_name = "resp", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(resp_features, "~/eICU/imputed_3_9_20/resp_features.Rds")
# 
# idias_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/idias_filled.Rds"), pts_name = "idias", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(idias_features, "~/eICU/imputed_3_9_20/idias_features.Rds")
# 
# isys_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/isys_filled.Rds"), pts_name = "isys", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(isys_features, "~/eICU/imputed_3_9_20/isys_features.Rds")
# 
# imean_features <- extract_PTS(patientunitstayids = readRDS("~/TBI/dec_19/pids/pid_TBI_4450.Rds")[,1], table_df = readRDS("~/eICU/imputed_3_9_20/imean_filled.Rds"), pts_name = "imean", lower_hr = 0, upper_hr = 24, per_pid_boolean = FALSE)
# 
# saveRDS(imean_features, "~/eICU/imputed_3_9_20/imean_features.Rds")
