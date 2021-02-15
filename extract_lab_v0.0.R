# extract all features from LAB table 1

#example
# patient <- fread("~/RO1_2020/ICP/propofol/ICP_meta_table.csv")
# patient <- patient %>% dplyr::select(patientunitstayid, drug_offset)
# colnames(patient)[2] <- "offset" #diagosis offset
# 
# #we want to only extract data from before the offset for each patient.
# 
# patient <- patient %>% group_by(patientunitstayid) %>% slice(which.min(offset))
# 
# pids <- patient$patientunitstayid
# 
# lab_df <- fread("/storage/eICU/lab.csv")
# lab_df <- lab_df[which(lab_df$patientunitstayid %in% unique(pids)), ]
# 
# lab_final <- rbind()
# for (i in 1:length(pids)) {
#   print(paste0("--------------", i, "--------------"))
# 
#   tempid <- pids[i]
#   temp_row <- patient[patient$patientunitstayid == tempid, ]
#   upper_hr <- temp_row$offset / 60
#   lower_hr <- upper_hr - 6
# 
#   temp_lab <- extract_lab(patientunitstayids = tempid, table_df = lab_df, lower_hr = lower_hr, upper_hr = upper_hr, per_pid_boolean = T, pre_stats = T, extras = T)
# 
#   lab_final <- rbind.fill(lab_final, temp_lab)
# }


extract_lab <- function(patientunitstayids, table_df, lower_hr, upper_hr, per_pid_boolean, pre_stats, extras) {
  source("./useful_functions.R")
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
  
  if (missing(extras)) {
    extras == FALSE
  }

  if (per_pid_boolean == TRUE) {
    patientunitstayids <- c(patientunitstayids)
  }
  
  
  
  # filter for pids and time frame
  table_df <- table_df[which(table_df$patientunitstayid %in% patientunitstayids), ] %>% droplevels(.)
  
  if (pre_stats == TRUE) {
    pre_stats_df <- table_df[which(table_df$labresultrevisedoffset > -Inf & table_df$labresultrevisedoffset <= (lower_hr * 60)), ] %>% droplevels(.)
    pre_stats_df$labname <- as.character(pre_stats_df$labname)
  }
  
  table_df <- table_df[which(table_df$labresultrevisedoffset > (lower_hr * 60) & table_df$labresultrevisedoffset < (upper_hr * 60)), ] %>% droplevels(.)
  
  table_df$labname <- as.character(table_df$labname)
  
  varnames <- as.character(unique(table_df$labname))
  
  lab_features <- as.data.frame(patientunitstayids)
  colnames(lab_features) <- c("patientunitstayid")
  id <- "patientunitstayid"
  
  #fix some of the varnames for uniformity by removing spaces, numbers, dots, and other unsavory variable names.
  table_df$labname <- as.character(table_df$labname)
  table_df$labname <- gsub("-", "", table_df$labname)
  table_df$labname <- gsub(" ", "_", table_df$labname)
  table_df$labname <- gsub("__", "_", table_df$labname)
  table_df$labname <- gsub("\\(|\\)", "", table_df$labname)
  table_df$labname <- gsub("\\%|\\)", "percent", table_df$labname)
  table_df$labname <- gsub("\\.", "", table_df$labname)
  
  table_df$labname[which(table_df$labname == "bedside_glucose")] <- "glucose"
  table_df$labname[which(table_df$labname == "HCO3")] <- "bicarbonate"

  
  varnames <- as.character(unique(table_df$labname))
  # View(varnames)
  
  for (i in 1:length(varnames)) {
    labname <- as.character(varnames[i])
    # print(i)

    #find pre-stats ----------------------------------------------------------------------------------------------
    if (pre_stats == TRUE) {
      idx <- which(pre_stats_df$labname == labname)
      temppre <- pre_stats_df[idx, ] %>% droplevels(.)
      premeanperpatient <- temppre %>% group_by(patientunitstayid) %>% dplyr::summarise(temppre_mean = mean(labresult, 
                                                                                                           na.rm = T))
      colnames(premeanperpatient) <- c(id, paste0(labname, "_pre_mean"))
    }
    
    #find simple statistics --------------------------------------------------------------------------------------
    idx <- which(table_df$labname == labname)
    templabs <- table_df[idx, ] %>% droplevels(.)
    templabs <- templabs[complete.cases(templabs), ]
    #mean
    meanperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::summarise(templabs_mean = mean(labresult, 
                                                                                                         na.rm = T))
    colnames(meanperpatient) <- c(id, paste0(labname, "_mean"))
    
    #sd
    sdperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::summarise(templabs_sd = sd(labresult, 
                                                                                                   na.rm = T))
    colnames(sdperpatient) <- c(id, paste0(labname, "_sd"))
    
    #min
    minperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::summarise(templabs_min = min(labresult, 
                                                                                                      na.rm = T))
    colnames(minperpatient) <- c(id, paste0(labname, "_min"))
    
    #max
    maxperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::summarise(templabs_max = max(labresult, 
                                                                                                      na.rm = T))
    colnames(maxperpatient) <- c(id, paste0(labname, "_max"))
    
    ### extras
    if (extras == T) {
      
      #first
      firstperpatient <- templabs %>% group_by(patientunitstayid) %>% slice(which.min(labresultrevisedoffset))
      firstperpatient <- firstperpatient %>% dplyr::select(patientunitstayid, labresult)
      colnames(firstperpatient) <- c(id, paste0(labname, "_first"))
      
      #last
      lastperpatient <- templabs %>% group_by(patientunitstayid) %>% slice(which.max(labresultrevisedoffset))
      lastperpatient <- lastperpatient %>% dplyr::select(patientunitstayid, labresult)
      colnames(lastperpatient) <- c(id, paste0(labname, "_last"))
      
    }
    
    #combine data and remove correlated features
    combined <- cbind(meanperpatient, minperpatient[,2], maxperpatient[,2], sdperpatient[,2])
    
    if (extras == T) {
      combined <- cbind(combined, firstperpatient[,2], lastperpatient[,2])
    }
    
    #REMOVED
    # #zeros replace all NAs that may come up due to not enough data points for sd. This is done to ensure correlation is still calculated despite there being NAs. Also, by replacing with ZERO it considers the correlation to every point and there is no need to do pairwise complete correlation, which then reduces samples down to only the complete cases. 
    # 
    # tmpcombined <- combined
    # tmpcombined[is.na(tmpcombined)] <- 0
    # 
    # #remove zero variance variables (same number = not predictive)
    # if (nrow(combined) > 1) {
    #   zero_var <- apply(tmpcombined, 2, function(x) length(unique(x)) == 1)
    #   combined <- combined[, !zero_var]
    #   #reset tmpcombined
    #   tmpcombined <- combined
    #   tmpcombined[is.na(tmpcombined)] <- 0
    # }
    # 
    # tmp_cor <- cor(tmpcombined, method = "kendall")
    # tmp_cor_varnames <- findCorrelation(tmp_cor, cutoff = 0.9, names = T)
    # if (length(tmp_cor_varnames != 0)) {
    #   combined <- combined[ , -which(names(combined) %in% tmp_cor_varnames)]
    # }
    
    if (pre_stats == TRUE) {
      combined <- merge(combined, premeanperpatient, by = "patientunitstayid", all = T)
    }
    
    lab_features <- merge(lab_features, combined, by = "patientunitstayid", all = T)
    lab_features <- replace_infinite(lab_features, NA)
    lab_features <- replace_nan(lab_features, NA)
  }
  
  return(lab_features)
  
}
