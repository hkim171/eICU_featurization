######  Statistical Lab Features from LAB table. #####

#V0.0 - 2/5/2021 - uploaded without final edits. 
#v0.1 - 3/5/2021 - Edited to contain proper parameters and functionality. 
#v0.2 - 3/11/2021 - wrapper that wraps in for loop for different start and stop times. 

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
#   temp_lab <- extract_lab(patientunitstayids = tempid, table_df = lab_df, lower_hr = lower_hr, upper_hr = upper_hr, per_pid = T, pre_stats = T, extras = T)
# 
#   lab_final <- rbind.fill(lab_final, temp_lab)
# }


######
#Dataframe format for patientunitstayid_dataframe:

#column 1: patientunitstayid
#column 2: start_obs_window in hours
#column 3: stop_obs_window in hours

#Even if they are from the same window, requires the same data format. 
#
#
#

#Example taken from pids from test_data
#
code_dir <- "/storage/eICU/eICU_feature_extract/eICU_featurization/"
data_dir <- paste0(code_dir, "/test_data/")

patientunitstayid_dataframe <- fread(paste0(data_dir, "/test_feature_space.csv"))[,1]
patientunitstayid_dataframe$start <- sample(seq(from = 0, to = 24, by = 0.01), size = nrow(patientunitstayid_dataframe), replace = T)
patientunitstayid_dataframe$end <- patientunitstayid_dataframe$start + 24

eicu_dir <- "/storage/eICU/"
per_pid <- TRUE


extract_lab <- function(code_dir, 
                        patientunitstayid_dataframe,
                        eicu_dir,
                        lab_table, #optional
                        per_pid #T/F
                        ) {
  
  #directory checks
  if (!dir.exists(code_dir)) {
    stop("the CODE directory is not valid and/or does not exist. please double check")
  }
  
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
      source("extract_lab.R")
    }
  }
  
  #check if eicu_dir exists
  if (!dir.exists(eicu_dir)) {
    stop("eicu_dir does not exist - please double check")
  }
  
  #check whether lab table directories have been specified, if not take eicu_dir and 
  #create directory calls within it for lab.csv and/or lab.Rds 
  #error and messages included to inform user. 
  if (missing(lab_table)) {
    message("lab_table not directly specified - assuming patient.csv exists in eicu_dir")
    lab_table_rds <- paste0(eicu_dir, "/lab.Rds")
    lab_table_csv <- paste0(eicu_dir, "/lab.csv")
    
    if(!file.exists(lab_table_rds)) {
      message("lab.Rds does not exist in the eicu_dir specified, will attempt to read in lab.csv if it exists.")
      
      if(!file.exists(lab_table_csv)) {
        stop("lab.csv and lab.Rds does not exist in the eicu_dir specified, please check again")
      } else {
        message("loading lab.csv")
        lab_table <- fread(lab_table_csv)
      }
      
    } else {
      message("loading lab.Rds")
      lab_table <- readRDS(lab_table_rds)
    }
    
  } else {
    message("Using specified lab file in lab_table directory")
  }
  
  #set as_binary to FALSE and labels_only to FALSE by default. 
  if (missing(per_pid)) {
    stop("Must specify whether to consider a per_pid = TRUE (different start and stop offset per pid)
         Or, per_pid = FALSE (the same start and stop offset per pid)")
  }
  
  package_list <- c("tidyverse", "tidyverse", "data.table", "doParallel")
  load_packages(package_list)
  
  
  start_unique_count <- length(unlist(unique(patientunitstayid_dataframe[,2])))
  stop_unique_count <- length(unlist(unique(patientunitstayid_dataframe[,3])))
  
  unique_counts <- start_unique_count + stop_unique_count
  
  if (unique_counts != 0 & per_pid) {
    message("SANITY CHECK: starts and stops do differ per pid")
  } else{
    stop("starts and stops do not differ per pid, or per_pid was set to F yet starts ands stops are different per pid.")
  }
  
  pids <- unlist(unique(patientunitstayid_dataframe[,1]))
  lab_table <- lab_table[which(lab_table$patientunitstayid %in% pids), ]
  lab_table <- lab_table %>% dplyr::select(patientunitstayid, labresultrevisedoffset, labname, labresult)

  lab_final <- rbind()
  for (i in 1:length(pids)) {
    print(paste0("--------------", i, "--------------"))

    tempid <- pids[i]
    temp_row <- patientunitstayid_dataframe[patientunitstayid_dataframe$patientunitstayid == tempid, ]
    lower_hr <- temp_row[[2]]
    upper_hr <- temp_row[[3]]


    lab_table_temp <- lab_table[which(lab_table$patientunitstayid == tempid), ]
    temp_lab <- extract_lab_base_function(patientunitstayids = tempid, table_df = lab_table_temp, lower_hr = lower_hr, upper_hr = upper_hr, per_pid = T, pre_stats = T, extras = T)

    lab_final <- rbind.fill(lab_final, temp_lab)
  }
  
  
  
  
}


extract_lab_base_function <- function(patientunitstayids, 
                                      table_df, 
                                      lower_hr, 
                                      upper_hr, 
                                      per_pid, 
                                      pre_stats, 
                                      extras) {
  
 #set proper directory then make sure all random forest trained files exist. 
  
  if (missing(lower_hr)) {
    lower_hr <- -Inf
  }
  
  if (missing(upper_hr)) {
    upper_hr <- Inf
  }
  
  if (missing(extras)) {
    extras == FALSE
  }

  if (per_pid == TRUE) {
    patientunitstayids <- c(patientunitstayids)
  }
  
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
  

  premeanperpatient_final <- cbind()
  if (pre_stats == TRUE) {
    varnames <- as.character(unique(pre_stats_df$labname))
    # View(varnames)
    
    for (i in 1:length(varnames)) {
      labname <- as.character(varnames[i])
      # print(i)
      
      #find pre-stats ----------------------------------------------------------------------------------------------
      idx <- which(pre_stats_df$labname == labname)
      temppre <- pre_stats_df[idx, ] %>% droplevels(.)
      premeanperpatient <- temppre %>% group_by(patientunitstayid) %>% dplyr::summarise(temppre_mean = mean(labresult, 
                                                                                                            na.rm = T))
      colnames(premeanperpatient) <- c(id, paste0(labname, "_pre_mean"))
    }
  }

  
  varnames <- as.character(unique(table_df$labname))
  
  for (i in 1:length(varnames)) {
    labname <- as.character(varnames[i])
    # print(i)
    
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
      firstperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::slice(which.min(labresultrevisedoffset))
      firstperpatient <- firstperpatient %>% dplyr::select(patientunitstayid, labresult)
      colnames(firstperpatient) <- c(id, paste0(labname, "_first"))
      
      #last
      lastperpatient <- templabs %>% group_by(patientunitstayid) %>% dplyr::slice(which.max(labresultrevisedoffset))
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
  
    
    lab_features <- merge(lab_features, combined, by = "patientunitstayid", all = T)
    lab_features <- replace_infinite(lab_features, NA)
    lab_features <- replace_nan(lab_features, NA)
  }
  
  return(lab_features)
  
}
