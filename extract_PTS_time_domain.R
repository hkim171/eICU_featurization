#### Usage: extracttime domain PTS features from PTS data
#### data table includes: Heart rate, pulse ox, respiratory rate, diastolic bp, systolic bp, and mean bp. 
#### Extract 30 features per signal per patient. 
####
####
# Author: Han
# Updated and Formatted for distribution: MAY 2021
# V1.0
# Fixes:
#  
#-------------------------------------------------------------------------------------------------------------------------------------#

#### input parameters: ####
# code_dir:                  Main github pull directory - ensures all required scripts are in 
#                            one place - "./eICU_featurization/" 
#                           
# patientunitstayid_time:    dataframe with patientunitstayid column, lower bound time column, and upper bound time column. 
# 
# eicu_dir:                  directory pointing to the raw eICU data tables- it expects all the csv datafiles to be present in the database. 
#                            This includes hr.csv
#                            Link to google drive folder with PTS csv files required. 
#                            https://drive.google.com/drive/folders/1CfrT0MEcgNk0UFvvFlfBfYAhjvpsRM4d?usp=sharing
#
# lower_bound_hours:         (optional) if lower_bound is specified, the entire PTS dataset will be thresholded by that time in hours.
#                            If not included, the data will be thresholded by the lower bound time column and upper bound time column included
#                            in the patientunitstayid_time dataframe. 
#
# upper_bound_hours:         (optional) if upper_bound is specified, the entire PTS dataset will be thresholded by that time in hours. 
#                            If not included, the data will be thresholded by the lower bound time column and upper bound time column included
#                            in the patientunitstayid_time dataframe. 
#
# PTS_signal_type:           (optional) can specifiy which signal you want to extract time domain features for. 
#
# binwidth:                  (optional) default to 1 hour. bin size determine initial stat extraction.   
#
#-------------------------------------------------------------------------------------------------------------------------------------#

#### Requied files description: ####
# patientunitstayid_time: dataframe with patientunitstayid column, lower bound time column, and upper bound time column
# 
# hr.csv, spo2.csv, resp.csv, dias.csv, sys.csv, mean.csv should be downloaded and put into the eicu_dir specified. 
# these raw files will be read in. 
#
#-------------------------------------------------------------------------------------------------------------------------------------#

# #### Code example ####
# library(data.table)
# code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")
# 
# patient <- fread(paste0(code_dir, "/test_data/test_feature_space.csv"))[, 1, drop = FALSE]
# patient$lower <- 0
# patient$upper <- 6
# 
# just_hr = extract_PTS_time_domain(code_dir = code_dir,
#                         patientunitstayid_time = patient,
#                         eicu_dir = "/storage/eICU/",
#                         lower_bound_hours = 0,
#                         upper_bound_hours = 6,
#                         PTS_signal_type = "hr",
#                         binwidth = 1)
# 
# all_signals = extract_PTS_time_domain(code_dir = code_dir,
#                                patientunitstayid_time = patient,
#                                eicu_dir = "/storage/eICU/",
#                                lower_bound_hours = 0,
#                                upper_bound_hours = 6,
#                                binwidth = 1)
# 
# all_signals_default = extract_PTS_time_domain(code_dir = code_dir,
#                                       patientunitstayid_time = patient,
#                                       eicu_dir = "/storage/eICU/",
#                                       binwidth = 1)

####Function

extract_PTS_time_domain <- function(code_dir, 
                            patientunitstayid_time,
                            eicu_dir,
                            lower_bound_hours,
                            upper_bound_hours,
                            PTS_signal_type,
                            binwidth) {
  
  #directory checks
  if (!dir.exists(code_dir)) {
    stop("the CODE directory is not valid and/or does not exist. please double check")
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
      source("extract_PTS_init.R")
      source("extract_PTS_aggregate.R")
    }
    
  }
  
  #check if eicu_dir exists
  if (!dir.exists(eicu_dir)) {
    stop("eicu_dir does not exist - please double check")
  }
  
  #check whether data table directories have been specified, if not take eicu_dir and check whether files exist. 

  possible_types = c("hr", "sao2", "resp", "dias", "sys", "mean")
  if (!missing(PTS_signal_type)) {
    
    if (!(PTS_signal_type %in% possible_types)) {
      stop("wrong type specified - check documentation before proceeding")
    } else{
      PTS_table <- paste0(eicu_dir, "/", PTS_signal_type, ".csv")
      
      if(!file.exists(PTS_table)) {
        stop(paste0(PTS_signal_type, ".csv does not exist in the eicu_dir specified"))
      } else {
        print(paste0(PTS_signal_type, ".csv exists: ", file.exists(PTS_table)))
      }
    } 
  } else {
    print("extracting all features for all 5 signal types - checking whether csv files for each signal exists ")
    
    check_data = c()
    for (type in possible_types) {

      PTS_table <- paste0(eicu_dir, "/", type, ".csv")
      print(paste0(type, ".csv exists: ", file.exists(PTS_table)))
      check_data = c(check_data, file.exists(PTS_table))
      
    }
    
    if (sum(check_data) != 6) {
      stop("not all csv files are in the eICU directory as required. please check documentation and re-run")
    }
    
  }
  
  print("checking whether patientunitstayid_time dataframe has required columns:")
  
  required_colnames <- c("patientunitstayid", "lower", "upper")
  if (sum(colnames(patientunitstayid_time) %in% required_colnames) != 3) {
    stop("patientunitstayid_time dataframe is missing the correct column naming convention. Please check the documentation")
  } else {
    print("all columns present")
  }
  

  package_list <- c("tidyverse", "chron", "data.table")
  load_packages(package_list)
  
  if (missing(PTS_signal_type)) {
    types <- c("hr", "resp", "sao2", "dias", "sys", "mean")
  } else {
    types <- c(PTS_signal_type)
  }
  
  if (missing(binwidth)) {
    binwidth = 1
  }
  
  pids = patientunitstayid_time$patientunitstayid
  PTS_features <- as.data.frame(pids)
  colnames(PTS_features) <- "patientunitstayid"
  
  for (type in types) {
    print(type)
    
    PTS_table <- fread(paste0(eicu_dir, "/", type, ".csv"))
    PTS_table <- PTS_table[which(PTS_table$patientunitstayid %in% pids), ]
    
    
    if (!missing(lower_bound_hours)){
      PTS_table <- PTS_table[which(PTS_table$offset >= (lower_bound_hours * 60)), ]
    }

    if (!missing(upper_bound_hours)){
      PTS_table <- PTS_table[which(PTS_table$offset <= (upper_bound_hours * 60)), ]
    }
    
    #start extraction
    # print(type)
    pts_init <- rbind()
    
    print(paste0(type, ": Time-domain - initial statistic extraction"))
    for (i in 1:length(pids)) {
      print(paste0(type, "(", which(types %in% type), "/", length(types), ") -- init -- ", i, "/", length(pids)))
      pid <- patientunitstayid_time[i, ]
      pts_temp <- extract_PTS_init_stats(patientunitstayids = pid$patientunitstayid, table_df = PTS_table, pts_name = type, lower_minute = (pid$lower * 60), upper_minute = (pid$upper * 60), binwidth = binwidth, per_pid_boolean = T)
      pts_init <- rbind(pts_init, pts_temp)
      
    }
    
    pts_init <- as.data.frame(pts_init)
    
    #remove those that have < 2 rows of data to analyze. 
    pts_temp <- pts_init[complete.cases(pts_init), ]
    pts_temp <- pts_temp %>% group_by(patientunitstayid) %>% dplyr::summarise(count = n())
    pts_temp <- pts_temp[which(pts_temp$count >= 2), ]
    hr_pid_2 <- pts_temp$patientunitstayid
    pts_init <- merge(pts_init, pts_temp, by = "patientunitstayid", all = T)
    pts_init$count[is.na(pts_init$count)] <- 0
    
    pts_init <- pts_init[which(pts_init$patientunitstayid %in% hr_pid_2), ]
    
    print(paste0(type, ": Time-domain - aggregate statistic extraction"))
    pts_agg <- extract_PTS_aggregate(init_stat_df = pts_init, pids = pids, identifier = type)
    tmp_df <- as.data.frame(pids)
    colnames(tmp_df) <- "patientunitstayid"
    pts_agg <- merge(tmp_df, pts_agg, by = "patientunitstayid", all = T) %>% distinct(.)
    
    PTS_features <- merge(PTS_features, pts_agg, by = "patientunitstayid", all = T) %>% distinct(.)
  }
  
  
  return(PTS_features)
}
