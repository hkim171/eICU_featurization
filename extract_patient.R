#### Usage: extract features from patient table
# Author: Han
# Updated and Formatted for distribution: FEB 2021
# V1.2
# Fixes:
# 2/18/21 - updates to function parameters
#         - usage of code_dir instead of function_dir to load in data.
#         - usage of eICU_dir instead of individual table calls to load in patient and hospital data
#         - integrates optional boolean IF patient and hospital tables are not located in the same directory
#         - renamed make_binary and label_boolean to as_binary and labels_only
# 3/8/21  - Updates to input parameters for simplicity. (less parameters)
#-------------------------------------------------------------------------------------------------------------------------------------#

#### input parameters: ####
# code_dir:                  Main github pull directory - ensures all required scripts are in 
#                            one place - "./eICU_featurization/" 
#                           
# patientunitstayids:        patientunitstayid identifier list or dataframe with patientunitstayid column
# 
# eicu_dir:                  directory pointing to the raw eICU data tables (if all in one place)
#                            extract patient utilizes patient and hospital tables in their original names: 
#                            hospital.csv and patient.csv
#
# patient_table:             [OPTIONAL] if and only if eicu_dir does not contain the patient table
#
# hospital_table:            [OPTIONAL] if and only if eicu_dir does not contain the hospital table
#
# as_binary:                 [OPTIONAL defualt = F] Boolean - if set to TRUE will make categorical variables 
#                            into one-hot-encoded (dummy) binary categories. 
#
# labels_only:               [OPTIONAL defualt = F] Boolean - if set to TRUE, only extracts features relavant 
#                            as labels (discuss with mentors to identify which clinical outcomes are of interest)
#
#-------------------------------------------------------------------------------------------------------------------------------------#

#### Requied files description: ####
# patient table (raw eICU downloaded patient table)
# hospital table (raw eICU downloaded hospital table)
# ApacheDX_dict.xlsx - used for categorizing apache admission dx - this file is still incomplete as there are more input required from clinicians. 
#   Current iteration version is looked at by Dr. Nadar Faraday - July 2020
#   Sending to Dr. Robert Stevens DEC 2020 to finalize and get consensus. 
#-------------------------------------------------------------------------------------------------------------------------------------#

#### preprocessing done: ####
# *note: All preprocessing is done with the entirety of the dataset then subsetted for population of interest. Therefore, 
#        the result of the preprocessing will be exactly the same no matter which subset of the population is selected. 
#
#age "> 89" is replaced with 90
#age < 18 removed
#gender: male female and NA
#admission height: height 3sd from mean are removed. 
#admission weight: implausible weights removed.
#admit time of day are converted to hours from start of day. 
#BMI and all other ideal body weight calculations are derived from weight and height. 
#categorizations of APACHE admission diagnosis. 
#-------------------------------------------------------------------------------------------------------------------------------------#
# 
# #### Code example ####
# library(data.table)
# dir <- "/storage/eICU/" #data_dir
# 
# patient <- fread(paste0(dir, "/patient.csv"))
# pids <- patient$patientunitstayid
# hospital <- fread(paste0(dir, "/hospital.csv"))
# dictionary_path = "/storage/eICU/eICU_feature_extract/ApacheDX_dict.xlsx"
# function_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/"
# 
# patient_labels <- extract_patient(patientunitstayid_list = pids, patient_table = patient, hospital_table = hospital, apachedx_dictionary_path = dictionary_path, as_binary = T, labels_only = T, function_dir = function_dir)
# 
# patient_features <- extract_patient(patientunitstayid_list = pids, patient_table = patient, hospital_table = hospital, apachedx_dictionary_path = dictionary_path, as_binary = T, labels_only = F, function_dir = function_dir)
# 
# code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")
# 
# patient_features <- extract_patient(code_dir = code_dir, eicu_dir = "/storage/eICU/", patientunitstayid_list = pids, as_binary = T, labels_only = F)


####Function

extract_patient <- function(code_dir, 
                            patientunitstayid_list,
                            eicu_dir,
                            patient_table, #optional
                            hospital_table, #optional
                            as_binary, #optional
                            labels_only ) { #optional
  
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
    }
    
    if (!file.exists("ApacheDX_dict.xlsx")) {
      stop("ApacheDX_dict.xlsx is not within the code dir. Make sure code_dir is correct.\n
            If correct, do not rename or move code out of the code_dir and ensure ApacheDX_dict.xlsx\n
            is in the directory")
    } else {
      print("found ApacheDX_dict.xlsx in code dir")
      apachedx_dictionary_path <- paste0(code_dir, "/ApacheDX_dict.xlsx")
    }
    
  }
  
  #check if eicu_dir exists
  if (!dir.exists(eicu_dir)) {
    stop("eicu_dir does not exist - please double check")
  }
  
  #check whether patient and hospital table directories have been specified, if not take eicu_dir and 
  #create directory calls within it for patient.csv and hospital.csv. 
  #error and messages included to inform user. 
  if (missing(patient_table)) {
    message("patient_table not directly specified - assuming patient.csv exists in eicu_dir")
    patient_table <- paste0(eicu_dir, "/patient.csv")
    
    if(!file.exists(patient_table)) {
      stop("patient.csv does not exist in the eicu_dir specified")
    } else {
      patient_table <- fread(patient_table)
    }
  } else {
    message("Using specified patient.csv in patient_table directory")
  }
  
  if (missing(hospital_table)) {
    message("hospital_table not directly specified - assuming hospital.csv exists in eicu_dir")
    hospital_table <- paste0(eicu_dir, "/hospital.csv")
    
    if(!file.exists(hospital_table)) {
      stop("hospital.csv does not exist in the eicu_dir specified")
    } else {
      hospital_table <- fread(hospital_table)
    }
  } else {
    message("Using specified hospital.csv in hospital_table directory")
  }
  
  #set as_binary to FALSE and labels_only to FALSE by default. 
  if (missing(as_binary)) {
    as_binary <- FALSE
  }
  
  if (missing(labels_only)) {
    labels_only <- FALSE
  }

  package_list <- c("tidyverse", "plotly", "chron", "openxlsx")
  load_packages(package_list)
  
  #Will preprocess age to filter by age under 18. 
  #age -------------------------------------------------------------------
  #description: randomly sample from decreasing tail of normal distribution to impute values encoded as >89. 
  set.seed(1206)
  A <- rnorm(1e+05, mean = 0, sd = 1)
  A <- A[A >= 0]
  A <- 10 * A/max(A) + 90
  patient_table$age <- as.character(patient_table$age)
  
  over89 <- patient_table[patient_table$age == "> 89", ]
  under90 <- patient_table[-which(patient_table$age == "> 89"), ]
  
  patient_table$age[patient_table$age == "> 89"] <- as.character(round(sample(x = A, size = nrow(over89), 
                                                                    replace = F)))
  patient_table$age <- as.numeric(patient_table$age)
  
  #main population of interest is adult yet eICU contains < 18. Remove anything under 18
  patient_table <- patient_table[which(patient_table$age >= 18), ]  
  
  if (labels_only == FALSE) {
    
    #will only dplyr::select the following raw data from table. 
    patient_table <- patient_table %>% dplyr::select(patientunitstayid, gender, age, ethnicity, 
                                    hospitalid, wardid, apacheadmissiondx, admissionheight, admissionweight, 
                                    hospitaladmittime24, hospitaladmitoffset, hospitaladmitsource, unittype, 
                                    unitadmittime24, unitadmitsource, unitvisitnumber)
    
    #### each feature will need some preprocessing to proceed. ####
    
    #gender ----------------------------------------------------------------
    #description: remove all non female and male and make binary: male = 1, female = 0
    patient_table$gender <- as.factor(patient_table$gender)
    levels(patient_table$gender) <- c(NA, "Female", "Male", NA, NA)
    patient_table$gender <- as.character(patient_table$gender)
    
    #encode male as 1, female as 0. 
    patient_table$gender <- 1 * (patient_table$gender == "Male")
    
    #ethnicity -------------------------------------------------------------
    patient_table$ethnicity <- as.factor(patient_table$ethnicity)
    levels(patient_table$ethnicity) <- c(NA, "African American", "Asian", "Caucasian", 
                                    "Hispanic", "Native American", NA)
    patient_table$ethnicity <- as.character(patient_table$ethnicity)
    
    #apacheadmissiondx -----------------------------------------------------
    patient_table$apacheadmissiondx <- as.character(patient_table$apacheadmissiondx)
    patient_table$apacheadmissiondx[which(patient_table$apacheadmissiondx == "")] <- NA
    
    ApacheDx_Lookup_Table  = read.xlsx(apachedx_dictionary_path, sheet = 'AllDxLookup')
    
    apache <- patient_table %>% dplyr::select(patientunitstayid, apacheadmissiondx)
    colnames(apache)[2] <- "ApacheDx"
    
    Dx_combined <- merge(apache, ApacheDx_Lookup_Table, by = "ApacheDx", all = T)
    Dx_combined <- Dx_combined[!is.na(Dx_combined$ApacheDx), ]
    
    Dx_combined <- Dx_combined[, -1]
    
    #make apache dx binary

    #requires which_character function - 
    #requires create_binary function - 
    Dx_combined <- create_binary(as.data.frame(Dx_combined))
    
    patient_table <- merge(patient_table, Dx_combined, by = "patientunitstayid", all = T)
    patient_table <- patient_table %>% dplyr::select(-apacheadmissiondx)
    
    #admission_weight ------------------------------------------------------
    #remove weight < 30 and weight > 300
    patient_table$admissionweight[which(patient_table$admissionweight < 30 | patient_table$admissionweight > 
                                     250)] <- NA
    
    #admission_height ------------------------------------------------------
    #remove all outside of 3 sd.
    lower <- mean(patient_table$admissionheight, na.rm = T) - (3 * sd(patient_table$admissionheight, 
                                                                 na.rm = T))
    upper <- mean(patient_table$admissionheight, na.rm = T) + (3 * sd(patient_table$admissionheight, 
                                                                 na.rm = T))
    patient_table$admissionheight[which(patient_table$admissionheight < lower | patient_table$admissionheight > 
                                     upper)] <- NA
    
    #Hospital Admit time ---------------------------------------------------
    #will encode this as hours
    patient_table$hospitaladmittime24 <- as.character(patient_table$hospitaladmittime24)
    
    patient_table$hospitaladmittime24 <- 60 * 24 * as.numeric(chron::times(patient_table$hospitaladmittime24))/60
    
    #Unit Admit time -------------------------------------------------------
    #will encode this as hours
    patient_table$unitadmittime24 <- as.character(patient_table$unitadmittime24)
    
    patient_table$unitadmittime24 <- 60 * 24 * as.numeric(chron::times(patient_table$unitadmittime24))/60
    
    #hospitalid ------------------------------------------------------------
    hospital <- hospital_table
    
    hospital$numbedscategory <- as.factor( hospital$numbedscategory)
    levels(hospital$numbedscategory) <- c(NA, "< 100", ">= 500", "100 - 249", 
                                          "250 - 499")
    hospital$bed_count <- as.character(hospital$numbedscategory)
    
    hospital$region <- as.factor(hospital$region)
    levels(hospital$region) <- c(NA, "Midwest", "Northeast", "South", "West")
    hospital$region <- as.character(hospital$region)
    
    hospital$teachingstatus <- (as.character(hospital$teachingstatus) == "True") * 
      1
    
    hospital <- hospital %>% dplyr::select(hospitalid, teachingstatus, bed_count, region)
    
    patient_table <- merge(patient_table, hospital, by = "hospitalid", all = T)
    
    #hospital_admit_source ------------------------------------------------
    patient_table$hospitaladmitsource <- as.character(patient_table$hospitaladmitsource)
    patient_table$hospitaladmitsource[which(patient_table$hospitaladmitsource == "")] <- NA
    
    #hospitaladmitoffset to hours------------------------------------------
    patient_table$hospitaladmitoffset <- patient_table$hospitaladmitoffset/60
    
    #unit_type ------------------------------------------------------------
    patient_table$unittype <- as.character(patient_table$unittype)
    
    #unit_admit_source ----------------------------------------------------
    patient_table$unitadmitsource <- as.character(patient_table$unitadmitsource)
    patient_table$unitadmitsource[which(patient_table$unitadmitsource == "")] <- NA
    
    
    #### DERIVED VARIABLES ####
    #BMI ------------------------------------------------------------------
    tempdf <- patient_table[which(!is.na(patient_table$admissionheight) & !is.na(patient_table$admissionweight)),]
    tempdf$BMI <- tempdf$admissionweight/(tempdf$admissionheight/100)^2
    
    #Ideal body weight (IBW) ----------------------------------------------
    #many different folumations. will include two. 
    
    #DEVINE 1974 ideal body weight based on height------------------- 
    tempdf$inches <- tempdf$admissionheight * 0.393701
    
    tempdf_male <- tempdf[which(tempdf$gender == 1), ]
    tempdf_female <- tempdf[which(tempdf$gender == 0), ]
    tempdf_NA <- tempdf[which(is.na(tempdf$gender)), ]
    
    tempdf_male$IBW_devine <- 50 + 2.3 * (tempdf_male$inches - (5 * 12))
    tempdf_female$IBW_devine <- 45.5 + 2.3 * (tempdf_female$inches - (5 * 12))
    tempdf_NA$IBW_devine <- NA
    
    tempdf <- rbind(tempdf_male, tempdf_female, tempdf_NA)
    
    #PETERSON 2017 ideal body weight based on height-----------------
    tempdf$IBW_peterson <- 2.2 * tempdf$BMI + 3.5 * tempdf$BMI * (tempdf$admissionheight/100 - 
                                                                    1.5)
    
    #Adjusted body weight (AjBW) ------------------------------------------
    #DEVINE ---------------------------------------------------------
    tempdf$AjBW_devine <- tempdf$IBW_devine + 0.4 * (tempdf$admissionweight - 
                                                       tempdf$IBW_devine)
    
    #PETERSON -------------------------------------------------------
    tempdf$AjBW_peterson <- tempdf$IBW_peterson + 0.4 * (tempdf$admissionweight - 
                                                           tempdf$IBW_peterson)
    
    tempdf <- tempdf %>% dplyr::select(patientunitstayid, BMI, IBW_devine, IBW_peterson, 
                                AjBW_devine, AjBW_peterson)
    
    patient_table <- merge(patient_table, tempdf, by = "patientunitstayid", all = T)
    
    
    #final patient subsetting and extraction ------------------------------
    patient_table <- patient_table[which(patient_table$patientunitstayid %in% patientunitstayid_list), 
                         ]
    
    #binarize
    if (as_binary) {
      patient_table <- create_binary(as.data.frame(patient_table))
    }
    
    return(patient_table)
    
  } else if (labels_only == TRUE) {
    
    #extracting labels ----------------------------------------------------
    patient_table$unit_LOS <- patient_table$unitdischargeoffset/60
    patient_table$hospital_LOS <- (abs(patient_table$hospitaladmitoffset) + abs(patient_table$hospitaldischargeoffset))/60
    
    patient_table <- patient_table %>% dplyr::select(patientunitstayid, hospitaldischargelocation,hospitaldischargestatus, 
                                    hospital_LOS, unitdischargelocation, unitdischargestatus, 
                                    unit_LOS)
    
    #locations----------------------------------------------------
    patient_table$hospitaldischargelocation <- as.character(patient_table$hospitaldischargelocation)
    patient_table$hospitaldischargelocation[which(patient_table$hospitaldischargelocation == "")] <- NA
    
    patient_table$unitdischargelocation <- as.character(patient_table$unitdischargelocation)
    patient_table$unitdischargelocation[which(patient_table$unitdischargelocation == "")] <- NA
    
    #discharge status---------------------------------------------
    patient_table$hospitaldischargestatus <- as.character(patient_table$hospitaldischargestatus)
    patient_table$hospitaldischargestatus[which(patient_table$hospitaldischargestatus ==  "")] <- NA
    
    patient_table$unitdischargestatus <- as.character(patient_table$unitdischargestatus)
    patient_table$unitdischargestatus[which(patient_table$unitdischargestatus == "")] <- NA
    
    
    #final patient subsetting and extraction ------------------------------
    patient_table <- patient_table[which(patient_table$patientunitstayid %in% patientunitstayid_list), ]
    
    #binarize
    if (as_binary) {
      patient_table <- create_binary(as.data.frame(patient_table))
    }
    
    return(patient_table)
  }
  
}

