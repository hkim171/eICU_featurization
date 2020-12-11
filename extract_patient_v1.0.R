#### Usage: extract features from patient table
# Author: Han
# Updated and Formatted for distribution: DEC 2020
# V1.0
# Fixes:
#-------------------------------------------------------------------------------------------------------------------------------------#

#### input parameters: ####
#patientunitstayid_list:    int/num list - A R list of patientunitstayids as numerics
#patient_table:             data.table/data.frame - RAW eICU patient table
#hospital_table:            data.table/data.frame - RAW eICU hospital table
#apachedx_dictionary_path:  STR - Path to an excel file required to get apachedx information.
#make_binary:               Boolean - if set to TRUE will make categorical variables into one-hot-encoded (dummy) binary categories. 
#label_boolean:             Boolean - if set to TRUE, only extracts features relavant as modeling end points
#
# *note: tables are loaded in prior to input into function. Not a real issue with patient function
#        but other larger tables require this so tables are not repeated called when function it called
#        in a loop. 
#
#-------------------------------------------------------------------------------------------------------------------------------------#

#### files required: ####
#patient table
#hospital table
#ApacheDX_dict.xlsx - used for categorizing apache admission dx - this file is still incomplete as there are more input required from clinicians. 
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

#### Code example is provide below. select lines 46-58 then click ctrl+shift+c to uncomment. ####

# dir <- "/storage/eICU/" #data_dir
# 
# # source <- source("/storage/eICU/eICU_feature_extract/extract_patient_v2.R") #this is an example. you dont want to source the file you are
# #running an example in - do not uncomment.
# 
# patient <- fread(paste0(dir, "/patient.csv"))
# pids <- patient$patientunitstayid
# hospital <- fread(paste0(dir, "/hospital.csv"))
# dictionary_path = "/storage/eICU/eICU_feature_extract/ApacheDX_dict.xlsx"
# 
# patient_labels <- extract_patient(patientunitstayid_list = pids, patient_table = patient, hospital_table = hospital, apachedx_dictionary_path = dictionary_path, make_binary = T, label_boolean = T)
# 
# patient_features <- extract_patient(patientunitstayid_list = pids, patient_table = patient, hospital_table = hospital, apachedx_dictionary_path = dictionary_path, make_binary = T, label_boolean = F)


####Function

extract_patient <- function(patientunitstayid_list, patient_table, hospital_table, apachedx_dictionary_path, make_binary, label_boolean) {
  require(tidyverse)
  require(plotly)
  require(chron)
  require(openxlsx)
  
  if (missing(make_binary)) {
    make_binary <- FALSE
  }
  
  if (missing(label_boolean)) {
    label_boolean <- FALSE
  }
  
  if (missing(apachedx_dictionary_path) & label_boolean) {
    apachedx_dictionary_path <- NA
  } else if (missing(apachedx_dictionary_path)) {
    stop("missing path to apachedx_dictionary_excel_file")
  }
  
  
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
  
  if (label_boolean == FALSE) {
    
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
    
    #gives column name of features that are characters
    which_character <- function(df) {
      df <- as.data.frame(df, stringasfactor = F)
      charvars <- c()
      varnames <- colnames(df)
      for (vars in varnames) {
        if (is.character(df[,which(names(df) == vars)]) == TRUE) {
          # print("true")
          charvars <- c(charvars, vars)
        }
      }
      return(charvars)
    }
    
    binaryvarnames <- which_character(Dx_combined)
    Dx_combined <- as.data.frame(Dx_combined)
    patient_binary <- Dx_combined[, which(colnames(Dx_combined) %in% c("patientunitstayid", binaryvarnames))]
    
    dummies <- fastDummies::dummy_cols(patient_binary, select_columns = c(binaryvarnames), remove_selected_columns = T, ignore_na = T)
    # dummies[is.na(dummies)] <- 0

    Dx_combined <- Dx_combined %>% dplyr::select(-System)
    Dx_combined <- merge(Dx_combined, dummies, by = "patientunitstayid")
    
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
    if (make_binary) {
      binaryvarnames <- which_character(patient_table)
      patient_table <- as.data.frame(patient_table)
      patient_binary <- patient_table[, which(colnames(patient_table) %in% c("patientunitstayid", binaryvarnames))]
      
      dummies <- fastDummies::dummy_cols(patient_binary, select_columns = c(binaryvarnames), remove_selected_columns = T, ignore_na = T)
      
      patient_table <- patient_table %>% dplyr::select(-all_of(binaryvarnames))
      patient_table <- merge(patient_table, dummies, by = "patientunitstayid", all = T)
    }
    
    return(patient_table)
    
  } else if (label_boolean == TRUE) {
    
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
    if (make_binary) {
      binaryvarnames <- which_character(patient_table)
      patient_table <- as.data.frame(patient_table)
      patient_binary <- patient_table[, which(colnames(patient_table) %in% c("patientunitstayid", binaryvarnames))]
      
      dummies <- fastDummies::dummy_cols(patient_binary, select_columns = c(binaryvarnames), remove_selected_columns = T, ignore_na = T)
      
      patient_table <- patient_table %>% dplyr::select(-all_of(binaryvarnames))
      patient_table <- merge(patient_table, dummies, by = "patientunitstayid", all = T)
    }
    
    return(patient_table)
  }
  
}

