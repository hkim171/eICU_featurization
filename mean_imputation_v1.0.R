#### Usage: extract features from patient table
# Author: Jocelyn
# Updated and Formatted for distribution: MAR 2021
# V1.0
# Fixes:
#-------------------------------------------------------------------------------------------------------------------------------------#

#### input parameters: ####
# data_file:                 data.table/data.frame - List of patientunitstayids as numerics
#
#-------------------------------------------------------------------------------------------------------------------------------------#

#### files required: ####
# data_file

#### returns: ####
# data_imputed               data.table - data_file with NA values imputed with mean of column OR 0 if column is binary 0/1 variable

#### Code example is provide below. select lines 20-58 then click ctrl+shift+c to uncomment. ####
source("~/scratch/eICU_featurization/mean_imputation_v1.0.R")
data_table <- mean_imputation("data_file.csv")

#### function ####

mean_imputation <- function(data_file) {
    
    if (!file.exists(data_file)) {
        stop("Input data file does not exist")
    }
    data <- read.csv(data_file)
    
    count_na <- 0
    num_pids <- length(data$patientunitstayid)
    data_imputed <- data
    for (i in 1:ncol(data_imputed)) {
        count_na = sum(is.na(data_imputed[i]))
        count <- 0
        sum <- 0
        if (count_na > 0) {
            for (j in 1:length(data_imputed[,i])) {
                if (!is.na(data_imputed[j, i])) {
                    sum = sum + data_imputed[j, i]
                    if (data_imputed[j, i] == 0 || data_imputed[j, i] == 1) {
                        count = count + 1
                    }
                }
            }
            avg <- sum / (num_pids - count_na);
            print(avg)
            for (j in 1:length(data_imputed[,i])) {
                if (is.na(data_imputed[j, i])) {
                    if (count == num_pids - count_na) {
                        data_imputed[j, i] = 0
                    } else {
                        data_imputed[j, i] = avg
                    }
                }
            }
        }
    }
    return(data_imputed)
    
}

