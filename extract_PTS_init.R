#EXTRACT initial statistics features applicable to all PTS data. 


#author:Han
#date: 1 23 2020
#updated: 4 30 20 -> commented

## arguments
#@args table_df : name of table data is extracted from
#@args patientunitstayid : this feeds in 1 patient unit at a time since lower and upper times may be different per pid. 
#@args lower_hr : lower time bound. In eICU database, this is represented as the offset from the start of ICU stay. 
#@args upper_hr : upper time bound. 
#@args binwidth : bin width to extract initial statistics
#@args ts_length_hour : length of the time series data. This is mainly used if there is a uniform number of bins per PID. 
#@args per_pid_boolean : Needs to be set to TRUE .

## output
# init_stats : c("hour_bin", "patientunitstayid", "mean", "min", "max", "sd", "range", "delta", "data_count")
# dataframe with X number of bin hours per patientunitstayid.

extract_PTS_init_stats <- function(patientunitstayids, table_df, pts_name, lower_hr, upper_hr, binwidth, ts_length_hour, per_pid_boolean) {
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
  
  if (missing(ts_length_hour)) {
    ts_length_hour == ceiling(((upper_hr - lower_hr) * 60))
  }
  
  if (per_pid_boolean == TRUE) {
    patientunitstayids <- c(patientunitstayids)
  }
  
  
  
  # filter for pids and time frame
  table_df <- table_df[which(table_df$patientunitstayid %in% patientunitstayids), ] %>% droplevels(.) %>% distinct(.)
  table_df <- table_df[which(table_df$offset >= (lower_hr) & table_df$offset <= (upper_hr)), ] %>% droplevels(.) %>% distinct(.)
  
  features <- as.data.frame(patientunitstayids)
  colnames(features) <- c("patientunitstayid")
  id <- "patientunitstayid"
  
  #extracting a dataframe of statistics per binwidth. 
  
  #computing number of bins to iterate through using ts_length_hour and binwidth
  bins <- ceiling(ts_length_hour / binwidth)
  
  #empty dataframe to append to
  init_stats <- rbind()
  
  for (i in 1:bins) {
    #for each bin, iterate through and get each binwidth time frame. 
    templower <- lower_hr + ((i - 1) * (binwidth * 60))
    tempupper <- lower_hr + ((i) * (binwidth * 60))
    tempdf <- table_df[which(table_df$offset >= (templower) & table_df$offset <= (tempupper)), ]
    tempdf <- tempdf[complete.cases(tempdf), ]
    #calculate delta
    
    if (nrow(tempdf) != 0) {
      delta <- (tempdf$value[which.max(tempdf$offset)] - tempdf$value[which.min(tempdf$offset)]) /
        (tempdf$offset[which.max(tempdf$offset)] - tempdf$offset[which.min(tempdf$offset)])
    } else {
      delta <- NA
    }

    tempdf <- tempdf$value

    tempstats <- c((i * 2), patientunitstayids, 
                   mean(tempdf, na.rm = T), 
                   min(tempdf, na.rm = T), 
                   max(tempdf, na.rm = T), 
                   sd(tempdf, na.rm = T), 
                   (max(tempdf, na.rm = T) - min(tempdf, na.rm = T)), 
                   delta, 
                   length(tempdf))
    
    tempstats[is.nan(tempstats)] <- NA
    tempstats[is.infinite(tempstats)] <- NA
    
    init_stats <- rbind(init_stats, tempstats)
    
  }
  
  colnames(init_stats) <- c("hour_bin", "patientunitstayid", "mean", "min", "max", "sd", "range", "delta", "data_count")
  
  
  return(init_stats)
}
  
