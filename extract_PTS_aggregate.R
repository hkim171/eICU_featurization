#Aggregate simple stats by time frame into aggregate feature. 
#1. Simple statistics_only

#author:Han
#date: 1 23 2020
#update: 4 30 20 : comments added. 

#@args : init_stat_df : a full list of statistics extracted out of extract_pts_init. 
#@identifier : used to name features. 

extract_PTS_aggregate <- function(init_stat_df, pids, identifier) {
  source("~/useful_functions.R")
  require(tidyverse)
  require(plotly)
  library(caret)
  library(doParallel)
  
  init_stat_df <- as.data.frame(init_stat_df)
  # pids <- unique(init_stat_df$patientunitstayid)
  pids <- as.data.frame(pids)
  colnames(pids) <- "patientunitstayid"

  aggregate_stats <- rbind()
  #iterates through each pids within the init_df to gather aggregate statistics for each init statistics. 
  for (i in 1:nrow(pids)) {
    print(i)
    pid <- pids$patientunitstayid[i]
    tempdf <- init_stat_df[which(init_stat_df$patientunitstayid == pid), ]
    tempdf <- tempdf[complete.cases(tempdf), ]
    
    varnames <- c("mean", "min", "max", "sd", "range", "delta")
    new_varnames <- c()
    agg_stats <- c()
    #iterate through all init statistics. 
    for (j in 3:8) {
      vals <- tempdf[,j]
      
      agg_names <- c("mean", "min", "max", "sd", "range")
      tempvarnames <- paste0(agg_names, "_", colnames(tempdf)[j], "_", identifier)
      new_varnames <- c(new_varnames, tempvarnames)
      
      if (is_empty(vals)) {
        agg_stats <- c(agg_stats, 
                       NA, 
                        NA, 
                       NA, 
                       NA, 
                       NA)
        
        agg_stats[is.nan(agg_stats)] <- NA
        agg_stats[is.infinite(agg_stats)] <- NA
        
        if (j == 8) {
          agg_stats <- c(agg_stats, NA)
        }
        
      } else {
        agg_stats <- c(agg_stats, 
                       mean(vals, na.rm = T), 
                       min(vals, na.rm = T), 
                       max(vals, na.rm = T), 
                       sd(vals, na.rm = T), 
                       (max(vals, na.rm = T) - min(vals, na.rm = T)))
        
        agg_stats[is.nan(agg_stats)] <- NA
        agg_stats[is.infinite(agg_stats)] <- NA
        
        if (j == 8) {
          agg_stats <- c(agg_stats, sum(diff(sign(vals)) != 0))
        }
      }
      
      
    }
    
    # #get sign change feature for delta. 
    # delta <- tempdf$delta
    # delta_sign_change <- sum(diff(sign(delta)) != 0)
    # 
    
    agg_stats <- c(pid, agg_stats)
    
    aggregate_stats <- rbind(aggregate_stats, agg_stats)
  
  }
  
  colnames(aggregate_stats) <- c("patientunitstayid", new_varnames, paste0("delta_sign_change_", identifier))

  return(aggregate_stats)
}
