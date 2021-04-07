

extract_last_GCS <- function(patientunitstayids, table_df) {
  source("~/useful_functions.R")
  require(tidyverse)
  require(plotly)
  library(caret)
  library(doParallel)

  
  #filter for patients and time frame
  table_df <- table_df[which(table_df$patientunitstayid %in% patientunitstayids), ]

  variables <- c("Verbal", "Motor", "Eyes", "GCS Total")
  
  pids <- as.data.frame(patientunitstayids)
  colnames(pids) <- "patientunitstayid"
  
  for (variable in variables) {
    tempdf <- table_df %>% select(patientunitstayid, nursingchartoffset, nursingchartcelltypevalname, 
                                  nursingchartvalue)
    
    tempdf <- tempdf[which(tempdf$nursingchartcelltypevalname %in% variable),] %>% droplevels(.)
    tempdf <- tempdf[complete.cases(tempdf),]
    
    if (variable == "GCS Total") {
      tempdf <- tempdf[which(tempdf$nursingchartvalue != "Unable to score due to medication"),]
    }
    
    tempdf$nursingchartvalue <-as.numeric(as.character(tempdf$nursingchartvalue))
    tempdf <- tempdf %>% select(patientunitstayid, nursingchartoffset, nursingchartvalue)
    
    last_val <- tempdf %>% group_by(patientunitstayid) %>% slice(which.max(nursingchartoffset))
    
    last_val <- last_val %>% select(patientunitstayid, nursingchartvalue)
    
    if (variable == "GCS Total") {
      variable <- "GCS_Total"
    }
    
    colnames(last_val) <- c("patientunitstayid", paste0("last_ICU_",variable))

    pids <- merge(pids, last_val, by = "patientunitstayid", all = T)
    
  }
  
  
  return(pids)
  
}
