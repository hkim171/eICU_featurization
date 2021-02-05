#assortment of functions

#finds which columns of a dataframe are character vectors
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

#converts character vectors into one hot encoded columns and removes original column.
create_binary <- function(df) {
  binaryvarnames <- which_character(df)
  df <- as.data.frame(df)
  patient_binary <- df[, which(colnames(df) %in% c("patientunitstayid", binaryvarnames))]
  
  dummies <- fastDummies::dummy_cols(patient_binary, select_columns = c(binaryvarnames), remove_selected_columns = T, ignore_na = T)
  # dummies[is.na(dummies)] <- 0
  
  df <- df %>% dplyr::select(-all_of(binaryvarnames))
  df <- merge(df, dummies, by = "patientunitstayid")
  
  return(df)
}
