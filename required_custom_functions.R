#assortment of functions

#function created to simplify the loading of eICU data based on table directory. 
#This will help to use this code on other machines. 
read.eicu <- function(directory, tablename, extension) {
  #checks for two extensions: variations of csv and rds
  #extension should represent the case of the extension of the file
  #to be loaded 
  
  if (grepl(".", extension, fixed=TRUE)) {
    extension <- gsub("[[:punct:]]", "", extension)
  }
  
  if(extension == "csv" | extension == "CSV" | extension == "Csv") {
    
    return(data.table::fread(paste0(directory, "/",tablename, ".", extension)))
    
  } else if (extension == "rds" | extension == "RDS" | extension == "Rds") {
    
    return(readRDS(paste0(directory, "/",tablename, ".", extension)))
  } else {
    stop("check extension")
  }
  
}



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
create_binary <- function(df, id_name) {
  
  if(missing(id_name)){
    id_name = "patientunitstayid"
  }
  
  binaryvarnames <- which_character(df)
  
  if (id_name %in% binaryvarnames) {
    
    binaryvarnames <- binaryvarnames[-which(binaryvarnames %in% id_name)]
    
  }
  
  df <- as.data.frame(df)
  # df <- df %>% select(-id_name)
  patient_binary <- df[, which(colnames(df) %in% c(id_name, binaryvarnames))]
  
  dummies <- fastDummies::dummy_cols(patient_binary, select_columns = c(binaryvarnames), remove_selected_columns = T, ignore_na = T)
  # dummies[is.na(dummies)] <- 0
  
  df <- df %>% dplyr::select(-all_of(binaryvarnames))
  df <- merge(df, dummies, by = id_name)
  
  return(df)
}

#replaces infinites in a dataframe with replace_with
replace_infinite <- function(df, replace_with) {
  df <- do.call(data.frame,lapply(df, function(x) replace(x, is.infinite(x),replace_with)))
  return(df)
}

#replaces NAN in a dataframe with replace_with
replace_nan <- function(df, replace_with) {
  df <- do.call(data.frame,lapply(df, function(x) replace(x, is.nan(x),replace_with)))
  return(df)
}

#replaces NAs in a dataframe with replace_with
replace_NA <- function(df, replace_with) { 
  df[is.na(df)] = replace_with
  return(df)
}

#Counts number of NAs in a dataframe
check_na <- function(df) {
  table(is.na(df))
}

#counts number of Inf and -Inf in a dataframe
check_infinite <- function(df){
  yes <- 0
  no <- 0
  for (i in 1:ncol(df)) {
    tempyes <- sum(is.infinite(df[,i]))
    yes <- yes + tempyes
    no <- no + (length(df[,1]) - yes)
  }
  df <- data.frame(no, yes, row.names = "")
  colnames("FALSE", "TRUE")
  print(as.data.frame(df))
}

#polvar function: 
#Measures the probability of obtaining a sequence of consecutive ones or zeros.
#---INPUTS:
#x, the input time series
#d, the symbolic coding (amplitude) difference,
#D, the word length (classically words of length 6).
#---OUPUT:
#p - probability of obtaining a sequence of consecutive ones/zeros

polvar <- function(x, d, D) {
  dx = abs(diff(x))
  N = length(dx)
  
  xsym <- (dx >= d) * 1
  zseq <- rep(0, D)
  oseq <- rep(1, D)
  
  #search for D consecutive zeros
  i = 1;
  pc = 0;
  
  while (i <= (N - D)) {
    xseq <- xsym[i:(i + D - 1)]
    
    if ((sum(xseq == zseq) == D) | (sum(xseq == oseq) == D) ) {
      pc = pc + 1
      i = i + D
    } else {
      i = i + 1;
    }
  }
  p = pc / N
  
  return(p)
}

# EN_PermEn     Permutation Entropy of a time series.
#
# "Permutation Entropy: A Natural Complexity Measure for Time Series"
# C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
#
#---INPUTS:
# y, the input time series
# m, the embedding dimension (or order of the permutation entropy)
# tau, the time-delay for the embedding
#
#---OUTPUT:
# Outputs the permutation entropy and normalized version computed according to
# different implementations
permEN <- function(x, m) {
  require(stats)
  require(statcomp)
  # tauseq <- seq(from = 1, to = tau, by = 1)
  # xtau <- x[-tauseq]
  xtau = x
  
  opd = ordinal_pattern_distribution(x = xtau, ndemb = m)
  
  if (sum(opd) == 0) {
    NA
  } else {
    return(permutation_entropy(opd))
  }
  
}

db <- function(x, filter_name) {
  require(wavelets)
  
  if (length(x) < 4){
    return(NA)
  } else {
    d <- dwt(x, filter = filter_name)
    coeff <- sapply(d@W, mean)
    return(coeff)
  }
  
}


## Remove binary features 
remove_binary <- function(df) {
  binary <- as.data.frame(apply(df, 2, function(x) length(unique(x))))
  binary <- tibble::rownames_to_column(binary, "feature")
  binary <- binary[which(binary$`apply(df, 2, function(x) length(unique(x)))` <= 3), ]
  binary <- binary$feature
  df_nonbinary <-df[, -which(names(df) %in% binary)]
  return(df_nonbinary)
}

# makes binary features into factors in a dataframe
binary_to_factor <- function(df) {
  binary <- as.data.frame(apply(df, 2, function(x) length(unique(x))))
  binary <- tibble::rownames_to_column(binary, "feature")
  binary <- which(binary$`apply(df, 2, function(x) length(unique(x)))` <= 3)
  
  df[, binary] <-  lapply(df[, binary] , as.factor)  
  return(df)
}

# gives the column names of binary features
which_binary <- function(df, name_boolean) {
  binary <- as.data.frame(apply(df, 2, function(x) length(unique(x))))
  binary <- tibble::rownames_to_column(binary, "feature")
  binary <- which(binary$`apply(df, 2, function(x) length(unique(x)))` <= 3)
  if (name_boolean == TRUE) {
    binary <- colnames(df)[binary]
    return(binary)
  } else {
    return(binary)
  }
}

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

#gives column names of features that are continuous
which_continous <- function(df) {
  source("~/useful_functions.R")
  contvars <- c()
  varnames <- colnames(df)
  
  binary <- which_binary(df, name_boolean = T)
  char <- which_character(df)
  
  varnames <- varnames[-which(varnames %in% c(binary, char))]
  return(varnames)
  
}

#converts character features to a factor
char_to_factor <- function(df) {
  charvars <- c()
  varnames <- colnames(df)
  for (vars in varnames) {
    if (is.character(df[,which(names(df) == vars)]) == TRUE) {
      charvars <- c(charvars, vars)
    }
  }
  
  df[, charvars] <-  lapply(df[, charvars] , as.factor)  
  return(df)
}

#converts an RDS files into CSV given file path. 
writeRDS_to_CSV <- function(path_filename, ...) {
  
  new_filename <- str_replace(path_filename, ".Rds", ".csv")
  write.csv(readRDS(path_filename), file = new_filename, row.names = F, ...)
}

unique_pid <- function(df) {
  return(length(unique(df$patientunitstayid)))
}

#shannon_entropy
shannon_entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

energy <- function(target) {
  target <- target[target > 0]
  sum(target^2)/length(target)
}

#isoreg fit
fit.isoreg <- function(iso, x0) 
{
  o = iso$o
  if (is.null(o)) 
    o = 1:length(x)
  x = iso$x[o]
  y = iso$yf
  ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
  min.x <- min(x)
  max.x <- max(x)
  adjusted.knots <- iso$iKnots[c(1, which(iso$yf[iso$iKnots] > 0))]
  fits = sapply(seq(along = x0), function(i) {
    j = ind[i]
    
    # Handles the case where unseen data is outside range of the training data
    if (is.na(j)) {
      if (x0[i] > max.x) j <- length(x)
      else if (x0[i] < min.x) j <- 1
    }
    
    # Find the upper and lower parts of the step
    upper.step.n <- min(which(adjusted.knots > j))
    upper.step <- adjusted.knots[upper.step.n]
    lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )
    
    # Pefrom a liner interpolation between the start and end of the step
    denom <- x[upper.step] - x[lower.step] 
    denom <- ifelse(denom == 0, 1, denom)
    val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
    
    # Ensure we bound the probabilities to [0, 1]
    val <- ifelse(val > 1, max.x, val)
    val <- ifelse(val < 0, min.x, val)
    val <- ifelse(is.na(val), max.x, val) # Bit of a hack, NA when at right extreme of distribution
    val
  })
  fits
}

load_packages <- function(package_list) {
  #   Checks if given library is installed. If it is not, installs and loads it.
  # 	If it is, loads it.
  # 
  # 	Args:
  # 	library_list: A list of vectors containing strings of the required
  # 					packages
  # 	Returns:
  # 	Null
  
  for (pkg in package_list) {
    if (!require(pkg, character.only=TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}