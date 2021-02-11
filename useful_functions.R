#helpful functions that gets used

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

find_max_bin_offset <- function(carvedilol_dataframe, hr_dataframe, threshold) {
  new_carvedilol_dataframe <- carvedilol_dataframe[1,]
  new_carvedilol_dataframe$ADR_onset <- NA
  new_carvedilol_dataframe <- new_carvedilol_dataframe[-1,]
  
  for (i in 1:nrow(carvedilol_dataframe)) {
    tempdf <- carvedilol_dataframe[i, ]
    temp_start <- carvedilol_dataframe$drugstartoffset[i]
    # print(carvedilol_dataframe$drugstartoffset[i])
    # print(carvedilol_dataframe$drugstartoffset[i] + (12 * 60))
    temp_hr <- hr_carvedilol[which(hr_carvedilol$patientunitstayid == carvedilol_dataframe$patientunitstayid[i] & 
                                     hr_carvedilol$offset > carvedilol_dataframe$drugstartoffset[i] & 
                                     hr_carvedilol$offset <= (carvedilol_dataframe$drugstartoffset[i] + (12*60)) &
                                     hr_carvedilol$value < threshold), ]
    
    count_per_bin <- c()
    start_offsets <- c()
    for (j in 1:12) {
      start_offsets[j] <- temp_start + (j-1) * 60
      count_per_bin[j] <- sum(temp_hr$offset > (temp_start + (j-1) * 60) & temp_hr$offset <= (temp_start + (60 * j)))
    }
    
    max_count_offset <- start_offsets[which.max(count_per_bin)]
    tempdf$ADR_onset <- max_count_offset
    
    new_carvedilol_dataframe <- rbind(new_carvedilol_dataframe, tempdf)
  }
  
  return(new_carvedilol_dataframe)
}

#gets the histogram bins and counts of a dataframe object. 
get_hist <- function(p) {
  d <- ggplot_build(p)$data[[1]]
  data.frame(x = d$x, xmin = d$xmin, xmax = d$xmax, y = d$y)
}

# builds a simple histogram with ggplot2 given parameters
build_gghist <- function(df, colname, binsize, lower, upper) {
  library(tidyverse)
  df <- df %>% select(colname)
  
  m <- mean(df[,1], na.rm = T)
  s <- sd(df[,1], na.rm = T)
  
  s1 <- c(-1,1)*s + m
  s2 <- c(-1,1)*s*2 + m
  s3 <- c(-1,1)*s*3 + m
  s4 <- c(-1,1)*s*4 + m
  s5 <- c(-1,1)*s*5 + m
  s6 <- c(-1,1)*s*6 + m
  
  lower <- -1*s*8 + m
  upper <- 1*s*8 + m
  
  ggplot(df) + 
    geom_histogram(aes(x = df[,1]), binwidth = binsize, color = "darkgray", fill = "lightgray", alpha = .5) +
    geom_vline(xintercept = m, linetype = 'dashed', color = 'red', alpha = 0.9) +
    geom_vline(xintercept = s1, linetype = 'dashed', color = 'blue', alpha = 0.7) +
    geom_vline(xintercept = s2, linetype = 'dashed', color = 'blue', alpha = 0.6) +
    geom_vline(xintercept = s3, linetype = 'dashed', color = 'blue', alpha = 0.5) +
    geom_vline(xintercept = s4, linetype = 'dashed', color = 'blue', alpha = 0.4) +
    geom_vline(xintercept = s5, linetype = 'dashed', color = 'blue', alpha = 0.3) +
    geom_vline(xintercept = s6, linetype = 'dashed', color = 'blue', alpha = 0.2) +
    xlim(lower, upper) +
    xlab(colname) +
    ggtitle(paste0(colname, " histogram with mean(r), sd(blue)"))
}

#plots quick PCA chart
pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}

#shows basic information about a selected variable of a dataframe. 
columnContent <- function(x, select_variable) {
  if (missing(select_variable)) {
    print(table(x))
    print(length(unique(x)))
    hist(x, breaks = length(x)/2)
  } else if (select_variable == "length") {
    print(length(unique(x)))
  } else if (select_variable == "table") {
    print(table(x))
  } else if (select_variable == "histogram") {
    hist(x, breaks = length(x)/2)
  } 
  
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


read.MIMIC <- function(table_name) {
  return(read.csv(paste0("/storage/MIMIC/", table_name, ".csv")))
}

job_complete <- function () {
  library(emayili)
  library(dplyr)
  
  email <- emayili::envelope() %>%
    emayili::from("hanbiehn@gmail.com") %>%
    emayili::to("hkim171@jhu.edu") %>%
    emayili::subject("Job Completed!") %>% 
    emayili::subject("go check the results!")
  
  smtp <- server(host = "smtp.gmail.com",
                 port = 465,
                 username = "hanbiehn@gmail.com",
                 password = "rbyiqdwvjvvmmtny")
  # Finally send the message.
  
  smtp(email, verbose = TRUE)
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