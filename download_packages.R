#download packages in R

## INSTRUCTIONS
#make sure use following commands in terminal within MARCC.
#ml R/3.6.1
#ml gcc/5.5

#run this script with Rscript download_packages.R in MARCC terminal


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

package_list <- c("tidyverse", "data.table", "ranger", "randomForest", "caret", "MLmetrics", "Metrics", 
                  "EvaluationMeasures", "Rmisc", "xgboost", "glmnet", "randomForestSRC")

load_packages(package_list)