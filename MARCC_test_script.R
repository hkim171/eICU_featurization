### MARCC TEST CODE

#MY DIR - for local machine
code_folder_dir <- "/storage/eICU/eICU_feature_extract/"
code_dir <- paste0(code_folder_dir, "/eICU_featurization")


#MARCC DIR
# code_folder_dir <- "/home-2/hkim171@jhu.edu/scratch/"
# code_dir <- paste0(code_folder_dir, "/eICU_featurization")

source(paste0(code_dir, "/model_guide.R"))

build_prototype(code_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/", 
                save_dir = "/storage/eICU/eICU_feature_extract/",
                
                merge_identifier = "patientunitstayid",
                feature_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_feature_space.csv",
                label_dir = "/storage/eICU/eICU_feature_extract/eICU_featurization/test_data/test_label.csv",
                outcomes = c("Expired", "Alive"),
                experiment_name = "test",
                Already_trained = TRUE
)
